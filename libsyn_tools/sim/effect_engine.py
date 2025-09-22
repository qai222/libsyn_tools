from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import List, Union, Iterable as It, Optional
from uuid import uuid4

import pandas as pd
import simpy
from loguru import logger
from pyshacl import validate
from rdflib import Graph, Literal, URIRef, Namespace, ConjunctiveGraph
from rdflib.namespace import XSD, SH
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph.physical_entities import BaseClass, MaterialContainer
from libsyn_tools.sim.operation import Operation, UnitaryEdit, UnitaryEditType, FilterStoreRegistry, get_runtime_state
from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE, _needs_runtime_tracking
from .effect_shacl import SHACLViolationRecord, _iter_validation_results, _first

LIB_SYN = Namespace("https://libsyn-sim/kg/")


class SHACLValidationError(RuntimeError):
    """
    Raised when the knowledge graph violates one or more user-supplied
    SHACL shapes *after* an Operation has been applied (policy violation).
    """
    pass


class EffectEngine:
    """
    Apply `UnitaryEdit`s for one operation as a **micro-commit**:
    - run *mechanical* hard checks first (abort if fail; no edits applied)
    - apply edits in-order (no yield; single SimPy microstep)
    - run SHACL over (data ⊎ overlay); record soft violations

    Note: We *do not* roll back SHACL violations. They are recorded and
    the simulation proceeds (unless `raise_shacl=True`).
    """

    def __init__(
            self,
            *,
            shapes_graph: ConjunctiveGraph | None = None,
            raise_shacl: bool = False,
            inference: str = "owlrl"
    ):
        """
        Parameters
        ----------
        shapes_graph
            RDF graph containing user-defined SHACL shapes.
        raise_shacl
            If True, raise on SHACL non-conformance (after recording).
        inference
            Inference mode for pySHACL (e.g., "none", "rdfs", "owlrl").
        """
        self.shapes_graph = shapes_graph
        self.inference = inference
        self._shacl_violations: list[SHACLViolationRecord] = []
        self.raise_shacl = raise_shacl

    # ------------------------------------------------------------------ #
    # Overlay construction (ephemeral, unioned with the data graph)
    # ------------------------------------------------------------------ #
    def _build_overlay_graph(self) -> Graph:
        """
        Construct an *ephemeral* rdflib.Graph that contains **derived**
        triples needed for SHACL validation. Not persisted.

        Currently: adds LIB_SYN.currentVolume for each present container.
        """
        g = Graph()
        for c in MaterialContainer.object_lookup.values():
            if c.is_present != {True}:  # skip annihilated objects
                continue
            vol = c.directly_contained_pom_volume
            g.add(
                (
                    URIRef(c.instance_iri),
                    LIB_SYN.currentVolume,
                    Literal(vol, datatype=XSD.double),
                )
            )
        return g

    # ------------------------------------------------------------------ #
    # SHACL helpers
    # ------------------------------------------------------------------ #
    def _collect_shacl_violations(
            self,
            *,
            env_now: float,
            operation_id: str,
            shacl_report_graph: Graph,
            batch_id: Optional[str],
            edit_fingerprints: Optional[list[str]],
            seed: Optional[int]
    ) -> list[SHACLViolationRecord]:
        """
        Transform the SHACL results graph into uniform violation records.
        """
        records: list[SHACLViolationRecord] = []
        for vr in _iter_validation_results(shacl_report_graph):
            rec = SHACLViolationRecord(
                sim_time=env_now,
                operation_id=operation_id,
                origin="SHACL",
                severity="soft",
                disposition="committed",
                batch_id=batch_id,
                edit_fingerprints=edit_fingerprints,
                seed=seed,
                shape_iri=_first(shacl_report_graph, vr, SH.sourceShape),
                focus_iri=_first(shacl_report_graph, vr, SH.focusNode),
                message=_first(shacl_report_graph, vr, SH.resultMessage),
                report_graph_ttl=shacl_report_graph.serialize(format="turtle"),
            )
            records.append(rec)
        return records

    def _run_shacl_validation(
            self,
            *,
            env: simpy.Environment,
            operation_id: str,
            batch_id: Optional[str],
            edit_fingerprints: Optional[list[str]],
            seed: Optional[int]
    ) -> None:
        """
        Validate the union graph (data + overlay). Record (and optionally raise)
        SHACL violations as *soft* outcomes.
        """
        if self.shapes_graph is None:
            return

        data_graph: Graph = KnowledgeGraph.graph()
        overlay_graph: Graph = self._build_overlay_graph()

        # merged view for validation
        union_graph = ConjunctiveGraph()
        for t in data_graph.triples((None, None, None)):
            union_graph.add(t)
        for t in overlay_graph.triples((None, None, None)):
            union_graph.add(t)

        conforms, shacl_report_graph, _ = validate(
            union_graph,
            shacl_graph=self.shapes_graph,
            ont_graph=self.shapes_graph,
            inference=self.inference,
            advanced=True,
            debug=False,
        )
        logger.debug(f"shapes conform: {conforms}")
        if not conforms:
            logger.error(shacl_report_graph.serialize(format="turtle"))
            self._shacl_violations.extend(
                self._collect_shacl_violations(
                    env_now=env.now,
                    operation_id=operation_id,
                    shacl_report_graph=shacl_report_graph,
                    batch_id=batch_id,
                    edit_fingerprints=edit_fingerprints,
                    seed=seed,
                )
            )
            if self.raise_shacl:
                raise SHACLValidationError(shacl_report_graph)

    def write_shacl_csv(
            self,
            filepath: Union[str, Path],
            include_ttl: bool = False,
            **to_csv_kwargs,
    ) -> None:
        """
        Export all captured SHACL violations to a CSV via pandas.

        Parameters
        ----------
        filepath : str | pathlib.Path
            Where the CSV will be written.
        include_ttl : bool, default False
            If True, include the full Turtle string in a column `report_graph_ttl`.
            This can make the file very large; off by default.
        **to_csv_kwargs : dict
            Extra keyword arguments forwarded to `DataFrame.to_csv`.
            (e.g. sep=';', encoding='utf-8', etc.)
        """
        records_as_dicts = [
            rec.model_dump() if include_ttl else rec.model_dump(exclude={"report_graph_ttl"})
            for rec in self._shacl_violations
        ]
        df = pd.DataFrame.from_records(records_as_dicts)
        df.to_csv(Path(filepath), index=False, **to_csv_kwargs)

    # ------------------------------------------------------------------ #
    # Edit application pipeline
    # ------------------------------------------------------------------ #
    @staticmethod
    def _fingerprint_edit(edit: UnitaryEdit) -> str:
        """
        Produce a stable, human-readable fingerprint for an edit.
        Used in violation records for replay/debug.
        """
        dv = edit.data_value
        dv_repr = str(dv) if dv is None or isinstance(dv, (int, float, str, bool)) else f"<{type(dv).__name__}>"
        return "|".join(
            [
                edit.type.value,
                edit.instance_1_iri or "",
                edit.property_iri or "",
                edit.instance_2_iri or "",
                dv_repr,
            ]
        )

    def _precheck_mechanical(
            self,
            *,
            edits: It[UnitaryEdit],
            locked_iris: Optional[It[str]] = None
    ) -> None:
        """
        Mechanical (hard) checks that must pass *before* any edits apply.

        Minimal, non-policy safety:
        - duplicate CREATE on same IRI in batch
        - referenced IRIs must exist or be created by this batch
          (treat get_object_from_lookup(...) returning None as "not exists")
        - CREATE on an already-present object is an error
        - coverage (optional): for object-property writes, both endpoints
          must be either (a) created in batch or (b) among `locked_iris`
        """
        locked = set(locked_iris or [])
        edits_list = list(edits)

        creates = [e.instance_1_iri for e in edits_list if e.type is UnitaryEditType.CREATE]
        if len(creates) != len(set(creates)):
            raise RuntimeError(f"Mechanical check failed: duplicate CREATE in batch: {creates}")

        creates_set = set(creates)

        def _exists(iri: Optional[str]) -> bool:
            if not iri:
                return False
            try:
                obj = KnowledgeGraph.get_object_from_lookup(iri)
            except Exception:
                return False
            return obj is not None  # IMPORTANT: treat None as "not exists"

        for e in edits_list:
            t = e.type

            if t is UnitaryEditType.CREATE:
                subj = KnowledgeGraph.get_object_from_lookup(e.instance_1_iri)
                # If subj is None here it's okay (node exists in registry, presence toggled by CREATE)
                # We only forbid CREATE on an already-present object:
                if subj is not None and getattr(subj, "is_present", {False}) == {True}:
                    raise RuntimeError(f"Mechanical check failed: CREATE on present object {e.instance_1_iri}")

            elif t in (
                    UnitaryEditType.ADD_DATA_PROPERTY,
                    UnitaryEditType.CHANGE_DATA_PROPERTY,
                    UnitaryEditType.ANNIHILATE,
            ):
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling subject {e.instance_1_iri}")

            elif t in (
                    UnitaryEditType.ADD_OBJECT_PROPERTY,
                    UnitaryEditType.REMOVE_OBJECT_PROPERTY,
            ):
                # subj & obj must exist or be created
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling subject {e.instance_1_iri}")
                if not (_exists(e.instance_2_iri) or e.instance_2_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling object {e.instance_2_iri}")

                # coverage: endpoints must be locked or created
                if locked:
                    for iri in (e.instance_1_iri, e.instance_2_iri):
                        if iri not in locked and iri not in creates_set:
                            raise RuntimeError(
                                "Mechanical check failed: write coverage requires lock or create; "
                                f"got edit on {iri!r} not in locked set {locked}"
                            )

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
        """
        Return a **defensive copy** of the staged edits that `action`
        prepared in its pre-act phase.
        """
        return action.operation_effects.copy()

    def apply(
            self,
            edits: Iterable[UnitaryEdit],
            env: simpy.Environment,
            *,
            operation_id: str,
            locked_iris: Optional[Iterable[str]] = None,
            seed: Optional[int] = None,
    ) -> None:
        """
        Apply `edits` in a single SimPy microstep with **mechanical pre-checks**
        and SHACL auditing.

        Parameters
        ----------
        edits : list of UnitaryEdit
            The staged edits for a single operation instance.
        env : simpy.Environment
            Current simulation environment (for timestamps).
        operation_id : str
            The identifier of the operation whose edits these are.
        locked_iris : Optional[Iterable[str]]
            IRIs locked by the operation (participants). Used for the
            *coverage* safety check for object-property writes.
        seed : Optional[int]
            Run seed (for uniform violation records).
        """
        edits = list(edits)  # we iterate multiple times
        batch_id = str(uuid4())
        fingerprints = [self._fingerprint_edit(e) for e in edits]

        # 1) mechanical hard checks (abort before commit on failure)
        self._precheck_mechanical(edits=edits, locked_iris=list(locked_iris or []))

        # 2) apply each edit (in-order). No yields ⇒ single microstep.
        for edit in edits:
            subj = KnowledgeGraph.get_object_from_lookup(edit.instance_1_iri)
            obj2 = (
                KnowledgeGraph.get_object_from_lookup(edit.instance_2_iri)
                if edit.type in (UnitaryEditType.ADD_OBJECT_PROPERTY,
                                 UnitaryEditType.REMOVE_OBJECT_PROPERTY) and edit.instance_2_iri
                else None
            )

            # ensure runtime artefacts for LabObjects
            self._register_if_new(subj, env)
            if edit.type is not UnitaryEditType.CREATE and obj2:
                self._register_if_new(obj2, env)

            logger.debug(f"Applying edit: {edit.type} – {subj.__class__.__name__}={edit.instance_1_iri}")
            edit.apply()

            # resource/FilterStore normalisation + history
            self._log_and_sync(subj, edit, env)
            self._log_and_sync(obj2, edit, env)

            if edit.type is UnitaryEditType.CREATE:
                # now present ⇒ ensure resource & pool
                self._register_if_new(subj, env)
            elif edit.type is UnitaryEditType.ANNIHILATE:
                # remove resource + filter entry only after logging
                self._unregister_object(subj)

        # 3) SHACL audit (soft violations; recorded and optionally raised)
        self._run_shacl_validation(
            env=env,
            operation_id=operation_id,
            batch_id=batch_id,
            edit_fingerprints=fingerprints,
            seed=seed,
        )

    # ------------------------------------------------------------------ #
    # Runtime artefact helpers
    # ------------------------------------------------------------------ #
    def _log_and_sync(self, obj: BaseClass, edit: UnitaryEdit, env: simpy.Environment) -> None:
        """
        Append *edit* to obj.runtime history and normalise FilterStore
        membership. No-op for non-LabObjects.
        """
        if obj is None or not _needs_runtime_tracking(obj):
            return
        get_runtime_state(obj, env).recent_edits.append(edit)
        self._sync_filter_stores(obj, env)

    @staticmethod
    def _register_if_new(obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        if obj.instance_iri not in _RESOURCE_MAP:
            _RESOURCE_MAP[obj.instance_iri] = simpy.Resource(env, capacity=1)
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
            logger.debug(f"auto register new object: {obj.__class__.__name__}={obj.instance_iri}")

    def _unregister_object(self, obj: BaseClass):
        """Remove *all* runtime artefacts for a vanished object."""
        if not _needs_runtime_tracking(obj):
            return
        _RESOURCE_MAP.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        _RUNTIME_CACHE.pop(obj.instance_iri, None)
        setattr(obj, "_runtime", None)  # clear pointer

    def _sync_filter_stores(self, obj: BaseClass, env: simpy.Environment):
        """
        After *any* edit we normalise FilterStore membership.
        """
        if not _needs_runtime_tracking(obj):
            return
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        if obj.is_present == {True}:
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)

    def validate_now(self):
        """
        Run SHACL validation at arbitrary times **without** capturing
        violation records or raising. Returns `(conforms, report_graph, _)`.
        """
        if self.shapes_graph is None:
            from rdflib import Graph
            return True, Graph(), ""
        data_graph = KnowledgeGraph.graph()
        overlay = self._build_overlay_graph()
        union = data_graph + overlay
        return validate(
            union,
            shacl_graph=self.shapes_graph,
            ont_graph=self.shapes_graph,
            inference=self.inference,
            advanced=True,
            debug=False,
        )
