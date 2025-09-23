# ### THIS IS THE START OF CONTENT OF libsyn_tools/sim/effect_engine.py ###
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import List, Union, Iterable as It, Optional, Callable
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
from .lifecycle import LifecycleCallbacks

LIB_SYN = Namespace("https://libsyn-sim/kg/")


class SHACLValidationError(RuntimeError):
    """
    Raised when the knowledge graph violates one or more user-supplied
    SHACL shapes *after* an Operation has been applied (policy violation).
    """
    pass


class EffectEngine:
    """
    Micro-commit of edits + SHACL audit. Supports pluggable overlay providers.
    """

    def __init__(
        self,
        *,
        shapes_graph: ConjunctiveGraph | None = None,
        raise_shacl: bool = False,
        inference: str = "owlrl",
        callbacks: LifecycleCallbacks | None = None,
    ):
        self.shapes_graph = shapes_graph
        self.inference = inference
        self._shacl_violations: list[SHACLViolationRecord] = []
        self.raise_shacl = raise_shacl
        self.callbacks = callbacks

        # registered overlay providers → functions that return an rdflib.Graph
        self._overlay_providers: list[Callable[[], Graph]] = []


    # --- overlay provider registry ---
    def register_overlay_provider(self, provider: Callable[[], Graph]) -> None:
        """
        Register a callable that returns an rdflib.Graph to be unioned into the overlay.
        Providers should be fast and side-effect free.
        """
        self._overlay_providers.append(provider)

    # --- overlay construction (ephemeral) ---
    def _build_overlay_graph(self) -> Graph:
        """
        Build the overlay by calling all registered providers.
        Providers must be fast and side-effect free.
        """
        g = Graph()
        for prov in self._overlay_providers:
            try:
                pg = prov()
                if isinstance(pg, Graph):
                    for t in pg.triples((None, None, None)):
                        g.add(t)
            except Exception as e:
                logger.error(f"Overlay provider failed: {e!r}")
        return g

    # --- SHACL helpers (unchanged) ---
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
        if self.shapes_graph is None:
            return

        data_graph: Graph = KnowledgeGraph.graph()
        overlay_graph: Graph = self._build_overlay_graph()

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
            recs = self._collect_shacl_violations(
                env_now=env.now,
                operation_id=operation_id,
                shacl_report_graph=shacl_report_graph,
                batch_id=batch_id,
                edit_fingerprints=edit_fingerprints,
                seed=seed,
            )
            self._shacl_violations.extend(recs)
            if self.callbacks:
                for rec in recs:
                    self.callbacks.emit_violation(rec)
            if self.raise_shacl:
                raise SHACLValidationError(shacl_report_graph)

    # --- CSV export (unchanged) ---
    def write_shacl_csv(self, filepath: Union[str, Path], include_ttl: bool = False, **to_csv_kwargs) -> None:
        records_as_dicts = [
            rec.model_dump() if include_ttl else rec.model_dump(exclude={"report_graph_ttl"})
            for rec in self._shacl_violations
        ]
        pd.DataFrame.from_records(records_as_dicts).to_csv(Path(filepath), index=False, **to_csv_kwargs)

    # --- Mechanical pre-checks + apply (unchanged semantics) ---
    @staticmethod
    def _fingerprint_edit(edit: UnitaryEdit) -> str:
        dv = edit.data_value
        dv_repr = str(dv) if dv is None or isinstance(dv, (int, float, str, bool)) else f"<{type(dv).__name__}>"
        return "|".join([edit.type.value, edit.instance_1_iri or "", edit.property_iri or "", edit.instance_2_iri or "", dv_repr])

    def _precheck_mechanical(self, *, edits: It[UnitaryEdit], locked_iris: Optional[It[str]] = None) -> None:
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
            return obj is not None

        for e in edits_list:
            t = e.type
            if t is UnitaryEditType.CREATE:
                subj = KnowledgeGraph.get_object_from_lookup(e.instance_1_iri)
                if subj is not None and getattr(subj, "is_present", {False}) == {True}:
                    raise RuntimeError(f"Mechanical check failed: CREATE on present object {e.instance_1_iri}")
            elif t in (UnitaryEditType.ADD_DATA_PROPERTY, UnitaryEditType.CHANGE_DATA_PROPERTY, UnitaryEditType.ANNIHILATE):
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling subject {e.instance_1_iri}")
            elif t in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY):
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling subject {e.instance_1_iri}")
                if not (_exists(e.instance_2_iri) or e.instance_2_iri in creates_set):
                    raise RuntimeError(f"Mechanical check failed: dangling object {e.instance_2_iri}")
                if locked:
                    for iri in (e.instance_1_iri, e.instance_2_iri):
                        if iri not in locked and iri not in creates_set:
                            raise RuntimeError(
                                "Mechanical check failed: write coverage requires lock or create; "
                                f"got edit on {iri!r} not in locked set {locked}"
                            )

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
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
        edits = list(edits)
        batch_id = str(uuid4())
        fingerprints = [self._fingerprint_edit(e) for e in edits]

        self._precheck_mechanical(edits=edits, locked_iris=list(locked_iris or []))

        for edit in edits:
            subj = KnowledgeGraph.get_object_from_lookup(edit.instance_1_iri)
            obj2 = (
                KnowledgeGraph.get_object_from_lookup(edit.instance_2_iri)
                if edit.type in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY)
                and edit.instance_2_iri else None
            )
            self._register_if_new(subj, env)
            if edit.type is not UnitaryEditType.CREATE and obj2:
                self._register_if_new(obj2, env)

            logger.debug(f"Applying edit: {edit.type} – {subj.__class__.__name__}={edit.instance_1_iri}")
            edit.apply()
            self._log_and_sync(subj, edit, env)
            self._log_and_sync(obj2, edit, env)

            if edit.type is UnitaryEditType.CREATE:
                self._register_if_new(subj, env)
            elif edit.type is UnitaryEditType.ANNIHILATE:
                self._unregister_object(subj)

        self._run_shacl_validation(
            env=env,
            operation_id=operation_id,
            batch_id=batch_id,
            edit_fingerprints=fingerprints,
            seed=seed,
        )

    # --- runtime artefacts (unchanged) ---
    def _log_and_sync(self, obj: BaseClass, edit: UnitaryEdit, env: simpy.Environment) -> None:
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
        if not _needs_runtime_tracking(obj):
            return
        _RESOURCE_MAP.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        _RUNTIME_CACHE.pop(obj.instance_iri, None)
        setattr(obj, "_runtime", None)

    def _sync_filter_stores(self, obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        if obj.is_present == {True}:
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)

    def validate_now(self):
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
