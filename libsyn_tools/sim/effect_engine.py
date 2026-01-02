# ### THIS IS THE START OF CONTENT OF libsyn_tools/sim/effect_engine.py ###
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
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

from libsyn_tools.sim.knowledge_graph import BaseClass, MaterialContainer, SimOntology
from libsyn_tools.sim.operation import Operation, UnitaryEdit, UnitaryEditType, FilterStoreRegistry, get_runtime_state
from libsyn_tools.sim.operation.runtime import get_runtime_context, _needs_runtime_tracking
from .effect_shacl import SHACLViolationRecord, _iter_validation_results, _first
from .graph_utils import union_view_many, union_view_for_shacl
from .policy import PolicyBundle
from .lifecycle import LifecycleCallbacks

LIB_SYN = Namespace("https://libsyn-sim/kg/")


class SHACLValidationError(RuntimeError):
    """
    Raised when the knowledge graph violates one or more user-supplied
    SHACL shapes *after* an Operation has been applied (policy violation).
    """
    pass


class EngineMechanicalError(RuntimeError):
    def __init__(
        self,
        message: str,
        violation: SHACLViolationRecord | list[SHACLViolationRecord] | None = None,
    ):
        super().__init__(message, violation)
        if violation is None:
            self.violations = []
        else:
            self.violations = violation if isinstance(violation, list) else [violation]


class ContractViolationError(RuntimeError):
    def __init__(self, result: TransactionResult):
        super().__init__("Transaction aborted due to SHACL contract violation", result)
        self.result = result


@dataclass(frozen=True)
class TransactionResult:
    batch_id: str
    edit_fingerprints: list[str]
    edit_descriptions: list[str]
    committed: bool
    violations: list[SHACLViolationRecord]


@dataclass
class _ObjectSnapshot:
    is_present: set
    fields: dict[str, set]
    runtime_registered: bool
    runtime_cache_present: bool
    runtime_in_filter_store: bool
    pool_type: str | None


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
        policy: PolicyBundle | None = None,
    ):
        self.shapes_graph = shapes_graph
        self.inference = inference
        self._shacl_violations: list[SHACLViolationRecord] = []
        self.raise_shacl = raise_shacl
        self.callbacks = callbacks
        self.policy = policy

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
    def _collect_overlay_graphs(self) -> list[Graph]:
        """
        Collect overlay graphs by calling all registered providers.
        Providers must be fast and side-effect free.
        """
        graphs: list[Graph] = []
        for prov in self._overlay_providers:
            try:
                pg = prov()
                if isinstance(pg, Graph):
                    graphs.append(pg)
            except Exception as e:
                logger.error(f"Overlay provider failed: {e!r}")
        return graphs

    def build_query_graph(self) -> Graph:
        """
        Build a union graph of the base KG plus all overlay providers.
        This is side-effect free and intended for query/selection.
        """
        data_graph: Graph = KnowledgeGraph.graph()
        overlay_graphs = self._collect_overlay_graphs()
        return union_view_many([data_graph, *overlay_graphs])

    # --- SHACL helpers (unchanged) ---
    def _collect_shacl_violations(
        self,
        *,
        env_now: float,
        operation_id: str,
        shacl_report_graph: Graph,
        batch_id: Optional[str],
        edit_fingerprints: Optional[list[str]],
        edit_descriptions: Optional[list[str]],
        seed: Optional[int]
    ) -> list[SHACLViolationRecord]:
        records: list[SHACLViolationRecord] = []
        for vr in _iter_validation_results(shacl_report_graph):
            shape_iri = _first(shacl_report_graph, vr, SH.sourceShape)
            policy = self.policy.rule_for_shape(shape_iri) if self.policy else None
            rec = SHACLViolationRecord(
                sim_time=env_now,
                operation_id=operation_id,
                origin="SHACL",
                severity=policy.severity if policy else "soft",
                disposition=policy.disposition if policy else "committed",
                batch_id=batch_id,
                edit_fingerprints=edit_fingerprints,
                edit_descriptions=edit_descriptions,
                seed=seed,
                shape_iri=shape_iri,
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
        edit_descriptions: Optional[list[str]],
        seed: Optional[int]
    ) -> None:
        if self.shapes_graph is None:
            return

        data_graph: Graph = KnowledgeGraph.graph()
        overlay_graphs = self._collect_overlay_graphs()
        union_graph = union_view_for_shacl([data_graph, *overlay_graphs])

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
                edit_descriptions=edit_descriptions,
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

    def get_violation_records(self) -> list[SHACLViolationRecord]:
        return list(self._shacl_violations)

    # --- Mechanical pre-checks + apply (unchanged semantics) ---
    @staticmethod
    def _fingerprint_edit(edit: UnitaryEdit) -> str:
        dv = edit.data_value
        dv_repr = str(dv) if dv is None or isinstance(dv, (int, float, str, bool)) else f"<{type(dv).__name__}>"
        return "|".join([edit.type.value, edit.instance_1_iri or "", edit.property_iri or "", edit.instance_2_iri or "", dv_repr])

    def _precheck_mechanical(
        self,
        *,
        edits: It[UnitaryEdit],
        locked_iris: Optional[It[str]] = None,
        env_now: float,
        operation_id: str,
        batch_id: str,
        edit_fingerprints: list[str],
        edit_descriptions: list[str],
        seed: Optional[int],
    ) -> None:
        locked = set(locked_iris or [])
        edits_list = list(edits)
        creates = [e.instance_1_iri for e in edits_list if e.type is UnitaryEditType.CREATE]
        if len(creates) != len(set(creates)):
            self._raise_mechanical(
                f"Mechanical check failed: duplicate CREATE in batch: {creates}",
                env_now=env_now,
                operation_id=operation_id,
                batch_id=batch_id,
                edit_fingerprints=edit_fingerprints,
                edit_descriptions=edit_descriptions,
                seed=seed,
            )
        creates_set = set(creates)

        def _exists(iri: Optional[str]) -> bool:
            if not iri:
                return False
            try:
                obj = KnowledgeGraph.get_object_from_lookup(iri)
            except Exception:
                return False
            return obj is not None

        def _is_runtime_tracked(iri: Optional[str]) -> bool:
            if not iri:
                return False
            try:
                obj = KnowledgeGraph.get_object_from_lookup(iri)
            except Exception:
                return False
            if obj is None:
                return False
            return _needs_runtime_tracking(obj)

        def _require_lock_if_runtime_tracked(iri: Optional[str]) -> None:
            if not iri or iri in creates_set:
                return
            if _is_runtime_tracked(iri) and iri not in locked:
                self._raise_mechanical(
                    "Mechanical check failed: write coverage requires lock or create; "
                    f"got edit on {iri!r} not in locked set {locked}",
                    env_now=env_now,
                    operation_id=operation_id,
                    batch_id=batch_id,
                    edit_fingerprints=edit_fingerprints,
                    edit_descriptions=edit_descriptions,
                    seed=seed,
                )

        for e in edits_list:
            t = e.type
            if t is UnitaryEditType.CREATE:
                subj = KnowledgeGraph.get_object_from_lookup(e.instance_1_iri)
                if subj is not None and getattr(subj, "is_present", {False}) == {True}:
                    self._raise_mechanical(
                        f"Mechanical check failed: CREATE on present object {e.instance_1_iri}",
                        env_now=env_now,
                        operation_id=operation_id,
                        batch_id=batch_id,
                        edit_fingerprints=edit_fingerprints,
                        edit_descriptions=edit_descriptions,
                        seed=seed,
                    )
            elif t in (
                UnitaryEditType.ADD_DATA_PROPERTY,
                UnitaryEditType.CHANGE_DATA_PROPERTY,
                UnitaryEditType.REMOVE_DATA_PROPERTY,
                UnitaryEditType.ANNIHILATE,
            ):
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    self._raise_mechanical(
                        f"Mechanical check failed: dangling subject {e.instance_1_iri}",
                        env_now=env_now,
                        operation_id=operation_id,
                        batch_id=batch_id,
                        edit_fingerprints=edit_fingerprints,
                        edit_descriptions=edit_descriptions,
                        seed=seed,
                    )
                _require_lock_if_runtime_tracked(e.instance_1_iri)
            elif t in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY):
                if not (_exists(e.instance_1_iri) or e.instance_1_iri in creates_set):
                    self._raise_mechanical(
                        f"Mechanical check failed: dangling subject {e.instance_1_iri}",
                        env_now=env_now,
                        operation_id=operation_id,
                        batch_id=batch_id,
                        edit_fingerprints=edit_fingerprints,
                        edit_descriptions=edit_descriptions,
                        seed=seed,
                    )
                if not (_exists(e.instance_2_iri) or e.instance_2_iri in creates_set):
                    self._raise_mechanical(
                        f"Mechanical check failed: dangling object {e.instance_2_iri}",
                        env_now=env_now,
                        operation_id=operation_id,
                        batch_id=batch_id,
                        edit_fingerprints=edit_fingerprints,
                        edit_descriptions=edit_descriptions,
                        seed=seed,
                    )
                _require_lock_if_runtime_tracked(e.instance_1_iri)
                _require_lock_if_runtime_tracked(e.instance_2_iri)

    def _raise_mechanical(
        self,
        message: str,
        *,
        env_now: float,
        operation_id: str,
        batch_id: str,
        edit_fingerprints: list[str],
        edit_descriptions: list[str],
        seed: Optional[int],
    ) -> None:
        violation = SHACLViolationRecord(
            sim_time=env_now,
            operation_id=operation_id,
            origin="ENGINE",
            severity="hard",
            disposition="aborted",
            batch_id=batch_id,
            edit_fingerprints=edit_fingerprints,
            edit_descriptions=edit_descriptions,
            seed=seed,
            message=message,
        )
        raise EngineMechanicalError(message, violation)

    def _snapshot_objects(
        self,
        *,
        env: simpy.Environment,
        edits: list[UnitaryEdit],
    ) -> dict[str, _ObjectSnapshot]:
        field_names_by_iri: dict[str, set[str]] = {}
        affected_iris: set[str] = set()
        for edit in edits:
            affected_iris.add(edit.instance_1_iri)
            if edit.instance_2_iri:
                affected_iris.add(edit.instance_2_iri)
            if edit.property_iri is None:
                continue
            if edit.type in (
                UnitaryEditType.ADD_DATA_PROPERTY,
                UnitaryEditType.CHANGE_DATA_PROPERTY,
                UnitaryEditType.REMOVE_DATA_PROPERTY,
            ):
                data_prop = SimOntology.data_property_lookup[edit.property_iri]
                field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
                field_names_by_iri.setdefault(edit.instance_1_iri, set()).add(field_name)
            elif edit.type in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY):
                obj_prop = SimOntology.object_property_lookup[edit.property_iri]
                field_name = obj_prop.__name__[0].lower() + obj_prop.__name__[1:]
                field_names_by_iri.setdefault(edit.instance_1_iri, set()).add(field_name)

        ctx = get_runtime_context(env)
        snapshots: dict[str, _ObjectSnapshot] = {}
        for iri in affected_iris:
            obj = KnowledgeGraph.get_object_from_lookup(iri)
            if obj is None:
                continue
            fields: dict[str, set] = {}
            for name in field_names_by_iri.get(iri, set()):
                fields[name] = set(getattr(obj, name))
            runtime_registered = False
            runtime_cache_present = False
            runtime_in_filter_store = False
            pool_type = None
            if _needs_runtime_tracking(obj):
                pool_type = next(iter(obj.has_pool_type), None)
                runtime_registered = obj.instance_iri in ctx.resource_map
                runtime_cache_present = obj.instance_iri in ctx.runtime_cache
                if pool_type:
                    store = ctx.filter_stores.get(pool_type)
                    runtime_in_filter_store = bool(store and obj in store.items)
            snapshots[iri] = _ObjectSnapshot(
                is_present=set(obj.is_present),
                fields=fields,
                runtime_registered=runtime_registered,
                runtime_cache_present=runtime_cache_present,
                runtime_in_filter_store=runtime_in_filter_store,
                pool_type=pool_type,
            )
        return snapshots

    def _restore_objects(
        self,
        *,
        env: simpy.Environment,
        snapshots: dict[str, _ObjectSnapshot],
    ) -> None:
        ctx = get_runtime_context(env)
        for iri, snap in snapshots.items():
            obj = KnowledgeGraph.get_object_from_lookup(iri)
            if obj is None:
                continue
            obj.is_present = set(snap.is_present)
            for name, values in snap.fields.items():
                setattr(obj, name, set(values))

            if not _needs_runtime_tracking(obj):
                continue

            if snap.runtime_registered:
                if obj.instance_iri not in ctx.resource_map:
                    ctx.resource_map[obj.instance_iri] = simpy.Resource(env, capacity=1)
            else:
                ctx.resource_map.pop(obj.instance_iri, None)

            if snap.runtime_cache_present:
                if obj.instance_iri not in ctx.runtime_cache:
                    if obj.instance_iri not in ctx.resource_map:
                        ctx.resource_map[obj.instance_iri] = simpy.Resource(env, capacity=1)
                    get_runtime_state(obj, env)
            else:
                ctx.runtime_cache.pop(obj.instance_iri, None)
                setattr(obj, "_runtime", None)

            if snap.pool_type:
                store = ctx.filter_stores.get(snap.pool_type)
                if snap.runtime_in_filter_store:
                    if store is None:
                        store = FilterStoreRegistry.get_filter_store(snap.pool_type, env)
                    if obj not in store.items:
                        store.items.append(obj)
                elif store and obj in store.items:
                    store.items.remove(obj)

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
        return action.operation_effects.copy()

    def apply_tx(
        self,
        edits: Iterable[UnitaryEdit],
        env: simpy.Environment,
        *,
        operation_id: str,
        locked_iris: Optional[Iterable[str]] = None,
        seed: Optional[int] = None,
    ) -> TransactionResult:
        """
        Apply edits transactionally: precheck -> apply -> SHACL audit -> commit/rollback.

        Default behavior is audit-only (committed=True) unless a policy marks any
        SHACL violation with disposition="aborted", which triggers rollback and
        returns committed=False.
        """
        edits = list(edits)
        batch_id = str(uuid4())
        fingerprints = [self._fingerprint_edit(e) for e in edits]
        descriptions = [e.describe() for e in edits]
        snapshots = self._snapshot_objects(env=env, edits=edits)

        try:
            self._precheck_mechanical(
                edits=edits,
                locked_iris=list(locked_iris or []),
                env_now=env.now,
                operation_id=operation_id,
                batch_id=batch_id,
                edit_fingerprints=fingerprints,
                edit_descriptions=descriptions,
                seed=seed,
            )
        except EngineMechanicalError as err:
            self._shacl_violations.extend(err.violations)
            if self.callbacks:
                for rec in err.violations:
                    self.callbacks.emit_violation(rec)
            raise

        before_len = len(self._shacl_violations)
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
                self._unregister_object(subj, env)

        self._run_shacl_validation(
            env=env,
            operation_id=operation_id,
            batch_id=batch_id,
            edit_fingerprints=fingerprints,
            edit_descriptions=descriptions,
            seed=seed,
        )
        new_violations = self._shacl_violations[before_len:]
        if any(v.disposition == "aborted" for v in new_violations):
            self._restore_objects(env=env, snapshots=snapshots)
            return TransactionResult(
                batch_id=batch_id,
                edit_fingerprints=fingerprints,
                edit_descriptions=descriptions,
                committed=False,
                violations=new_violations,
            )
        return TransactionResult(
            batch_id=batch_id,
            edit_fingerprints=fingerprints,
            edit_descriptions=descriptions,
            committed=True,
            violations=new_violations,
        )

    def apply(
        self,
        edits: Iterable[UnitaryEdit],
        env: simpy.Environment,
        *,
        operation_id: str,
        locked_iris: Optional[Iterable[str]] = None,
        seed: Optional[int] = None,
    ) -> None:
        result = self.apply_tx(
            edits,
            env,
            operation_id=operation_id,
            locked_iris=locked_iris,
            seed=seed,
        )
        if not result.committed:
            raise ContractViolationError(result)

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
        ctx = get_runtime_context(env)
        if obj.instance_iri not in ctx.resource_map:
            ctx.resource_map[obj.instance_iri] = simpy.Resource(env, capacity=1)
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
            logger.debug(f"auto register new object: {obj.__class__.__name__}={obj.instance_iri}")

    def _unregister_object(self, obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        ctx = get_runtime_context(env, create=False)
        ctx.resource_map.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj, env)
        ctx.runtime_cache.pop(obj.instance_iri, None)
        setattr(obj, "_runtime", None)

    def _sync_filter_stores(self, obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        FilterStoreRegistry.remove_obj_from_filter_store(obj, env)
        if obj.is_present == {True}:
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)

    def validate_now(self):
        if self.shapes_graph is None:
            from rdflib import Graph
            return True, Graph(), ""
        data_graph = KnowledgeGraph.graph()
        overlay_graphs = self._collect_overlay_graphs()
        union = union_view_for_shacl([data_graph, *overlay_graphs])
        return validate(
            union,
            shacl_graph=self.shapes_graph,
            ont_graph=self.shapes_graph,
            inference=self.inference,
            advanced=True,
            debug=False,
        )
