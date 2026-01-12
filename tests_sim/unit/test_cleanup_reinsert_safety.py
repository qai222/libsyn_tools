from __future__ import annotations

import inspect
import simpy
from loguru import logger
from pydantic import Field
import pytest
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD

from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject, MaterialContainer, Has_capacity
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.runtime import get_runtime_context, get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import ChangeDataProperty, Create, UnitaryEdit
from libsyn_tools.sim.policy import PolicyBundle, PolicyRule


class _OpTwoLocks(Operation):
    participant_a: str = Field(...)
    participant_b: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_cleanup_releases_all_locks_even_if_pool_reinsert_raises(
    env: simpy.Environment,
) -> None:
    obj_a = LabObject()
    obj_a.is_present = {True}
    obj_a.has_pool_type.add("POOL_CLEANUP_A")
    KnowledgeGraph.get_object_from_lookup(obj_a.identifier)
    Create(instance_1_iri=obj_a.identifier).apply()

    obj_b = LabObject()
    obj_b.is_present = {True}
    obj_b.has_pool_type.add("POOL_CLEANUP_B")
    KnowledgeGraph.get_object_from_lookup(obj_b.identifier)
    Create(instance_1_iri=obj_b.identifier).apply()

    engine = EffectEngine()
    engine._register_if_new(obj_a, env)
    engine._register_if_new(obj_b, env)

    FilterStoreRegistry.put_obj_into_filter_store(obj_a, env)
    FilterStoreRegistry.put_obj_into_filter_store(obj_b, env)

    op = _OpTwoLocks(participant_a=obj_a.identifier, participant_b=obj_b.identifier)
    env.run(op.pre_act(env))

    obj_b.has_pool_type.add("POOL_CLEANUP_B_EXTRA")

    records: list[str] = []
    sink_id = logger.add(lambda msg: records.append(str(msg)), level="WARNING")

    op.cleanup(env)

    logger.remove(sink_id)

    rs_a = get_runtime_state(obj_a, env)
    rs_b = get_runtime_state(obj_b, env)
    assert rs_a.lock.count == 0
    assert rs_b.lock.count == 0
    assert not rs_a.lock.queue
    assert not rs_b.lock.queue

    ctx = get_runtime_context(env, create=False)
    assert all(obj_b not in store.items for store in ctx.filter_stores.values())
    assert any("FilterStore reinsertion failed" in record for record in records)


def test_rollback_does_not_strand_locks_when_reinsert_fails(
    env: simpy.Environment,
    monkeypatch,
) -> None:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.RollbackShapeReinsertSafe)
    ttl = f"""
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:RollbackShapeReinsertSafe a sh:NodeShape ;
       sh:targetSubjectsOf lib:has_capacity ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select \"\"\"
    SELECT ?this WHERE {{
      ?this lib:has_capacity ?cap .
      FILTER(xsd:double(?cap) > 0)
    }}
    \"\"\" ;
       ] .
    """
    shape_graph = Graph().parse(data=ttl, format="turtle")

    container = MaterialContainer(identifier="rollback-reinsert-container")
    container.is_present = {True}
    container.has_pool_type.add("POOL_ROLLBACK_REINSERT")
    container.has_capacity.add(1.0)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    engine = EffectEngine(shapes_graph=shape_graph)
    engine._register_if_new(container, env)
    get_runtime_state(container, env)
    FilterStoreRegistry.put_obj_into_filter_store(container, env)
    engine.policy = PolicyBundle(
        per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")}
    )

    original_put = FilterStoreRegistry.put_obj_into_filter_store

    def _boom(*args, **kwargs) -> None:
        if any(frame.function == "safe_put_obj_into_filter_store" for frame in inspect.stack()):
            raise RuntimeError("boom")
        return original_put(*args, **kwargs)

    monkeypatch.setattr(FilterStoreRegistry, "put_obj_into_filter_store", _boom)

    bad_edit = ChangeDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_capacity.predicate_iri,
        data_value=0.0,
    )
    records: list[str] = []
    sink_id = logger.add(lambda msg: records.append(str(msg)), level="WARNING")

    def _apply_fail() -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(bad_edit, "apply", _apply_fail)

    with pytest.raises(EngineMechanicalError):
        engine.apply_tx(
            edits=[bad_edit],
            env=env,
            operation_id="rollback-reinsert",
            locked_iris=[container.identifier],
        )

    logger.remove(sink_id)

    rs = get_runtime_state(container, env)
    assert rs.lock.count == 0
    assert not rs.lock.queue
    assert any("FilterStore reinsertion failed" in record for record in records)
