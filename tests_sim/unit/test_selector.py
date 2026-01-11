# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_selector.py ###
from __future__ import annotations

import simpy
import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import LabObject, Has_interrupt_events
from libsyn_tools.sim.operation.runtime import get_runtime_context, get_runtime_state
from libsyn_tools.sim.operation.selector import (
    AttributeSelector,
    HistorySelector,
    LiteralSelector,
    FilterStoreRegistry,
)
from libsyn_tools.sim.operation.unitary_edit import AddDataProperty, Create


def _make_pool_obj(pool: str) -> LabObject:
    o = LabObject()
    o.has_pool_type.add(pool)
    return o


def test_literal_selector_locks(env: simpy.Environment):
    pool = "VIAL"
    o = _make_pool_obj(pool)
    # Register resource + put in pool
    eng = EffectEngine()
    KnowledgeGraph.get_object_from_lookup(o.identifier)
    Create(instance_1_iri=o.identifier).apply()
    eng._register_if_new(o, env)

    sel = LiteralSelector(o.identifier)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    iri, req = proc.value

    assert iri == o.identifier
    # lock held
    rs = get_runtime_state(o, env)
    assert req.resource is rs.lock
    rs.lock.release(req)


def test_attribute_selector_from_pool(env: simpy.Environment):
    pool = "VIAL"
    eng = EffectEngine()

    o1 = _make_pool_obj(pool)
    o2 = _make_pool_obj(pool)
    for o in (o1, o2):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
        eng._register_if_new(o, env)

    # Predicate: accept anything (keep it simple)
    def pred(obj: LabObject) -> bool:
        return True

    sel = AttributeSelector(pool_type=pool, predicate=pred)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    iri, req = proc.value

    assert iri in {o1.identifier, o2.identifier}
    # Clean up: release; post_act elsewhere would reinsert into pool as needed
    rs = get_runtime_state(KnowledgeGraph.get_object_from_lookup(iri), env)
    rs.lock.release(req)

    # Confirm pool exists and contains at least one instance
    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert (o1 in store.items) or (o2 in store.items)


def test_selector_predicate_exception_restores_pool(env: simpy.Environment) -> None:
    pool = "VIAL_PREDICATE_ERROR"
    eng = EffectEngine()

    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng._register_if_new(obj, env)

    def pred(_obj: LabObject) -> bool:
        raise RuntimeError("predicate failure")

    sel = AttributeSelector(pool_type=pool, predicate=pred)
    proc = env.process(sel.resolve(env))
    with pytest.raises(RuntimeError, match="predicate failure"):
        env.run(proc)

    store = FilterStoreRegistry.get_filter_store(pool, env)
    rs = get_runtime_state(obj, env)
    assert obj in store.items
    assert rs.lock.count == 0


def test_filter_store_isolated_between_envs():
    pool = "VIAL"
    eng = EffectEngine()

    env_one = simpy.Environment()
    env_two = simpy.Environment()

    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng._register_if_new(obj, env_one)

    store_one = FilterStoreRegistry.get_filter_store(pool, env_one)
    store_two = FilterStoreRegistry.get_filter_store(pool, env_two)

    assert store_one is not store_two
    assert obj in store_one.items
    assert obj not in store_two.items

    ctx_one = get_runtime_context(env_one, create=False)
    ctx_two = get_runtime_context(env_two, create=False)
    assert ctx_one is not ctx_two


def test_filter_store_membership_tracks_lock_state(env: simpy.Environment):
    pool = "VIAL"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert obj in store.items

    sel = LiteralSelector(obj.identifier)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    _, req = proc.value

    assert obj not in store.items
    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj not in store.items

    rs = get_runtime_state(obj, env)
    rs.lock.release(req)
    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj in store.items


def test_effect_engine_sync_skips_locked_objects(env: simpy.Environment):
    pool = "VIAL"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    sel = LiteralSelector(obj.identifier)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    _, req = proc.value

    edit = AddDataProperty(
        instance_1_iri=obj.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="LOCKED",
    )
    eng.apply([edit], env, operation_id="lock-test", locked_iris=[obj.identifier])

    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert obj not in store.items

    rs = get_runtime_state(obj, env)
    rs.lock.release(req)
    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj in store.items


def test_filter_store_rejects_queueing_locks(env: simpy.Environment) -> None:
    pool = "POOL_QUEUE_LOCK"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert obj in store.items
    store.items.remove(obj)

    rs = get_runtime_state(obj, env)
    req1 = rs.lock.request()
    req2 = rs.lock.request()
    env.run(until=req1)
    rs.lock.release(req1)

    assert rs.lock.count == 0
    assert rs.lock.queue

    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj not in store.items

    req2.cancel()


def test_literal_selector_heals_missing_store_entry(env: simpy.Environment) -> None:
    pool = "POOL_LITERAL_HEAL"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    store = FilterStoreRegistry.get_filter_store(pool, env)
    store.items.remove(obj)

    sel = LiteralSelector(obj.identifier)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    iri, req = proc.value

    assert iri == obj.identifier
    rs = get_runtime_state(obj, env)
    rs.lock.release(req)
    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj in store.items


def test_history_selector_accepts_env_predicate(env: simpy.Environment) -> None:
    pool = "POOL_HISTORY_ENV"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    def pred(candidate: LabObject, sim_env: simpy.Environment) -> bool:
        assert sim_env is env
        _ = get_runtime_state(candidate, sim_env)
        return True

    sel = HistorySelector(pool_type=pool, predicate=pred)
    proc = env.process(sel.resolve(env))
    env.run(proc)
    iri, req = proc.value

    assert iri == obj.identifier
    rs = get_runtime_state(obj, env)
    rs.lock.release(req)
    FilterStoreRegistry.put_obj_into_filter_store(obj, env)
    assert obj in FilterStoreRegistry.get_filter_store(pool, env).items


def test_predicate_exception_returns_object_to_pool(env: simpy.Environment) -> None:
    pool = "POOL_PRED_RAISE"
    obj = _make_pool_obj(pool)
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    eng = EffectEngine()
    eng._register_if_new(obj, env)

    calls = {"count": 0}

    def pred(candidate: LabObject) -> bool:
        calls["count"] += 1
        if calls["count"] == 1:
            return True
        raise RuntimeError("predicate boom")

    sel = AttributeSelector(pool_type=pool, predicate=pred)
    proc = env.process(sel.resolve(env))
    with pytest.raises(RuntimeError):
        env.run(proc)

    rs = get_runtime_state(obj, env)
    assert rs.lock.count == 0
    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert obj in store.items
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_selector.py ###
