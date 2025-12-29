# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_selector.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import LabObject
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import (
    AttributeSelector,
    LiteralSelector,
    FilterStoreRegistry,
)


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
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_selector.py ###
