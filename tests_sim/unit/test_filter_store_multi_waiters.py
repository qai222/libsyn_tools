"""Concurrency regression tests for FilterStore + Selector.

These focus on multi-waiter behavior to catch deadlocks / starvation
introduced by direct FilterStore.items manipulation.
"""

import simpy

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_context, get_runtime_state
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import Create
from twa.data_model.base_ontology import KnowledgeGraph


def _make_present_container(*, pool: str) -> MaterialContainer:
    c = MaterialContainer()
    c.has_pool_type.add(pool)
    c.has_capacity.add(100.0)
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    Create(instance_1_iri=c.identifier).apply()
    return c


def test_attribute_selector_two_waiters_get_distinct_resources_and_restore_store():
    env = simpy.Environment()
    get_runtime_context(env)

    engine = EffectEngine(shapes_graph=None)

    pool = "POOL_MULTI_WAITER"
    c1 = _make_present_container(pool=pool)
    c2 = _make_present_container(pool=pool)
    engine._register_if_new(c1, env)
    engine._register_if_new(c2, env)

    store = FilterStoreRegistry.get_filter_store(pool, env)
    assert len(store.items) >= 2  # may include other pooled objects from other tests

    sel = AttributeSelector(pool, predicate=lambda obj: True)

    results: list[tuple[str, simpy.events.Event]] = []

    def worker():
        iri, req = yield env.process(sel.resolve(env))
        results.append((iri, req))

    # Two concurrent waiters.
    env.process(worker())
    env.process(worker())
    env.run()

    assert len(results) == 2
    assert results[0][0] != results[1][0]

    # Both selected objects are removed from the store while locked.
    assert len([o for o in store.items if o.identifier in {results[0][0], results[1][0]}]) == 0

    # Release locks and return them to the store.
    for iri, req in results:
        obj = KnowledgeGraph.get_object_from_lookup(iri)
        rs = get_runtime_state(obj, env)
        assert rs.lock.count == 1  # held
        req.resource.release(req)
        assert rs.lock.count == 0
        FilterStoreRegistry.put_obj_into_filter_store(obj, env)

    assert len([o for o in store.items if o.identifier in {results[0][0], results[1][0]}]) == 2
