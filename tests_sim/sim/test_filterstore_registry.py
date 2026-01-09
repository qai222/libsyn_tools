from __future__ import annotations

import simpy

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry


def test_put_respects_queued_lock_request() -> None:
    obj = MaterialContainer(identifier="queue-object")
    obj.has_pool_type.add("POOL_QUEUE")
    obj.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    sim = Simulation([])
    rs = get_runtime_state(obj, sim.env)

    def holder(env: simpy.Environment) -> simpy.events.Event:
        req = rs.lock.request()
        yield req
        yield env.timeout(1.0)
        rs.lock.release(req)

    sim.env.process(holder(sim.env))
    sim.env.run(until=0.1)
    FilterStoreRegistry.remove_obj_from_filter_store(obj, sim.env)

    pending = rs.lock.request()
    assert len(rs.lock.queue) > 0

    FilterStoreRegistry.put_obj_into_filter_store(obj, sim.env)
    store = FilterStoreRegistry.get_filter_store("POOL_QUEUE", sim.env)
    assert obj not in store.items

    pending.cancel()


def test_put_dedupes_store_entries() -> None:
    obj = MaterialContainer(identifier="dedupe-object")
    obj.has_pool_type.add("POOL_DEDUPE")
    obj.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    sim = Simulation([])

    FilterStoreRegistry.put_obj_into_filter_store(obj, sim.env)
    FilterStoreRegistry.put_obj_into_filter_store(obj, sim.env)

    store = FilterStoreRegistry.get_filter_store("POOL_DEDUPE", sim.env)
    assert store.items.count(obj) == 1


def test_remove_safe_with_pending_get() -> None:
    obj = MaterialContainer(identifier="remove-object")
    obj.has_pool_type.add("POOL_REMOVE")
    obj.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    sim = Simulation([])
    store = FilterStoreRegistry.get_filter_store("POOL_REMOVE", sim.env)
    get_ev = store.get(filter=lambda candidate: candidate is obj)
    FilterStoreRegistry.remove_obj_from_filter_store(obj, sim.env)
    FilterStoreRegistry.put_obj_into_filter_store(obj, sim.env)
    sim.env.run(until=0.1)
    assert get_ev.triggered
