from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry, LiteralSelector
from libsyn_tools.sim.operation.unitary_edit import Create


def test_literal_selector_retries_when_target_becomes_non_present_before_lock() -> None:
    pool = "VIAL"
    container = MaterialContainer(identifier="https://libsyn-sim/kg/literal-race")
    container.has_pool_type.add(pool)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    sim = Simulation([])
    rs = get_runtime_state(container, sim.env)
    held = rs.lock.request()
    sim.env.run(until=held)

    selector = LiteralSelector(container.identifier)

    def _flip_presence() -> None:
        yield sim.env.timeout(0.1)
        container.is_present = {False}
        FilterStoreRegistry.remove_obj_from_filter_store(container, sim.env)
        rs.lock.release(held)
        yield sim.env.timeout(0.1)
        container.is_present = {True}
        FilterStoreRegistry.put_obj_into_filter_store(container, sim.env)

    proc = sim.env.process(selector.resolve(sim.env))
    sim.env.process(_flip_presence())
    sim.env.run(until=1.0)

    assert proc.triggered
    iri, req = proc.value
    assert iri == container.identifier
    req.resource.release(req)
