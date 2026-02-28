from __future__ import annotations

import pytest
from rdflib import Namespace
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import Has_interrupt_events, MaterialContainer
from libsyn_tools.sim.operation.selector import FilterStoreRegistry, KgQuerySelector
from libsyn_tools.sim.operation.unitary_edit import AddDataProperty, Create


def test_kg_query_selector_fail_fast_when_empty() -> None:
    pool = "VIAL"
    container = MaterialContainer(identifier="https://libsyn-sim/kg/kgq-fail-fast")
    container.has_pool_type.add(pool)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    sim = Simulation([])
    lib = Namespace("https://libsyn-sim/kg/")
    sparql = f"""
    PREFIX lib: <{lib}>
    SELECT ?s WHERE {{
      ?s lib:has_interrupt_events "ready" .
    }}
    """
    selector = KgQuerySelector(
        pool_type=pool,
        sparql=sparql,
        empty_result_policy="fail_fast",
    )

    proc = sim.env.process(selector.resolve(sim.env))
    with pytest.raises(RuntimeError, match="query returned no candidates"):
        sim.env.run(proc)


def test_kg_query_selector_wait_mode_refreshes_candidates() -> None:
    pool = "VIAL"
    container = MaterialContainer(identifier="https://libsyn-sim/kg/kgq-wait")
    container.has_pool_type.add(pool)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    sim = Simulation([])
    lib = Namespace("https://libsyn-sim/kg/")
    sparql = f"""
    PREFIX lib: <{lib}>
    SELECT ?s WHERE {{
      ?s lib:has_interrupt_events "ready" .
    }}
    """
    selector = KgQuerySelector(
        pool_type=pool,
        sparql=sparql,
        empty_result_policy="wait",
        refresh_interval=0.05,
    )

    def _enable_candidate() -> None:
        yield sim.env.timeout(0.2)
        AddDataProperty(
            instance_1_iri=container.identifier,
            property_iri=Has_interrupt_events.predicate_iri,
            data_value="ready",
        ).apply()
        FilterStoreRegistry.put_obj_into_filter_store(container, sim.env)

    proc = sim.env.process(selector.resolve(sim.env))
    sim.env.process(_enable_candidate())
    sim.env.run(until=1.0)

    assert proc.triggered
    iri, req = proc.value
    assert iri == container.identifier
    req.resource.release(req)
