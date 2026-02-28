from __future__ import annotations

import pytest
from rdflib import Namespace
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.selector import KgQuerySelector
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class _SelectTargetByKg(Operation):
    participant_target: StrOrSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_run_until_none_rejects_wait_mode_kg_query_selector() -> None:
    pool = "VIAL"
    container = MaterialContainer(identifier="https://libsyn-sim/kg/kgq-wait-unsat")
    container.has_pool_type.add(pool)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    lib = Namespace("https://libsyn-sim/kg/")
    sparql = f"""
    PREFIX lib: <{lib}>
    SELECT ?s WHERE {{
      ?s lib:has_interrupt_events "never-satisfied" .
    }}
    """
    op = _SelectTargetByKg(
        identifier="kgq-wait-unsat-op",
        participant_target=KgQuerySelector(
            pool_type=pool,
            sparql=sparql,
            empty_result_policy="wait",
            refresh_interval=0.05,
        ),
        temporal_cost=0.0,
    )
    sim = Simulation([op])

    with pytest.raises(ValueError, match="does not allow KgQuerySelector\\(empty_result_policy='wait'\\)"):
        sim.run(until=None)
