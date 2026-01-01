# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_wait_preset.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.unitary_edit import Create
from libsyn_tools.sim.operation_preset.wait import Wait


def _make_module(identifier: str) -> MaterialContainer:
    module = MaterialContainer(identifier=identifier)
    module.has_pool_type.add("MODULE")
    KnowledgeGraph.get_object_from_lookup(module.identifier)
    Create(instance_1_iri=module.identifier).apply()
    return module


def _timestamp_for(sim: Simulation, op_id: str, event_type: str) -> float:
    for rec in sim.history_log:
        if rec.operation_id == op_id and rec.event_type == event_type:
            return rec.timestamp
    raise AssertionError(f"Missing {event_type} for {op_id}")


def test_wait_serializes_on_same_module() -> None:
    module = _make_module("M1")
    op_a = Wait(identifier="wait-a", participant_resource=module.identifier, temporal_cost=1.0)
    op_b = Wait(identifier="wait-b", participant_resource=module.identifier, temporal_cost=1.0)
    sim = Simulation([op_a, op_b])
    sim.run()

    start_a = _timestamp_for(sim, "wait-a", "OPERATION_START")
    end_a = _timestamp_for(sim, "wait-a", "OPERATION_END")
    start_b = _timestamp_for(sim, "wait-b", "OPERATION_START")

    assert start_a == 0.0
    assert end_a == 1.0
    assert start_b >= end_a


def test_wait_without_participant_does_not_block() -> None:
    module = _make_module("M2")
    free_wait = Wait(identifier="wait-free", temporal_cost=2.0)
    locked_wait = Wait(identifier="wait-locked", participant_resource=module.identifier, temporal_cost=1.0)
    sim = Simulation([free_wait, locked_wait])
    sim.run()

    start_free = _timestamp_for(sim, "wait-free", "OPERATION_START")
    start_locked = _timestamp_for(sim, "wait-locked", "OPERATION_START")

    assert start_free == start_locked == 0.0
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_wait_preset.py ###
