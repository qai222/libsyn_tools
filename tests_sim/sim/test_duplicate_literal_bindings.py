# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_duplicate_literal_bindings.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer, canonical_iri
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import LiteralSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class DuplicateLiteralOp(Operation):
    participant_left: str | LiteralSelector
    participant_right: str | LiteralSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _world_one_container() -> MaterialContainer:
    container = MaterialContainer(identifier="dup-container")
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()
    return container


def test_duplicate_literal_strings_do_not_block(env: simpy.Environment):
    container = _world_one_container()
    full_iri = str(canonical_iri(container.identifier))

    op = DuplicateLiteralOp(
        identifier="dup-op-strings",
        participant_left=full_iri,
        participant_right=container.identifier,
    )
    sim = Simulation([op])
    sim.run(until=1.0)

    end_events = [
        r for r in sim.history_log
        if r.operation_id == op.identifier and r.event_type == "OPERATION_END"
    ]
    assert end_events


def test_duplicate_literal_selectors_do_not_block(env: simpy.Environment):
    container = _world_one_container()
    full_iri = str(canonical_iri(container.identifier))

    op = DuplicateLiteralOp(
        identifier="dup-op-selectors",
        participant_left=LiteralSelector(full_iri),
        participant_right=LiteralSelector(full_iri),
    )
    sim = Simulation([op])
    sim.run(until=1.0)

    end_events = [
        r for r in sim.history_log
        if r.operation_id == op.identifier and r.event_type == "OPERATION_END"
    ]
    assert end_events


def test_duplicate_resources_deduped_for_provenance(env: simpy.Environment):
    container = _world_one_container()

    op = DuplicateLiteralOp(
        identifier="dup-op-resources",
        participant_left=container.identifier,
        participant_right=container.identifier,
    )
    sim = Simulation([op])
    sim.run(until=1.0)

    assert op.resources == [container.identifier]
    runtime_state = get_runtime_state(container, sim.env)
    assert list(runtime_state.recent_operations).count(op) == 1
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_duplicate_literal_bindings.py ###
