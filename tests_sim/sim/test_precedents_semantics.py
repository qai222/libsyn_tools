# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_precedents_semantics.py ###
from __future__ import annotations

import simpy
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer, Is_directly_contained_by
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import AddObjectProperty, Create, UnitaryEdit


class AbortOp(Operation):
    participant_container: str
    missing_iri: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddObjectProperty(
                instance_1_iri=self.participant_container,
                instance_2_iri=self.missing_iri,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


class NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_precedent_abort_does_not_block_dependents():
    base = "https://libsyn-sim/kg/"
    container = MaterialContainer(identifier=f"{base}abort-container")
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    op_abort = AbortOp(
        identifier="abort-op",
        participant_container=container.identifier,
        missing_iri=f"{base}missing-object",
    )
    op_dependent = NoOp(identifier="dependent-op", required_precedents=["abort-op"])

    sim = Simulation([op_abort, op_dependent])
    sim.run()

    assert any(e.event_type == "OPERATION_ABORT" and e.operation_id == "abort-op" for e in sim.history_log)
    assert any(e.event_type == "OPERATION_END" and e.operation_id == "dependent-op" for e in sim.history_log)


class LongOp(Operation):
    temporal_cost: float = Field(default=1.0)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_precedent_interrupt_does_not_block_dependents():
    op_long = LongOp(identifier="long-op", temporal_cost=1.0)
    op_dependent = NoOp(identifier="dependent-interrupt", required_precedents=["long-op"])

    sim = Simulation([op_long, op_dependent])
    proc = sim.operation_registry["long-op"]

    def _interrupt(env: simpy.Environment) -> simpy.events.Event:
        yield env.timeout(0.1)
        proc.simpy_process.interrupt("boom")

    sim.env.process(_interrupt(sim.env))
    sim.run()

    assert any(e.event_type == "OPERATION_INTERRUPT" and e.operation_id == "long-op" for e in sim.history_log)
    assert any(e.event_type == "OPERATION_END" and e.operation_id == "dependent-interrupt" for e in sim.history_log)
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_precedents_semantics.py ###
