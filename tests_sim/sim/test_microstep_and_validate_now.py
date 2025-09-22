# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_microstep_and_validate_now.py ###
from __future__ import annotations

from pydantic import Field
from rdflib import Graph
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, Has_interrupt_events
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, AddDataProperty


class InstantOp(Operation):
    """Zero-duration operation that writes a small data property."""
    participant: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddDataProperty(
                instance_1_iri=self.participant,
                property_iri=Has_interrupt_events.predicate_iri,
                data_value="instant",
            )
        ]


def test_apply_happens_in_single_microstep_timestamp_equal():
    obj = LabObject()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op = InstantOp(identifier="I", participant=obj.identifier, temporal_cost=0.0)
    sim = Simulation([op])
    sim.run()

    # START and END at the same sim time for zero-duration op
    start_t = next(r.timestamp for r in sim.history_log if r.operation_id == "I" and r.event_type == "OPERATION_START")
    end_t = next(r.timestamp for r in sim.history_log if r.operation_id == "I" and r.event_type == "OPERATION_END")
    assert start_t == end_t


def test_validate_now_fast_path_returns_true_when_no_shapes():
    sim = Simulation([])
    conforms, report, _ = sim.effect_engine.validate_now()
    assert conforms is True
    assert isinstance(report, Graph)
    # When no shapes are configured, report should be empty/near-empty
    assert len(list(report.triples((None, None, None)))) == 0
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_microstep_and_validate_now.py ###
