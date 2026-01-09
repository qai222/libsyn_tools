# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_remove_object_property_idempotent.py ###
from __future__ import annotations

from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject, Is_directly_contained_by
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import Create, RemoveObjectProperty, UnitaryEdit


class RemoveAbsentRelation(Operation):
    participant_src: str = Field(...)
    participant_dst: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            RemoveObjectProperty(
                instance_1_iri=self.participant_src,
                instance_2_iri=self.participant_dst,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


def test_remove_object_property_is_idempotent():
    src = LabObject()
    dst = LabObject()
    for obj in (src, dst):
        obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    op = RemoveAbsentRelation(participant_src=src.identifier, participant_dst=dst.identifier)
    sim = Simulation([op])
    sim.run()

    types = [r.event_type for r in sim.history_log if r.operation_id == op.identifier]
    assert "OPERATION_END" in types
    assert "OPERATION_ABORT" not in types
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_remove_object_property_idempotent.py ###
