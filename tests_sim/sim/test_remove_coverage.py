# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_remove_coverage.py ###
from __future__ import annotations

import pytest
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    LabObject,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, RemoveObjectProperty


class BadRemove(Operation):
    """
    Attempt to REMOVE an object property where only one endpoint is locked
    and neither endpoint is created in this batch -> should fail coverage.
    """
    participant_src: str = Field(...)
    dst_iri: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            RemoveObjectProperty(
                instance_1_iri=self.participant_src,
                instance_2_iri=self.dst_iri,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


def test_remove_object_property_requires_lock_coverage():
    a = LabObject()
    b = LabObject()
    KnowledgeGraph.get_object_from_lookup(a.identifier)
    KnowledgeGraph.get_object_from_lookup(b.identifier)

    op = BadRemove(participant_src=a.identifier, dst_iri=b.identifier)
    sim = Simulation([op])

    with pytest.raises(RuntimeError, match="Mechanical check failed"):
        sim.run()
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_remove_coverage.py ###
