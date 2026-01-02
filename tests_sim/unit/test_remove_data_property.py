# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_remove_data_property.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject, Has_interrupt_events
from libsyn_tools.sim.operation.unitary_edit import AddDataProperty, RemoveDataProperty, Create


def test_remove_data_property_removes_value() -> None:
    obj = LabObject(identifier="OBJ_REMOVE")
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()

    AddDataProperty(
        instance_1_iri=obj.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="DIRTY",
    ).apply()
    assert "DIRTY" in obj.has_interrupt_events

    RemoveDataProperty(
        instance_1_iri=obj.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="DIRTY",
    ).apply()
    assert "DIRTY" not in obj.has_interrupt_events
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_remove_data_property.py ###
