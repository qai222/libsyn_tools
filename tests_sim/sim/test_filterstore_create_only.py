# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_filterstore_create_only.py ###
from __future__ import annotations

from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, Create


class CreateOnly(Operation):
    participant_obj: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [Create(instance_1_iri=self.participant_obj)]


def test_filterstore_contains_object_after_create():
    o = LabObject()
    o.has_pool_type.add("VIAL")
    KnowledgeGraph.get_object_from_lookup(o.identifier)

    op = CreateOnly(identifier="C", participant_obj=o.identifier)
    sim = Simulation([op])
    sim.run()

    store = FilterStoreRegistry.get_filter_store("VIAL", sim.env)
    assert o in store.items
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_filterstore_create_only.py ###
