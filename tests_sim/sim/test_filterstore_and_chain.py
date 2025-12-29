# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_filterstore_and_chain.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, Create, Annihilate


class CreateObj(Operation):
    participant_obj: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [Create(instance_1_iri=self.participant_obj)]


class AnnihilateObj(Operation):
    participant_obj: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [Annihilate(instance_1_iri=self.participant_obj)]


def test_filterstore_normalization_create_then_annihilate():
    # Object with a pool type should appear in FilterStore after CREATE, and be removed after ANNIHILATE.
    o = LabObject()
    o.has_pool_type.add("VIAL")
    KnowledgeGraph.get_object_from_lookup(o.identifier)

    op1 = CreateObj(identifier="C", participant_obj=o.identifier, temporal_cost=0.0)
    op2 = AnnihilateObj(identifier="D", participant_obj=o.identifier, temporal_cost=0.0, required_precedents=["C"])

    sim = Simulation([op1, op2])
    sim.run()

    store = FilterStoreRegistry.get_filter_store("VIAL", sim.env)
    assert o not in store.items  # annihilated objects are removed from pool
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_filterstore_and_chain.py ###
