# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_functional_set_validation.py ###
from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import Create


def test_multi_pool_type_raises_value_error() -> None:
    sim = Simulation([])
    container = MaterialContainer(identifier="pool-type-multi")
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()
    container.is_present = {True}
    container.has_pool_type.update({"POOL_TYPE_MULTI_A", "POOL_TYPE_MULTI_B"})

    with pytest.raises(ValueError, match=r"has_pool_type.*found 2 values"):
        FilterStoreRegistry.put_obj_into_filter_store(container, sim.env)


def test_missing_capacity_raises_value_error() -> None:
    container = MaterialContainer(identifier="missing-capacity")
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    with pytest.raises(ValueError, match=r"has_capacity.*missing-capacity"):
        _ = container.capacity
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_functional_set_validation.py ###
