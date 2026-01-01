# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_mix_preset.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
    LabObject,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.mix import MixInContainer

_EPS = 1e-6


def _create_container_with_poms():
    container = MaterialContainer(identifier="C1")
    pom_a = PortionOfMaterial(identifier="P1")
    pom_b = PortionOfMaterial(identifier="P2")
    pom_a.add_chemical(Chemical(mass=2.0, density=1.0))
    pom_b.add_chemical(Chemical(mass=3.0, density=1.0))

    for obj in (container, pom_a, pom_b):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    prop_iri = Is_directly_contained_by.predicate_iri
    for pom in (pom_a, pom_b):
        AddObjectProperty(
            instance_1_iri=pom.identifier,
            instance_2_iri=container.identifier,
            property_iri=prop_iri,
        ).apply()

    return container


def test_mix_in_container_merges_poms() -> None:
    container = _create_container_with_poms()
    volume_before = container.directly_contained_pom_volume

    op = MixInContainer(identifier="mix", participant_container=container.identifier)
    sim = Simulation([op])
    sim.run()

    poms = LabObject.get_directly_contained_individuals(container, PortionOfMaterial, only_present=True)
    assert len(poms) == 1
    volume_after = container.directly_contained_pom_volume
    assert abs(volume_after - volume_before) <= _EPS
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_mix_preset.py ###
