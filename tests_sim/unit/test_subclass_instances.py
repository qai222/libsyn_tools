from __future__ import annotations

from rdflib import Namespace
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph import (
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
    canonical_iri,
)
from libsyn_tools.sim.operation.unitary_edit import AddObjectProperty, Create
from libsyn_tools.sim.overlay.chemistry_overlay import ChemistryOverlayProvider
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider


class MyContainer(MaterialContainer):
    pass


class MyPOM(PortionOfMaterial):
    pass


def test_subclass_instances_appear_in_overlays_and_queries() -> None:
    container = MyContainer()
    container.has_pool_type.add("POOL_SUBCLASS_CONTAINER")
    pom = MyPOM()
    pom.add_chemical(Chemical(mass=3.0, density=1.0))

    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    current_volume_graph = CurrentVolumeOverlayProvider().snapshot()
    lib = Namespace("https://libsyn-sim/kg/")
    container_iri = canonical_iri(container.identifier)
    assert (container_iri, lib.currentVolume, None) in current_volume_graph

    chemistry_graph = ChemistryOverlayProvider().snapshot()
    pom_iri = canonical_iri(pom.identifier)
    assert (pom_iri, lib.hasIngredient, None) in chemistry_graph

    assert container.directly_contained_pom_volume == 3.0
    contained = LabObject.get_directly_contained_individuals(container, PortionOfMaterial, only_present=True)
    assert pom in contained
    subclass_contained = LabObject.get_directly_contained_individuals(container, MyPOM, only_present=True)
    assert pom in subclass_contained
