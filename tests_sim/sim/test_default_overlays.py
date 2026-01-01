from __future__ import annotations

from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty


def test_simulation_registers_overlays_for_shacl_current_volume():
    lib = Namespace("https://libsyn-sim/kg/")
    ttl = f"""
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:CapacityCheckShape a sh:NodeShape ;
        sh:targetSubjectsOf lib:currentVolume ;
        sh:property [
            sh:path lib:currentVolume ;
            sh:lessThanOrEquals lib:has_capacity ;
        ] .
    """
    shape_graph = Graph().parse(data=ttl, format="turtle")

    base = "https://libsyn-sim/kg/"
    container = MaterialContainer(identifier=f"{base}container")
    container.has_capacity.add(1.0)
    pom = PortionOfMaterial(identifier=f"{base}pom")
    pom.add_chemical(Chemical(mass=1.2, density=1.0))

    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    sim = Simulation([], shacl_shapes=shape_graph)
    conforms, _, _ = sim.effect_engine.validate_now()

    assert conforms is False
