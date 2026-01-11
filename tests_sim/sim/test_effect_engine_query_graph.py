from __future__ import annotations

from rdflib import Namespace, URIRef
from rdflib.namespace import RDF
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.env_utils import get_effect_engine
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
    canonical_iri,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty


def test_build_query_graph_includes_current_volume_overlay():
    lib = Namespace("https://libsyn-sim/kg/")
    base = "https://libsyn-sim/kg/"

    container = MaterialContainer(identifier=f"{base}container")
    pom = PortionOfMaterial(identifier=f"{base}pom")
    pom.add_chemical(Chemical(mass=1.0, density=1.0))

    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    data_graph = KnowledgeGraph.graph()
    base_triple = next(
        data_graph.triples((URIRef(container.identifier), RDF.type, None)),
        None,
    )
    assert base_triple is not None

    sim = Simulation([])
    engine = get_effect_engine(sim.env)
    union_graph = engine.build_query_graph()

    overlay_triples = list(
        union_graph.triples(
            (canonical_iri(container.instance_iri), lib.currentVolume, None)
        )
    )
    assert len(overlay_triples) >= 1

    base_triples = list(union_graph.triples(base_triple))
    assert base_triples == [base_triple]
