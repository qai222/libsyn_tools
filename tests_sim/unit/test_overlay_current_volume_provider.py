from __future__ import annotations

from rdflib import Namespace
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph   import (
    MaterialContainer, PortionOfMaterial, Is_directly_contained_by
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider

LIB = Namespace("https://libsyn-sim/kg/")


def test_current_volume_overlay_emits_volume_triple():
    c = MaterialContainer()
    p = PortionOfMaterial()
    p.add_chemical(Chemical(mass=12.0, density=1.0))  # → 12 mL

    for o in (c, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()

    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=c.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    g = CurrentVolumeOverlayProvider().snapshot()

    triples = list(g.triples((None, LIB.currentVolume, None)))
    assert len(triples) == 1
    s, _, o = triples[0]
    assert str(s).endswith(c.identifier)
    assert float(o) == 12.0
