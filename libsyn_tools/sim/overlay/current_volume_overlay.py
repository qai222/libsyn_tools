from __future__ import annotations

"""
Current-volume overlay provider.

This emits ephemeral triples:
    <container> lib:currentVolume "double"

Notes
-----
• Only present MaterialContainers are reported.
• Lives entirely in the overlay; does not mutate the base KG.
"""

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import XSD
from libsyn_tools.sim.knowledge_graph.physical_entities import MaterialContainer

LIB = Namespace("https://libsyn-sim/kg/")


class CurrentVolumeOverlayProvider:
    def snapshot(self) -> Graph:
        g = Graph()
        for c in MaterialContainer.object_lookup.values():
            if c.is_present != {True}:
                continue
            vol = c.directly_contained_pom_volume
            g.add((URIRef(c.instance_iri), LIB.currentVolume, Literal(vol, datatype=XSD.double)))
        return g
