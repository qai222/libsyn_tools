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

from rdflib import Graph, Namespace, Literal, URIRef
from rdflib.namespace import XSD, OWL
from libsyn_tools.sim.knowledge_graph import MaterialContainer, canonical_iri

LIB = Namespace("https://libsyn-sim/kg/")


class CurrentVolumeOverlayProvider:
    def snapshot(self) -> Graph:
        g = Graph()
        for c in MaterialContainer.object_lookup.values():
            if c.is_present != {True}:
                continue
            vol = c.directly_contained_pom_volume
            canon = canonical_iri(c.instance_iri)
            g.add((canon, LIB.currentVolume, Literal(vol, datatype=XSD.double)))

            # Backwards-compat: some legacy SHACL shapes target the raw identifier
            # as a URIRef (e.g., <v1>) rather than the canonical base_url IRI.
            # Emit owl:sameAs so owlrl inference can bridge the alias without
            # duplicating the currentVolume triple.
            raw = URIRef(c.instance_iri)
            if raw != canon:
                g.add((raw, OWL.sameAs, canon))
        return g
