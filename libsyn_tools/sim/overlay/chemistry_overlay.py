from __future__ import annotations

"""
Chemistry overlay provider.

This emits ephemeral triples:
    <pom> lib:hasIngredient <ingredient_node>
    <ingredient_node> lib:smiles "..."
    <ingredient_node> lib:mass "..."
"""

import json

from rdflib import Graph, Namespace, Literal
from rdflib.namespace import XSD

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph import PortionOfMaterial, canonical_iri, identifier_from_iri, SimOntology

LIB = Namespace(SimOntology.base_url)


class ChemistryOverlayProvider:
    def snapshot(self) -> Graph:
        g = Graph()
        for pom in PortionOfMaterial.all_instances():
            if pom.is_present != {True}:
                continue
            pom_iri = canonical_iri(pom.identifier)
            ingredient_blobs = sorted(pom.has_ingredient)
            for idx, blob in enumerate(ingredient_blobs):
                chemical = Chemical(**json.loads(blob))
                ingredient_id = f"ingredient/{identifier_from_iri(pom.identifier)}/{idx}"
                ingredient_iri = canonical_iri(ingredient_id)
                g.add((pom_iri, LIB.hasIngredient, ingredient_iri))
                if chemical.smiles:
                    g.add((ingredient_iri, LIB.smiles, Literal(chemical.smiles)))
                if chemical.mass is not None:
                    g.add((ingredient_iri, LIB.mass, Literal(float(chemical.mass), datatype=XSD.double)))
        return g
