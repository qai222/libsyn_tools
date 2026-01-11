# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_chemistry_overlay.py ###
from __future__ import annotations

from pydantic import Field
import hashlib

from rdflib import Namespace, Literal, URIRef
from rdflib.namespace import OWL
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
    canonical_iri,
    identifier_from_iri,
)
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.selector import KgQuerySelector
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty, UnitaryEdit
from libsyn_tools.sim.overlay.chemistry_overlay import ChemistryOverlayProvider
from libsyn_tools.sim.packs.chemistry import install_chemistry_pack


class SelectByChemistry(Operation):
    participant_target: StrOrSelector
    name: str = Field("select_by_chemistry")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _build_container_with_pom(smiles: str, label: str) -> MaterialContainer:
    base = "https://libsyn-sim/kg/"
    container = MaterialContainer(identifier=f"{base}chem_container_{label}")
    container.has_pool_type.add("VIAL")
    pom = PortionOfMaterial(identifier=f"{base}chem_pom_{label}")
    pom.add_chemical(Chemical(smiles=smiles, mass=2.0, density=1.0))

    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return container


def test_chemistry_overlay_emits_smiles_triples() -> None:
    _build_container_with_pom("CCO", "overlay")
    sim = Simulation([])
    install_chemistry_pack(sim)
    graph = sim.effect_engine.build_query_graph()
    lib = Namespace("https://libsyn-sim/kg/")
    assert (None, lib.smiles, Literal("CCO")) in graph


def test_kg_query_selector_can_filter_by_smiles() -> None:
    container = _build_container_with_pom("CCO", "select")
    lib = Namespace("https://libsyn-sim/kg/")

    sparql = f"""
    PREFIX lib: <{lib}>
    SELECT ?s WHERE {{
      ?pom lib:is_directly_contained_by ?s .
      ?pom lib:hasIngredient ?ing .
      ?ing lib:smiles "CCO" .
    }}
    """

    op = SelectByChemistry(
        participant_target=KgQuerySelector(pool_type="VIAL", sparql=sparql),
        temporal_cost=0.0,
    )
    sim = Simulation([op])
    install_chemistry_pack(sim)
    sim.run()

    assert op.resolved_resources["target"] == container.identifier


def test_ingredient_iris_are_stable_and_hashed() -> None:
    pom = PortionOfMaterial(identifier="pom-hash")
    pom.add_chemical(Chemical(smiles="CCO", mass=1.0, density=1.0))
    pom.add_chemical(Chemical(smiles="O", mass=2.0, density=1.0))

    KnowledgeGraph.get_object_from_lookup(pom.identifier)
    Create(instance_1_iri=pom.identifier).apply()

    provider = ChemistryOverlayProvider()
    first_graph = provider.snapshot()
    lib = Namespace("https://libsyn-sim/kg/")

    ingredient_iris_first = {
        obj for _, _, obj in first_graph.triples((canonical_iri(pom.identifier), lib.hasIngredient, None))
    }
    assert ingredient_iris_first

    ingredient_blobs = list(pom.has_ingredient)
    pom.has_ingredient = {ingredient_blobs[1], ingredient_blobs[0]}

    second_graph = provider.snapshot()
    ingredient_iris_second = {
        obj for _, _, obj in second_graph.triples((canonical_iri(pom.identifier), lib.hasIngredient, None))
    }
    assert ingredient_iris_first == ingredient_iris_second

    expected_hash = hashlib.sha256(ingredient_blobs[0].encode("utf-8")).hexdigest()
    expected_ingredient_id = f"ingredient/{identifier_from_iri(pom.identifier)}/{expected_hash}"
    expected_ingredient_iri = canonical_iri(expected_ingredient_id)
    assert expected_ingredient_iri in ingredient_iris_first
    assert (URIRef(pom.identifier), OWL.sameAs, canonical_iri(pom.identifier)) in first_graph
    assert (URIRef(expected_ingredient_id), OWL.sameAs, expected_ingredient_iri) in first_graph
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_chemistry_overlay.py ###
