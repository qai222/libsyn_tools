from __future__ import annotations

from pydantic import Field
from rdflib import Namespace
from rdflib.namespace import XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.selector import KgQuerySelector
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty, UnitaryEdit


class SelectTargetByKg(Operation):
    participant_target: StrOrSelector
    name: str = Field("select_target")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_kg_query_selector_uses_overlay_current_volume():
    lib = Namespace("https://libsyn-sim/kg/")
    base = "https://libsyn-sim/kg/"
    pool = "VIAL"

    container_ok = MaterialContainer(identifier=f"{base}with_pom")
    container_ok.has_pool_type.add(pool)
    container_empty = MaterialContainer(identifier=f"{base}empty")
    container_empty.has_pool_type.add(pool)

    pom = PortionOfMaterial(identifier=f"{base}pom")
    pom.add_chemical(Chemical(mass=1.0, density=1.0))

    for obj in (container_ok, container_empty, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container_ok.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    sparql = f"""
    PREFIX lib: <{lib}>
    PREFIX xsd: <{XSD}>
    SELECT ?s WHERE {{
      ?s lib:currentVolume ?v .
      FILTER(xsd:double(?v) > 0)
    }}
    """

    op = SelectTargetByKg(
        participant_target=KgQuerySelector(pool_type=pool, sparql=sparql),
        temporal_cost=0.0,
    )
    sim = Simulation([op])
    sim.run()

    assert op.resolved_resources["target"] == container_ok.identifier
