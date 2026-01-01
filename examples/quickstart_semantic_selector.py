from __future__ import annotations

from pathlib import Path

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


class SelectTarget(Operation):
    participant_target: StrOrSelector
    name: str = Field("select_target")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def run(out_dir: str | Path | None = None):
    base = "https://libsyn-sim/kg/"
    lib = Namespace(base)
    pool = "VIAL"

    container = MaterialContainer(identifier=f"{base}selector-container")
    container.has_pool_type.add(pool)
    pom = PortionOfMaterial(identifier=f"{base}selector-pom")
    pom.add_chemical(Chemical(mass=1.0, density=1.0))

    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
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

    op = SelectTarget(
        participant_target=KgQuerySelector(pool_type=pool, sparql=sparql),
        temporal_cost=0.0,
    )
    sim = Simulation([op])
    if out_dir:
        sim.run_and_report(out_dir)
    else:
        sim.run()
    return sim


if __name__ == "__main__":
    run()
