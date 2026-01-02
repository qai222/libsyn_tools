from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def build_world() -> tuple[MaterialContainer, MaterialContainer, MaterialContainer]:
    src = MaterialContainer(identifier="SRC")
    dst = MaterialContainer(identifier="DST")
    pip = MaterialContainer(identifier="PIP")
    pom = PortionOfMaterial(identifier="POM")
    pom.add_chemical(Chemical(mass=2.0, density=1.0))
    pom.is_directly_contained_by.add(src)

    for obj in (src, dst, pip, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return src, dst, pip


def main() -> None:
    src, dst, pip = build_world()
    op = TransferMaterialByPortionSize(
        identifier="transfer-demo",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=pip.identifier,
        portion_size=0.5,
        temporal_cost=0.0,
    )
    sim = Simulation([op])
    sim.run()
    report = sim.build_report()
    print(report.summary)


if __name__ == "__main__":
    main()
