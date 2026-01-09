from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    Is_directly_contained_by,
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
)
from libsyn_tools.sim.operation.unitary_edit import AddObjectProperty, Create
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def test_transfer_removes_containment_link_from_annihilated_pom() -> None:
    src = MaterialContainer()
    dst = MaterialContainer()
    dev = MaterialContainer()

    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=10.0, density=1.0))

    for o in (src, dst, dev, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    op = TransferMaterialByPortionSize(
        identifier="transfer_cleanup",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        portion_size=0.4,
        temporal_cost=0.0,
    )

    sim = Simulation([op])
    sim.run()

    assert pom.is_present == {False}
    assert src not in pom.is_directly_contained_by

    all_src_poms = LabObject.get_directly_contained_individuals(
        src, PortionOfMaterial, only_present=False
    )
    assert pom not in all_src_poms

    src_poms = LabObject.get_directly_contained_individuals(
        src, PortionOfMaterial, only_present=True
    )
    dst_poms = LabObject.get_directly_contained_individuals(
        dst, PortionOfMaterial, only_present=True
    )

    assert len(src_poms) == 1
    assert len(dst_poms) == 1
    assert src_poms[0].volume == pytest.approx(6.0)
    assert dst_poms[0].volume == pytest.approx(4.0)
