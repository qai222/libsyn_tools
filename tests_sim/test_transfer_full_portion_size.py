from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
    LabObject,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def test_transfer_portion_size_one_moves_all_and_leaves_no_residual() -> None:
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
        identifier="transfer_all",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )

    sim = Simulation([op])
    sim.run()

    # Source should have no present POMs after a full transfer.
    assert LabObject.get_directly_contained_individuals(src, PortionOfMaterial, only_present=True) == []

    # Destination should have exactly one present POM with the full original volume.
    dst_poms = LabObject.get_directly_contained_individuals(dst, PortionOfMaterial, only_present=True)
    assert len(dst_poms) == 1
    assert dst_poms[0].volume == pytest.approx(pom.volume)