from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject, MaterialContainer
from libsyn_tools.sim.operation_preset.transfer import (
    TransferMaterialByPortionSize,
    TransferMaterialByVolume,
)


def test_transfer_portion_size_missing_source_raises() -> None:
    op = TransferMaterialByPortionSize(
        participant_source="missing-source",
        participant_destination="dst",
        participant_device="dev",
        portion_size=0.5,
    )
    with pytest.raises(RuntimeError, match="missing-source"):
        op.get_operation_effects()


def test_transfer_portion_size_none_source_raises_clear_error() -> None:
    with pytest.raises(RuntimeError, match="source container is required"):
        TransferMaterialByPortionSize._get_operation_effects(
            None, "dst", "dev", 0.5
        )


def test_transfer_volume_missing_source_raises() -> None:
    op = TransferMaterialByVolume(
        participant_source="missing-source-volume",
        participant_destination="dst",
        participant_device="dev",
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="missing-source-volume"):
        op.get_operation_effects()


def test_transfer_volume_none_source_raises_clear_error() -> None:
    op = TransferMaterialByVolume.model_construct(
        participant_source=None,
        participant_destination="dst",
        participant_device="dev",
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="source container is required"):
        op.get_operation_effects()


def test_transfer_volume_non_container_source_raises() -> None:
    obj = LabObject()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op = TransferMaterialByVolume(
        participant_source=obj.identifier,
        participant_destination="dst",
        participant_device="dev",
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="MaterialContainer"):
        op.get_operation_effects()


def test_transfer_volume_none_destination_raises_clear_error() -> None:
    src = MaterialContainer()
    dev = MaterialContainer()
    for obj in (src, dev):
        obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op = TransferMaterialByVolume.model_construct(
        participant_source=src.identifier,
        participant_destination=None,
        participant_device=dev.identifier,
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="destination container is required"):
        op.get_operation_effects()


def test_transfer_volume_non_container_destination_raises() -> None:
    src = MaterialContainer()
    dst = LabObject()
    dev = MaterialContainer()
    for obj in (src, dst, dev):
        if isinstance(obj, MaterialContainer):
            obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op = TransferMaterialByVolume(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="destination container"):
        op.get_operation_effects()
