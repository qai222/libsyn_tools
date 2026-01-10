from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject
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


def test_transfer_volume_missing_source_raises() -> None:
    op = TransferMaterialByVolume(
        participant_source="missing-source-volume",
        participant_destination="dst",
        participant_device="dev",
        transfer_volume=1.0,
    )
    with pytest.raises(RuntimeError, match="missing-source-volume"):
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
