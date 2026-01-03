from __future__ import annotations

from typing import List

from loguru import logger
from pydantic import Field, field_validator
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import MaterialContainer
from libsyn_tools.sim.knowledge_graph import LabObject, PortionOfMaterial, Is_directly_contained_by
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.unitary_edit import (
    UnitaryEdit,
    Create,
    Annihilate,
    AddObjectProperty,
    RemoveObjectProperty,
)

_VOLUME_EPS = 1e-6


class TransferMaterialByPortionSize(Operation):
    """Transfer the same fraction of each POM in a source container to a destination.

    This models transfer as:
        - annihilate the original POM;
        - (optionally) create a residual POM staying in the source;
        - create a transferred POM that moves via the device to the destination.

    Correctness detail:
        When `portion_size` is ~1.0, the residual fraction is ~0.0. We therefore
        skip creating a zero-volume residual POM.
    """

    participant_source: StrOrSelector
    participant_destination: StrOrSelector
    participant_device: StrOrSelector

    portion_size: float = Field(
        ...,
        description=(
            "Fraction (0 < value ≤ 1) of every POM in the source container that will be moved "
            "to the destination."
        ),
    )

    @field_validator("portion_size")
    def _check_portion(cls, v: float) -> float:
        if not 0.0 < v <= 1.0:
            raise ValueError("portion_size must be in the range (0, 1].")
        return v

    @staticmethod
    def _get_operation_effects(src_iri: str, dst_iri: str, dev_iri: str, portion_size: float) -> List[UnitaryEdit]:
        src = KnowledgeGraph.get_object_from_lookup(src_iri)

        if not isinstance(src, MaterialContainer):
            raise RuntimeError(
                f"source container: {src.__class__.__name__} {getattr(src, 'instance_iri', src_iri)} is not a {MaterialContainer.__name__}"
            )

        edits: List[UnitaryEdit] = []
        prop_iri = Is_directly_contained_by.predicate_iri

        poms = LabObject.get_directly_contained_individuals(src, PortionOfMaterial, only_present=True)
        logger.debug(f"get directly contained pom: {poms}")

        for pom in poms:
            edits.append(Annihilate(instance_1_iri=pom.identifier))

            residual_fraction = max(0.0, 1.0 - portion_size)
            if residual_fraction > _VOLUME_EPS:
                residual_pom = pom.get_portion(residual_fraction)
                edits.extend(
                    [
                        Create(instance_1_iri=residual_pom.identifier),
                        AddObjectProperty(
                            instance_1_iri=residual_pom.identifier,
                            instance_2_iri=src_iri,
                            property_iri=prop_iri,
                        ),
                    ]
                )

            transfer_pom = pom.get_portion(portion_size)
            edits.extend(
                [
                    Create(instance_1_iri=transfer_pom.identifier),
                    AddObjectProperty(
                        instance_1_iri=transfer_pom.identifier,
                        instance_2_iri=dev_iri,
                        property_iri=prop_iri,
                    ),
                    RemoveObjectProperty(
                        instance_1_iri=transfer_pom.identifier,
                        instance_2_iri=dev_iri,
                        property_iri=prop_iri,
                    ),
                    AddObjectProperty(
                        instance_1_iri=transfer_pom.identifier,
                        instance_2_iri=dst_iri,
                        property_iri=prop_iri,
                    ),
                ]
            )

        return edits

    def get_operation_effects(self) -> List[UnitaryEdit]:
        return self._get_operation_effects(
            self.participant_source,
            self.participant_destination,
            self.participant_device,
            self.portion_size,
        )

    class Config:
        arbitrary_types_allowed = True


class TransferMaterialByVolume(Operation):
    """Transfer a specified total volume from a source container to a destination.

    The operation computes a `portion_size = transfer_volume / current_source_volume`
    and delegates to `TransferMaterialByPortionSize`.

    Correctness detail:
        If the source currently contains no POM volume, we raise a clear error instead
        of dividing by zero.
    """

    participant_source: StrOrSelector
    participant_destination: StrOrSelector
    participant_device: StrOrSelector

    transfer_volume: float = Field(
        ...,
        description=(
            "Total volume (same unit as Chemical.volume) that must be moved from *source* to *destination*."
        ),
    )

    @field_validator("transfer_volume")
    def _check_volume(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("transfer_volume must be > 0.")
        return v

    def get_operation_effects(self) -> List[UnitaryEdit]:
        src = KnowledgeGraph.get_object_from_lookup(self.participant_source)
        if not isinstance(src, MaterialContainer):
            raise RuntimeError("TransferMaterialByVolume: source must be a MaterialContainer")

        src_v = src.directly_contained_pom_volume
        if src_v <= _VOLUME_EPS:
            raise RuntimeError(f"{src.identifier} contains no POM volume to transfer")

        if src_v + _VOLUME_EPS < self.transfer_volume:
            raise RuntimeError(
                f"{src.identifier} contains only {src_v:.3g} but {self.transfer_volume:.3g} requested in volume transfer."
            )

        portion_size = min(self.transfer_volume / src_v, 1.0)
        return TransferMaterialByPortionSize._get_operation_effects(
            self.participant_source,
            self.participant_destination,
            self.participant_device,
            portion_size,
        )

    class Config:
        arbitrary_types_allowed = True
