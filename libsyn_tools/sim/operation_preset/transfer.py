from __future__ import annotations

from typing import List

from loguru import logger
from pydantic import Field, field_validator
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import MaterialContainer
from libsyn_tools.sim.knowledge_graph import (
    LabObject,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.unitary_edit import (
    UnitaryEdit,
    Create,
    Annihilate,
    AddObjectProperty,
    RemoveObjectProperty,
)

_volume_eps = 1e-6


class TransferMaterialByPortionSize(Operation):
    # Participant roles ------------------------------------------------
    participant_source: StrOrSelector
    participant_destination: StrOrSelector
    participant_device: StrOrSelector

    # Parameter --------------------------------------------------------
    portion_size: float = Field(
        ...,
        description="Fraction (0 < value ≤ 1) of every POM in the source "
                    "container that will be moved to the destination.",
    )

    @field_validator("portion_size")
    def _check_portion(cls, v: float) -> float:
        if not 0.0 < v <= 1.0:
            raise ValueError("portion_size must be in the range (0, 1].")
        return v

    @staticmethod
    def _get_operation_effects(src_iri, dst_iri, dev_iri, portion_size):
        src = KnowledgeGraph.get_object_from_lookup(src_iri)
        if src is None:
            raise RuntimeError(f"source container {src_iri!r} not found")

        if not isinstance(src, MaterialContainer):
            raise RuntimeError(
                f"source container {src_iri!r} is not a {MaterialContainer.__name__}"
            )

        edits: List[UnitaryEdit] = []
        prop_iri = Is_directly_contained_by.predicate_iri

        poms = LabObject.get_directly_contained_individuals(
            src, PortionOfMaterial, only_present=True
        )

        logger.debug(f"get directly contained pom: {poms}")

        for pom in poms:
            # Split the original POM into a transferred portion and (optionally)
            # a residual portion. Avoid creating/creating+linking a zero-volume
            # residual when portion_size == 1 (or numerically very close).
            transfer_pom = pom.get_portion(portion_size)
            residual_fraction = 1 - portion_size
            residual_pom = None
            if residual_fraction > _volume_eps:
                residual_pom = pom.get_portion(residual_fraction)

            src_instance_iri = getattr(src, "instance_iri", src.identifier)
            unlink_targets = {src.identifier, src_instance_iri}
            link_present = src in pom.is_directly_contained_by or any(
                target in pom.is_directly_contained_by for target in unlink_targets
            )
            if link_present:
                edits.append(
                    RemoveObjectProperty(
                        instance_1_iri=pom.identifier,
                        instance_2_iri=src.identifier,
                        property_iri=prop_iri,
                    )
                )
            edits.append(Annihilate(instance_1_iri=pom.identifier))

            if residual_pom is not None:
                edits.append(Create(instance_1_iri=residual_pom.identifier))
                edits.append(
                    AddObjectProperty(
                        instance_1_iri=residual_pom.identifier,
                        instance_2_iri=src_iri,
                        property_iri=prop_iri,
                    )
                )

            edits.append(Create(instance_1_iri=transfer_pom.identifier))
            edits.append(
                AddObjectProperty(
                    instance_1_iri=transfer_pom.identifier,
                    instance_2_iri=dev_iri,
                    property_iri=prop_iri,
                )
            )

            edits.append(
                RemoveObjectProperty(
                    instance_1_iri=transfer_pom.identifier,
                    instance_2_iri=dev_iri,
                    property_iri=prop_iri,
                )
            )
            edits.append(
                AddObjectProperty(
                    instance_1_iri=transfer_pom.identifier,
                    instance_2_iri=dst_iri,
                    property_iri=prop_iri,
                )
            )
        return edits

    def get_operation_effects(self) -> List[UnitaryEdit]:
        return self._get_operation_effects(
            self.participant_source, self.participant_destination, self.participant_device, self.portion_size
        )

    # just to make IDE happy
    class Config:
        arbitrary_types_allowed = True


class TransferMaterialByVolume(Operation):
    participant_source: StrOrSelector
    participant_destination: StrOrSelector
    participant_device: StrOrSelector

    transfer_volume: float = Field(
        ...,
        description="Total volume (same unit as Chemical.volume) that "
                    "must be moved from *source* to *destination*.",
    )

    @field_validator("transfer_volume")
    def _check_volume(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("transfer_volume must be > 0.")
        return v

    def get_operation_effects(self) -> List[UnitaryEdit]:
        src = KnowledgeGraph.get_object_from_lookup(self.participant_source)
        if src is None:
            raise RuntimeError(f"source container {self.participant_source!r} not found")
        if not isinstance(src, MaterialContainer):
            raise RuntimeError(
                f"source container {self.participant_source!r} is not a {MaterialContainer.__name__}"
            )

        src_v = src.directly_contained_pom_volume

        # Avoid division-by-zero and provide a clear error message.
        if src_v <= _volume_eps:
            raise RuntimeError(
                f"{src.identifier} is empty (no pom volume) but "
                f"{self.transfer_volume:.3g} requested in volume transfer."
            )

        if src_v + _volume_eps < self.transfer_volume:
            raise RuntimeError(
                f"{src.identifier} has insufficient material; contains only {src_v:.3g} but "
                f"{self.transfer_volume:.3g} requested in volume transfer."
            )

        portion_size = self.transfer_volume / src_v
        portion_size = min(portion_size, 1.0)
        return TransferMaterialByPortionSize._get_operation_effects(
            self.participant_source, self.participant_destination, self.participant_device, portion_size
        )

    class Config:
        arbitrary_types_allowed = True
