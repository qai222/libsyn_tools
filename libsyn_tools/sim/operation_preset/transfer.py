from __future__ import annotations

from typing import List

from pydantic import Field, field_validator

from libsyn_tools.sim.knowledge_graph.physical_entities import (
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

    def get_operation_effects(self) -> List[UnitaryEdit]:
        src = LabObject.object_lookup[self.participant_source]
        dst = LabObject.object_lookup[self.participant_destination]
        dev = LabObject.object_lookup[self.participant_device]

        edits: List[UnitaryEdit] = []
        prop_iri = Is_directly_contained_by.predicate_iri

        poms = LabObject.get_directly_contained_individuals(
            src, PortionOfMaterial, only_present=True
        )

        src_iri = src.instance_iri
        dst_iri = dst.instance_iri
        dev_iri = dev.instance_iri

        for pom in poms:
            transfer_pom = pom.get_portion(self.portion_size)
            residual_pom = pom.get_portion(1 - self.portion_size)

            edits += [
                Annihilate(instance_1_iri=pom.identifier),

                Create(instance_1_iri=residual_pom.identifier),
                AddObjectProperty(instance_1_iri=residual_pom.identifier, instance_2_iri=src_iri, property_iri=prop_iri),

                Create(instance_1_iri=transfer_pom.identifier),
                AddObjectProperty(instance_1_iri=transfer_pom.identifier, instance_2_iri=dev_iri, property_iri=prop_iri),

                RemoveObjectProperty(instance_1_iri=transfer_pom.identifier, instance_2_iri=dev_iri, property_iri=prop_iri),
                AddObjectProperty(instance_1_iri=transfer_pom.identifier, instance_2_iri=dst_iri, property_iri=prop_iri),
            ]
        return edits

    # just to make IDE happy
    class Config:
        arbitrary_types_allowed = True
