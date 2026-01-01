from __future__ import annotations

from typing import List

from loguru import logger
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import MaterialContainer, PortionOfMaterial
from libsyn_tools.sim.knowledge_graph import Is_directly_contained_by, LabObject
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.unitary_edit import (
    UnitaryEdit,
    AddObjectProperty,
    RemoveObjectProperty,
    Annihilate,
    Create,
)


class MixInContainer(Operation):
    """
    Mix all directly contained POMs inside a container into one.
    """

    participant_container: StrOrSelector

    def get_operation_effects(self) -> List[UnitaryEdit]:
        container = KnowledgeGraph.get_object_from_lookup(self.participant_container)
        if not isinstance(container, MaterialContainer):
            raise RuntimeError("MixInContainer: container must be a MaterialContainer")

        poms = LabObject.get_directly_contained_individuals(container, PortionOfMaterial, only_present=True)
        if len(poms) <= 1:
            logger.info(f"MixInContainer: {container.identifier} has <= 1 POM; no mix needed")
            return []

        mixed_pom = PortionOfMaterial()
        for pom in poms:
            for chemical in pom.get_ingredients():
                mixed_pom.add_chemical(chemical)

        prop_iri = Is_directly_contained_by.predicate_iri
        edits: list[UnitaryEdit] = []
        for pom in poms:
            edits.extend(
                [
                    RemoveObjectProperty(
                        instance_1_iri=pom.identifier,
                        instance_2_iri=container.identifier,
                        property_iri=prop_iri,
                    ),
                    Annihilate(instance_1_iri=pom.identifier),
                ]
            )
        edits.extend(
            [
                Create(instance_1_iri=mixed_pom.identifier),
                AddObjectProperty(
                    instance_1_iri=mixed_pom.identifier,
                    instance_2_iri=container.identifier,
                    property_iri=prop_iri,
                ),
            ]
        )
        return edits
