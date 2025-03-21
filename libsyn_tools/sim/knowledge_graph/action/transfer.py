from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, PortionOfMaterial
from .core import Action, UnitaryEdit, UnitaryEditType


class TransferMaterialByPortionSize(Action):
    """ transfer materials quantified by a given proportion """

    source_iri: str
    """ the iri of the source container """

    destination_iri: str
    """ the iri of the destination container """

    transfer_device_iri: str
    """ the iri of the transfer device """

    portion_size: float
    """ the portion size of the transferred materials from source container """

    def get_resources(self) -> list[str]:
        """ a list of iris of the resources, they are assumed to be `LabObject` instances """
        resources = [self.source_iri, self.destination_iri, self.transfer_device_iri]
        for iri in resources:
            resource = KnowledgeGraph.get_object_from_lookup(iri=iri)
            resources += [part.instance_iri for part in LabObject.get_parts(resource)]
        return resources

    def get_action_effects(self) -> list[UnitaryEdit]:
        unitary_edits = []

        source_container = KnowledgeGraph.get_object_from_lookup(self.source_iri)
        destination_container = KnowledgeGraph.get_object_from_lookup(self.destination_iri)
        source_container: LabObject
        destination_container: LabObject

        poms_source = LabObject.get_directly_contained_individuals(source_container, PortionOfMaterial)
        for pom_source in poms_source:
            pom_transfer = pom_source.get_portion(portion_size=self.portion_size)
            pom_transfer.is_directly_contained_by.add(destination_container)
            pom_remain = pom_source.get_portion(portion_size=1 - self.portion_size)
            pom_remain.is_directly_contained_by.add(source_container)
            create_pom_transfer = UnitaryEdit(
                type=UnitaryEditType.CREATE,
                instance_1_iri=pom_transfer.instance_iri,
            )
            create_pom_remain = UnitaryEdit(
                type=UnitaryEditType.CREATE,
                instance_1_iri=pom_remain.instance_iri,
            )
            annihilate_pom_source = UnitaryEdit(
                type=UnitaryEditType.ANNIHILATE,
                instance_1_iri=pom_source.instance_iri,
            )
            unitary_edits += [create_pom_transfer, create_pom_remain, annihilate_pom_source]

        return unitary_edits
