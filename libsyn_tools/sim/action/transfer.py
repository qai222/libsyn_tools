from abc import ABC

from loguru import logger
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, PortionOfMaterial
from .core import Action, UnitaryEdit, UnitaryEditType


# TODO presumptions

class TransferBase(Action, ABC):
    """ transfer materials quantitatively """

    source_iri: str
    """ the iri of the source container """

    destination_iri: str
    """ the iri of the destination container """

    transfer_device_iri: str
    """ the iri of the transfer device """

    def get_resources(self) -> list[str]:
        """ a list of iris of the resources, they are assumed to be `LabObject` instances """
        resources = [self.source_iri, self.destination_iri, self.transfer_device_iri]
        for iri in resources:
            resource = KnowledgeGraph.get_object_from_lookup(iri=iri)
            resources += [part.instance_iri for part in LabObject.get_all_parts(resource)]
        return resources

    def post_act(self):
        for iri in self.resources:
            instance = KnowledgeGraph.get_object_from_lookup(iri=iri)
            logger.debug(instance)
            logger.debug(LabObject.get_directly_contained_individuals(instance, PortionOfMaterial))


class TransferMaterialByPortionSize(TransferBase):
    """ transfer materials quantified by a given proportion """

    portion_size: float
    """ the portion size of the transferred materials from source container """

    def get_action_effects(self) -> list[UnitaryEdit]:
        unitary_edits = []

        source_container = KnowledgeGraph.get_object_from_lookup(self.source_iri)
        destination_container = KnowledgeGraph.get_object_from_lookup(self.destination_iri)
        transfer_device = KnowledgeGraph.get_object_from_lookup(self.transfer_device_iri)
        source_container: LabObject
        destination_container: LabObject
        transfer_device: LabObject

        poms_source = LabObject.get_directly_contained_individuals(source_container, PortionOfMaterial)
        assert len(poms_source), f"transferring from an empty container: '{self.identifier}' from '{source_container}'"

        for pom_source in poms_source:
            # annihilate pom source
            annihilate_pom_source = UnitaryEdit(type=UnitaryEditType.ANNIHILATE,
                                                instance_1_iri=pom_source.instance_iri, )

            # create transfer
            pom_transfer = pom_source.get_portion(portion_size=self.portion_size)
            pom_transfer.is_directly_contained_by.add(transfer_device)
            create_pom_transfer = UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=pom_transfer.instance_iri, )

            # create pom left in source
            pom_source_new = pom_source.get_portion(portion_size=1 - self.portion_size)
            pom_source_new.is_directly_contained_by.add(source_container)
            create_pom_source_new = UnitaryEdit(type=UnitaryEditType.CREATE,
                                                instance_1_iri=pom_source_new.instance_iri, )

            # create pom in destination
            pom_destination = pom_source.get_portion(portion_size=self.portion_size)
            pom_destination.is_directly_contained_by.add(destination_container)
            create_pom_destination = UnitaryEdit(type=UnitaryEditType.CREATE,
                                                 instance_1_iri=pom_destination.instance_iri, )

            # annihilate pom in transfer device
            annihilate_pom_transfer = UnitaryEdit(type=UnitaryEditType.ANNIHILATE,
                                                  instance_1_iri=pom_transfer.instance_iri, )

            unitary_edits += [
                annihilate_pom_source, create_pom_transfer, create_pom_source_new, create_pom_destination,
                annihilate_pom_transfer
            ]
        return unitary_edits


# TODO DRY
class TransferMaterialByVolume(TransferBase):
    """ transfer materials quantified by a given volume (mL) """

    transfer_volume: float
    """ the volume of the transferred materials from source container """

    def get_action_effects(self) -> list[UnitaryEdit]:
        unitary_edits = []

        source_container = KnowledgeGraph.get_object_from_lookup(self.source_iri)
        destination_container = KnowledgeGraph.get_object_from_lookup(self.destination_iri)
        transfer_device = KnowledgeGraph.get_object_from_lookup(self.transfer_device_iri)
        source_container: LabObject
        destination_container: LabObject
        transfer_device: LabObject

        poms_source = LabObject.get_directly_contained_individuals(source_container, PortionOfMaterial)
        assert len(poms_source), f"transferring from an empty container: '{self.identifier}' from '{source_container}'"

        portion_size = self.transfer_volume / source_container.current_volume

        for pom_source in poms_source:
            # annihilate pom source
            annihilate_pom_source = UnitaryEdit(type=UnitaryEditType.ANNIHILATE,
                                                instance_1_iri=pom_source.instance_iri, )

            # create transfer
            pom_transfer = pom_source.get_portion(portion_size=portion_size)
            pom_transfer.is_directly_contained_by.add(transfer_device)
            create_pom_transfer = UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=pom_transfer.instance_iri, )

            # create pom left in source
            pom_source_new = pom_source.get_portion(portion_size=1 - portion_size)
            pom_source_new.is_directly_contained_by.add(source_container)
            create_pom_source_new = UnitaryEdit(type=UnitaryEditType.CREATE,
                                                instance_1_iri=pom_source_new.instance_iri, )

            # create pom in destination
            pom_destination = pom_source.get_portion(portion_size=portion_size)
            pom_destination.is_directly_contained_by.add(destination_container)
            create_pom_destination = UnitaryEdit(type=UnitaryEditType.CREATE,
                                                 instance_1_iri=pom_destination.instance_iri, )

            # annihilate pom in transfer device
            annihilate_pom_transfer = UnitaryEdit(type=UnitaryEditType.ANNIHILATE,
                                                  instance_1_iri=pom_transfer.instance_iri, )

            unitary_edits += [
                annihilate_pom_source, create_pom_transfer, create_pom_source_new, create_pom_destination,
                annihilate_pom_transfer
            ]
        return unitary_edits
