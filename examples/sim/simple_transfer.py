import os.path
from typing import Any

from loguru import logger
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph.action import UnitaryEdit, UnitaryEditType, Action
from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, PortionOfMaterial

"""
a simple transfer action between two containers
"""


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

    def model_post_init(self, __context: Any) -> None:
        """ populate action effects and resources """
        self.action_effects = self.get_action_effects()
        self.resources = self.get_resources()
        for iri in self.resources:
            resource = KnowledgeGraph.get_object_from_lookup(iri=iri)
            assert resource.is_present == {True, }

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


def init_world():
    # lab objects
    beaker_1 = LabObject(identifier="beaker_1")
    beaker_2 = LabObject(identifier="beaker_2")
    pipette_1 = LabObject(identifier="pipette_1")

    # init materials
    water = Chemical.make_up_from_smiles("O")
    water.quantify_by_moles(0.2)
    water_pom = PortionOfMaterial()
    water_pom.add_chemical(water)
    water_pom.is_directly_contained_by.add(beaker_1)

    # creation
    logger.info("create world")
    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=beaker_1.instance_iri, ).apply()
    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=beaker_2.instance_iri, ).apply()
    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=water_pom.instance_iri, ).apply()
    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=pipette_1.instance_iri, ).apply()
    return beaker_1, beaker_2, pipette_1


if __name__ == '__main__':
    beaker_1, beaker_2, pipette_1 = init_world()
    transfer = TransferMaterialByPortionSize(
        source_iri=beaker_1.instance_iri,
        destination_iri=beaker_2.instance_iri,
        transfer_device_iri=pipette_1.instance_iri,
        portion_size=0.3
    )
    transfer.execute()
    g = KnowledgeGraph.graph()
    g.serialize(destination=f"{os.path.basename(__file__)[:-3]}.ttl", format="turtle")
