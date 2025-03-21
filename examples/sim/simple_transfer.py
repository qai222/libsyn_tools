import os.path

from libsyn_tools.sim import *

"""
a simple transfer action between two containers
"""


def init_world():
    # lab objects
    beaker_1 = LabObject(identifier="beaker_1")
    beaker_2 = LabObject(identifier="beaker_2")
    pipette_1 = LabObject(identifier="pipette_3")

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
