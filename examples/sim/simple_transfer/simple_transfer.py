import os

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import LabObject, Chemical, PortionOfMaterial, logger, Create, Simulation, MaterialContainer
from libsyn_tools.sim.operation_preset import TransferMaterialByPortionSize

"""
a simple transfer action between two containers
"""


def init_world():
    # lab objects
    beaker_1 = MaterialContainer(identifier="beaker_1")
    beaker_2 = MaterialContainer(identifier="beaker_2")
    pipette_1 = LabObject(identifier="pipette_3")

    # init materials
    water = Chemical.make_up_from_smiles("O")
    water.quantify_by_moles(0.2)
    water_pom = PortionOfMaterial()
    water_pom.add_chemical(water)
    water_pom.is_directly_contained_by.add(beaker_1)

    # creation
    logger.info("create world")
    for thing in [beaker_1, beaker_2, pipette_1, water_pom]:
        Create(instance_1_iri=thing.instance_iri).apply()
    return beaker_1, beaker_2, pipette_1


if __name__ == '__main__':
    beaker_1, beaker_2, pipette_1 = init_world()
    transfer = TransferMaterialByPortionSize(
        participant_source=beaker_1.instance_iri,
        participant_destination=beaker_2.instance_iri,
        participant_device=pipette_1.instance_iri,
        portion_size=0.3
    )
    sim = Simulation(operations=[transfer])
    sim.run()
    sim.export_instance_history(f"{os.path.basename(__file__)[:-3]}_instance_history.csv")
    g = KnowledgeGraph.graph()
    g.serialize(destination=f"{os.path.basename(__file__)[:-3]}.ttl", format="turtle")
