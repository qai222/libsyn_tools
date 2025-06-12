import random

from loguru import logger

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import LabObject, KnowledgeGraph, PortionOfMaterial, UnitaryEdit, UnitaryEditType, \
    TransferMaterialByPortionSize, ActionSimulation

"""
This example is a simulation of multiple transfer actions with the following precedence relations

transfer_1 -- transfer_2 -- transfer 3

transfer_4 -- transfer_5 -- transfer 6
transfer_7 -- |
"""


def init_world():
    """
    Create a few lab objects and water portion, all set to 'present',
    by applying CREATE edits to the global KnowledgeGraph.
    """
    logger.info("Initializing the world.")
    beaker_1 = LabObject(identifier="beaker_1")
    beaker_2 = LabObject(identifier="beaker_2")
    beaker_3 = LabObject(identifier="beaker_3")
    beaker_4 = LabObject(identifier="beaker_4")
    pipette_1 = LabObject(identifier="pipette_1")
    pipette_2 = LabObject(identifier="pipette_2")

    # Create a chemical portion for water
    for lab_object in [beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2]:
        UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=lab_object.instance_iri).apply()

    for beaker in [beaker_1, beaker_2, beaker_3, beaker_4]:
        water = Chemical.make_up_from_smiles("O")  # e.g., H2O
        water.quantify_by_moles(0.2)
        water_pom = PortionOfMaterial()
        water_pom.add_chemical(water)
        water_pom.is_directly_contained_by.add(beaker)
        UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=water_pom.instance_iri).apply()

    return beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2


def get_transfer(index, precedents, source, target, pipette, portion_size, time_cost):
    return TransferMaterialByPortionSize(
        identifier=f"transfer_{index}",
        source_iri=source.instance_iri,
        destination_iri=target.instance_iri,
        transfer_device_iri=pipette.instance_iri,
        portion_size=portion_size,
        temporal_cost=time_cost,  # Takes 3 time units
        scheduled_start_time=0.0,  # Earliest start
        required_precedents=precedents,
    )


def example_simulation():
    # 1) Initialize the KnowledgeGraph world
    beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2 = init_world()
    beakers = [beaker_1, beaker_2, beaker_3, beaker_4]
    pipettes = [pipette_1, pipette_2]
    random.seed(42)
    transfers = dict()
    for i in range(1, 8):
        source, target = random.sample(beakers, 2)
        portion_size = random.random()
        time_cost = random.uniform(1, 15)
        t = get_transfer(
            i, [],
            source, target,
            random.choice(pipettes), portion_size, time_cost)
        transfers[i] = t

    transfers[2].required_precedents = [transfers[1].identifier]
    transfers[3].required_precedents = [transfers[2].identifier]
    transfers[5].required_precedents = [transfers[4].identifier, transfers[7].identifier]
    transfers[6].required_precedents = [transfers[5].identifier]

    sim = ActionSimulation(actions=list(transfers.values()))

    sim.run()

    out_file = __file__.replace(".py", ".ttl")
    KnowledgeGraph.graph().serialize(destination=out_file, format="turtle")
    logger.info(f"Exported updated knowledge graph to {out_file}")
    sim.export_event_log("action_history.csv")


if __name__ == "__main__":
    example_simulation()
