from pathlib import Path

from loguru import logger

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import LabObject, KnowledgeGraph, PortionOfMaterial, Simulation, \
    Create, MaterialContainer
from libsyn_tools.sim.operation_preset import TransferMaterialByPortionSize

"""
This example is a simulation of multiple transfer actions with the following precedence relations

transfer_1 -- transfer_2 -- transfer 3

transfer_4 -- transfer_5 -- transfer 6
transfer_7 -- |
"""
HERE = Path(__file__).resolve().parent


def init_world():
    """
    Create a few lab objects and water portion, all set to 'present',
    by applying CREATE edits to the global KnowledgeGraph.
    """
    logger.info("Initializing the world.")
    beaker_1 = MaterialContainer(identifier="beaker_1")
    beaker_2 = MaterialContainer(identifier="beaker_2")
    beaker_3 = MaterialContainer(identifier="beaker_3")
    beaker_4 = MaterialContainer(identifier="beaker_4")
    pipette_1 = LabObject(identifier="pipette_1")
    pipette_2 = LabObject(identifier="pipette_2")

    # Create a chemical portion for water
    for obj in [beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2]:
        Create(instance_1_iri=obj.instance_iri).apply()

    for beaker in [beaker_1, beaker_2, beaker_3, beaker_4]:
        water = Chemical.make_up_from_smiles("O")  # e.g., H2O
        water.quantify_by_moles(0.2)
        water_pom = PortionOfMaterial()
        water_pom.add_chemical(water)
        water_pom.is_directly_contained_by.add(beaker)
        Create(instance_1_iri=water_pom.instance_iri).apply()

    return beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2


def get_transfer(index, precedents, source, target, pipette, portion_size, time_cost):
    return TransferMaterialByPortionSize(
        identifier=f"transfer_{index}",
        participant_source=source.instance_iri,
        participant_destination=target.instance_iri,
        participant_device=pipette.instance_iri,
        portion_size=portion_size,
        temporal_cost=time_cost,  # Takes 3 time units
        scheduled_start_time=0.0,  # Earliest start
        required_precedents=precedents,
    )


def build_transfers(beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2):
    transfers = {
        1: get_transfer(1, [], beaker_1, beaker_2, pipette_1, 0.35, 3.0),
        2: get_transfer(2, ["transfer_1"], beaker_2, beaker_3, pipette_1, 0.30, 4.0),
        3: get_transfer(3, ["transfer_2"], beaker_3, beaker_4, pipette_1, 0.25, 2.0),
        4: get_transfer(4, ["transfer_1"], beaker_1, beaker_3, pipette_2, 0.20, 5.0),
        7: get_transfer(7, ["transfer_1"], beaker_2, beaker_4, pipette_2, 0.15, 4.0),
        5: get_transfer(5, ["transfer_4", "transfer_7"], beaker_3, beaker_4, pipette_2, 0.10, 3.0),
        6: get_transfer(6, ["transfer_5"], beaker_4, beaker_1, pipette_1, 0.10, 2.0),
    }
    return [transfers[i] for i in [1, 2, 3, 4, 7, 5, 6]]


def run(out_dir: str | Path | None = None) -> Simulation:
    # 1) Initialize the KnowledgeGraph world
    beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2 = init_world()
    sim = Simulation(
        operations=build_transfers(
            beaker_1, beaker_2, beaker_3, beaker_4, pipette_1, pipette_2
        )
    )

    sim.run()

    output_dir = Path(out_dir) if out_dir else HERE / "_generated" / "dag_transfer"
    output_dir.mkdir(parents=True, exist_ok=True)
    state_file = output_dir / "state.ttl"
    KnowledgeGraph.graph().serialize(destination=state_file, format="turtle")
    logger.info(f"Exported updated knowledge graph to {state_file}")
    sim.export_event_log(output_dir / "event_log.csv")
    return sim


if __name__ == "__main__":
    run()
