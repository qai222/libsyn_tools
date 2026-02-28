from pathlib import Path

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import (
    Chemical,
    Create,
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
    Simulation,
    logger,
)
from libsyn_tools.sim.operation_preset import TransferMaterialByPortionSize

HERE = Path(__file__).resolve().parent


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


def run(out_dir: str | Path | None = None) -> Simulation:
    beaker_1, beaker_2, pipette_1 = init_world()
    transfer = TransferMaterialByPortionSize(
        participant_source=beaker_1.instance_iri,
        participant_destination=beaker_2.instance_iri,
        participant_device=pipette_1.instance_iri,
        portion_size=0.3
    )
    sim = Simulation(operations=[transfer])
    sim.run()
    output_dir = Path(out_dir) if out_dir else HERE / "_generated" / "simple_transfer"
    output_dir.mkdir(parents=True, exist_ok=True)
    sim.export_instance_history(output_dir / "instance_history.csv")
    KnowledgeGraph.graph().serialize(destination=output_dir / "state.ttl", format="turtle")
    return sim


if __name__ == "__main__":
    run()
