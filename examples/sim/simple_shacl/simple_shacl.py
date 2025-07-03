import os
from pathlib import Path

from loguru import logger

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import (
    Simulation,
    PortionOfMaterial,
    MaterialContainer,
    LabObject,
    Create,
    KnowledgeGraph
)
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize

HERE = Path(__file__).parent


def init_world():
    rack = LabObject(identifier="rack")

    v1 = MaterialContainer(identifier="material_container_v1")
    v2 = MaterialContainer(identifier="material_container_v2")
    for v in (v1, v2):
        v.has_capacity.add(10)
        v.is_directly_contained_by.add(rack)
    reservoir = MaterialContainer(identifier="reservoir", has_capacity={600, })
    pipette_1 = LabObject(identifier="pipette_1")

    # init materials
    water = Chemical.make_up_from_smiles("O")
    water.quantify_by_volume(60)
    water_pom = PortionOfMaterial()
    water_pom.add_chemical(water)
    water_pom.is_directly_contained_by.add(reservoir)

    # creation
    logger.info("create world")
    for thing in [v1, v2, pipette_1, water_pom, reservoir, rack]:
        Create(instance_1_iri=thing.instance_iri).apply()
    return v1, v2, pipette_1, reservoir


def main():
    v1, v2, pipette_1, reservoir = init_world()
    op1 = TransferMaterialByPortionSize(
        identifier="op1",
        participant_source=reservoir.identifier,
        participant_destination=v1.identifier,
        participant_device=pipette_1.identifier,
        portion_size=0.5,
        temporal_cost=2,
    )

    op2 = TransferMaterialByPortionSize(
        identifier="op2",
        participant_source=reservoir.identifier,
        participant_destination=v2.identifier,
        participant_device=pipette_1.identifier,
        portion_size=0.5,
        required_precedents=[op1.identifier],
        temporal_cost=4
    )

    sim = Simulation(
        [op1, op2],
        simulation_speed_factor=1.0,
        shacl_shapes=HERE / "constraints.ttl"
    )

    sim.run()

    overlay = sim.effect_engine._build_overlay_graph()  # noqa: SLF001
    outfile = HERE / "overlay.ttl"
    overlay.serialize(outfile, format="turtle")
    logger.info(f"Overlay graph exported → {outfile}")
    g = KnowledgeGraph.graph()
    g.serialize(destination=f"{os.path.basename(__file__)[:-3]}.ttl", format="turtle")


if __name__ == "__main__":
    main()
