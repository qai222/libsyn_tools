from __future__ import annotations
import random
from pathlib import Path
from textwrap import dedent
from rdflib import Graph, Namespace
from loguru import logger
from rdflib.namespace import SH
from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import (
    LabObject, MaterialContainer, PortionOfMaterial, KnowledgeGraph,
    Simulation, Create,
)
from libsyn_tools.sim.operation.unitary_edit import AddDataProperty
from libsyn_tools.sim.operation.operation import Operation, _OpState
from libsyn_tools.sim.spawner import KGInspectorSpawner, ProcessInterruptSpawner
from libsyn_tools.sim.knowledge_graph import Has_interrupt_events
from libsyn_tools.sim.operation_preset import TransferMaterialByPortionSize, DrainExcess

LIB_SYN = Namespace("https://libsyn-sim/kg/")
HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# 0. very small helper operations -------------------------------------------
# ---------------------------------------------------------------------------
class LogNote(Operation):
    participant_target: str
    note: str

    def get_operation_effects(self):
        return [
            AddDataProperty(
                instance_1_iri=self.participant_target,
                property_iri=Has_interrupt_events.predicate_iri,
                data_value=f"NOTE:{self.note}@t={self.temporal_cost}",
            )
        ]

class OverflowAlarm(Operation):
    participant_beaker: str
    def get_operation_effects(self):
        logger.warning(f"Overflow detected in {self.participant_beaker}")
        return []

class ServicePipette(Operation):
    participant_device: str
    temporal_cost: float = 5.0
    def get_operation_effects(self):
        logger.info(f"Servicing {self.participant_device}")
        return []                      # no KG effect for demo

# ---------------------------------------------------------------------------
# 1. build the world (unchanged) --------------------------------------------
# ---------------------------------------------------------------------------
def init_world():
    beakers = [MaterialContainer(identifier=f"beaker_{i}") for i in range(1, 5)]
    pipettes = [LabObject(identifier=f"pipette_{i}") for i in (1, 2)]
    for obj in beakers + pipettes:
        Create(instance_1_iri=obj.identifier).apply()

    for b in beakers:
        w = Chemical.make_up_from_smiles("O"); w.quantify_by_volume(10)      # 10 mL
        pom = PortionOfMaterial(); pom.add_chemical(w)
        pom.is_directly_contained_by.add(b); Create(instance_1_iri=pom.identifier).apply()

    waste = MaterialContainer(identifier="waste_1")
    Create(instance_1_iri=waste.identifier).apply()
    return beakers, pipettes, waste

# ---------------------------------------------------------------------------
# 2. SHACL subset for KGInspectorSpawner ------------------------------------
# ---------------------------------------------------------------------------
shacl_ttl = dedent(f"""
@prefix sh:  <http://www.w3.org/ns/shacl#> .
@prefix lib: <{LIB_SYN}> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

lib:BeakerOverflowShape
    a sh:NodeShape ;
    sh:targetClass lib:MaterialContainer ;
    sh:property [
        sh:path lib:currentVolume ;
        sh:maxExclusive 15 ;   # mL
        sh:datatype xsd:double
    ] .
""")
overflow_shapes = Graph().parse(data=shacl_ttl, format="turtle")

# ---------------------------------------------------------------------------
# 3. transfer network identical to previous example -------------------------
# ---------------------------------------------------------------------------
def make_transfers(beakers, pipettes):
    rng = random.Random(42)
    transfers = {}
    for i in range(1, 8):
        src, dst = rng.sample(beakers, 2)
        transfers[i] = TransferMaterialByPortionSize(
            identifier=f"transfer_{i}",
            participant_source=src.identifier,
            participant_destination=dst.identifier,
            participant_device=rng.choice(pipettes).identifier,
            portion_size=rng.uniform(0.2, 0.6),
            temporal_cost=rng.uniform(1, 8),
        )
    # precedence edges
    transfers[2].required_precedents = [transfers[1].identifier]
    transfers[3].required_precedents = [transfers[2].identifier]
    transfers[5].required_precedents = [transfers[4].identifier, transfers[7].identifier]
    transfers[6].required_precedents = [transfers[5].identifier]
    return list(transfers.values())

# ---------------------------------------------------------------------------
# 4. helper that **causes** random interrupts so we can demo the spawner ----
# ---------------------------------------------------------------------------
def random_interrupts(sim: Simulation, probability=0.3):
    env = sim.env
    yield env.timeout(2)                       # let things start
    while True:
        for proc in list(sim.operation_registry.values()):
            if proc.operation.sim_state is _OpState.RUNNING:
                if random.random() < probability:
                    proc.simpy_process.interrupt("maintenance")
        yield env.timeout(1)

# ---------------------------------------------------------------------------
# 5. build & run the simulation with the three spawners ---------------------
# ---------------------------------------------------------------------------
def example_with_spawners():
    rng = random.Random(42)
    beakers, pipettes, waste = init_world()
    sim = Simulation(
        operations=make_transfers(beakers, pipettes),
        shacl_shapes=overflow_shapes,          # needed for overlay volume
    )

    # # TimerSpawner – periodic status note on beaker_1
    # TimerSpawner(
    #     op_factory=lambda _:
    #         LogNote(participant_target=beakers[0].identifier,
    #                 note="heartbeat",
    #                 temporal_cost=0.0),
    #     interval=8.0
    # ).attach(sim)

    # KGInspectorSpawner – overflow check after every op
    KGInspectorSpawner(
        shape_dispatch={
            str(LIB_SYN.BeakerOverflowShape):
                lambda iri: DrainExcess(
                    identifier=f"drain-{iri.rsplit('/', 1)[-1]}",
                    participant_source=iri,
                    participant_destination="waste_1",
                    participant_device=rng.choice(pipettes).identifier,
                    target_volume=14.99,  # small ε below 15 mL
                )
        },
        inspect_interval=0.0
    ).attach(sim)

    # ProcessInterruptSpawner – schedule pipette service
    ProcessInterruptSpawner(
        interrupt_dispatch={
            "maintenance": lambda proc, _:
                ServicePipette(participant_device=proc.operation.participant_device)
        }
    ).attach(sim)

    # helper coroutine that actually throws the interrupts:
    sim.env.process(random_interrupts(sim))

    logger.info("=== Simulation with spawners – start ===")
    sim.run(until=200)
    logger.info("=== Simulation end, exporting artefacts ===")

    output_dir = HERE / "_generated" / "dag_transfer_spawners"
    output_dir.mkdir(parents=True, exist_ok=True)
    KnowledgeGraph.graph().serialize(output_dir / "state.ttl", format="turtle")
    sim.export_event_log(output_dir / "event_log.csv")
    return sim

if __name__ == "__main__":
    example_with_spawners()
