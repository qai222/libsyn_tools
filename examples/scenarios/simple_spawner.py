from libsyn_tools.sim import Simulation, UnitaryEdit, Create, LabObject
from libsyn_tools.sim.spawner import TimerSpawner, Operation

class CalibrateBalance(Operation):
    temporal_cost: float = 10
    participant_balance: str = "balance"
    def get_operation_effects(self) -> list[UnitaryEdit]: return []

def init_world():
    balance = LabObject(identifier="balance")
    Create(instance_1_iri=balance.instance_iri).apply()
init_world()
sim = Simulation.compile_actions()           # empty recipe for test

TimerSpawner(
    op_factory=lambda s: CalibrateBalance(),
    interval=12*3600,                        # float is fine
).attach(sim)

sim.run(until=24*3600)
