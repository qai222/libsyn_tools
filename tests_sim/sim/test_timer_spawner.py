# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_timer_spawner.py ###
from __future__ import annotations

import simpy
from pydantic import Field

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import TimerSpawner


class Calibrate(Operation):
    """No participants, no edits; just burns sim-time and logs lifecycle events."""
    label: str = Field(default="calib")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_timer_spawns_periodic_ops(env: simpy.Environment):
    sim = Simulation.compile_actions()  # start empty

    TimerSpawner(
        op_factory=lambda s: Calibrate(),
        interval=5.0,
    ).attach(sim)

    # Run for 20 time units: spawns at t=0,5,10,15 ⇒ expect 4 completions
    sim.run(until=20.0)

    # Count END events for Calibrate operations
    end_count = sum(1 for r in sim.history_log if r.event_type == "OPERATION_END")
    assert end_count == 4


def test_timer_spawns_scaled_by_speed_factor():
    sim = Simulation.compile_actions(simulation_speed_factor=2.0)
    TimerSpawner(
        op_factory=lambda s: Calibrate(),
        interval=5.0,
    ).attach(sim)

    # With speed_factor=2.0, interval=5.0 => spawns at t=0,10,20
    sim.run(until=20.1)
    end_count = sum(1 for r in sim.history_log if r.event_type == "OPERATION_END")
    assert end_count == 3
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_timer_spawner.py ###
