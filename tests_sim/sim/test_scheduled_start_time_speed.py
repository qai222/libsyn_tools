# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_scheduled_start_time_speed.py ###
from __future__ import annotations

import math
import pytest

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import TimerSpawner


class ScheduledOp(Operation):
    """No participants, no edits; used to validate scheduled start timing."""

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_scheduled_start_time_uses_base_time_with_speed_factor():
    speed_factor = 2.0
    scheduled_base = 5.0
    spawn_offset_base = 2.0

    sim = Simulation.compile_actions(simulation_speed_factor=speed_factor)

    def _factory(_sim: Simulation) -> Operation:
        return ScheduledOp(
            identifier="scheduled-op",
            scheduled_start_time=scheduled_base,
        )

    TimerSpawner(
        op_factory=_factory,
        interval=100.0,
        start_offset=spawn_offset_base,
    ).attach(sim)

    sim.run(until=10.1)

    start_time = next(
        r.timestamp
        for r in sim.history_log
        if r.operation_id == "scheduled-op" and r.event_type == "OPERATION_START"
    )

    expected_start_sim = scheduled_base * speed_factor
    assert math.isclose(start_time, expected_start_sim, rel_tol=0, abs_tol=1e-6)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), -1.0])
def test_scheduled_start_time_invalid_init_rejected(value: float) -> None:
    with pytest.raises(ValueError, match="scheduled_start_time"):
        ScheduledOp(scheduled_start_time=value)
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_scheduled_start_time_speed.py ###
