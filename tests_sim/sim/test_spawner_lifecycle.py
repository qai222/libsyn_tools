from __future__ import annotations

import pytest

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import KGInspectorSpawner, TimerSpawner


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_periodic_spawner_requires_until() -> None:
    sim = Simulation.compile_actions()
    KGInspectorSpawner(shape_dispatch={}, inspect_interval=1.0).attach(sim)

    with pytest.raises(ValueError, match="periodic spawner"):
        sim.run()


def test_spawn_operation_start_immediately_runs_once() -> None:
    sim = Simulation.compile_actions()
    op = _NoOp(identifier="once")

    sim.spawn_operation(op, start_immediately=True)
    sim.run(until=0.1)

    starts = [r for r in sim.history_log if r.event_type == "OPERATION_START" and r.operation_id == "once"]
    assert len(starts) == 1


def test_detach_timer_spawner_stops_future_spawns() -> None:
    sim = Simulation.compile_actions()

    spawner = TimerSpawner(op_factory=lambda s: _NoOp(), interval=1.0)
    spawner.attach(sim)

    def _detach_after(env):
        yield env.timeout(2.5)
        spawner.detach()

    sim.env.process(_detach_after(sim.env))
    sim.run(until=5.0)

    end_count = sum(1 for r in sim.history_log if r.event_type == "OPERATION_END")
    assert end_count == 3


def test_timer_spawner_factory_error_does_not_crash() -> None:
    sim = Simulation.compile_actions()

    def _bad_factory(_sim: Simulation) -> Operation:
        raise RuntimeError("boom")

    TimerSpawner(op_factory=_bad_factory, interval=1.0).attach(sim)

    sim.run(until=2.5)
    assert len(sim.history_log) == 0
