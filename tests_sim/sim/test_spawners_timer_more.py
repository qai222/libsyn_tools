# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_spawners_timer_more.py ###
from __future__ import annotations

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import TimerSpawner


class Ping(Operation):
    """No participants, no edits; just logs lifecycle."""

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _count_end_any(sim: Simulation) -> int:
    """Count all OPERATION_END events in this run (timer-only tests have only Ping)."""
    return sum(1 for r in sim.history_log if r.event_type == "OPERATION_END")


def test_timer_spawns_with_start_offset():
    sim = Simulation.compile_actions()
    t = TimerSpawner(op_factory=lambda s: Ping(), interval=5.0, start_offset=7.0)
    t.attach(sim)

    # Spawns at 7,12,17 within 20
    sim.run(until=20.0)
    assert _count_end_any(sim) == 3


def test_timer_cancel_stops_future_spawns():
    sim = Simulation.compile_actions()
    t = TimerSpawner(op_factory=lambda s: Ping(), interval=5.0)
    t.attach(sim)

    # Cancel after first spawn at t=0
    def canceller():
        yield sim.env.timeout(0.1)
        t.cancel()

    sim.env.process(canceller())

    sim.run(until=20.0)
    # Only the initial spawn should have run
    assert _count_end_any(sim) == 1


def test_timer_determinism_same_seed_same_log(tmp_path):
    # Run 1
    sim1 = Simulation.compile_actions()
    TimerSpawner(op_factory=lambda s: Ping(), interval=4.0).attach(sim1)
    sim1.run(until=16.0)
    log1 = [(r.event_type, r.timestamp) for r in sim1.history_log]  # ignore random op ids

    # Run 2 (fresh world, identical setup)
    sim2 = Simulation.compile_actions()
    TimerSpawner(op_factory=lambda s: Ping(), interval=4.0).attach(sim2)
    sim2.run(until=16.0)
    log2 = [(r.event_type, r.timestamp) for r in sim2.history_log]

    assert log1 == log2


def test_timer_factory_error_does_not_crash():
    sim = Simulation.compile_actions()
    calls = {"count": 0}

    def _factory(_sim: Simulation) -> Ping:
        calls["count"] += 1
        if calls["count"] == 1:
            raise RuntimeError("boom")
        return Ping()

    TimerSpawner(op_factory=_factory, interval=1.0).attach(sim)

    sim.run(until=3.1)

    assert _count_end_any(sim) >= 2
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_spawners_timer_more.py ###
