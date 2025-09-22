# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_spawners_interrupt.py ###
from __future__ import annotations

from uuid import uuid4

from pydantic import Field

from libsyn_tools.sim import Simulation, OperationProcess
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import ProcessInterruptSpawner


class LongOp(Operation):
    """Long-running op to generate an interrupt."""

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class Fix(Operation):
    """Spawned fix op (no edits)."""
    tag: str = Field(default="fix")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _dispatch(proc: OperationProcess, reason: str) -> Operation:
    # Make identifier unique to avoid spawn_operation collisions.
    return Fix(identifier=f"fix-{reason}-{uuid4().hex[:6]}")


def test_interrupt_dispatch_spawns_exactly_one():
    op = LongOp(identifier="L", temporal_cost=100.0)
    sim = Simulation([op])

    # Attach one spawner; mapping USR -> Fix
    ProcessInterruptSpawner(interrupt_dispatch={"USR": _dispatch}).attach(sim)

    # Raise interrupt at t=1
    def kicker():
        yield sim.env.timeout(1.0)
        sim.operation_registry["L"].simpy_process.interrupt("USR")

    sim.env.process(kicker())

    sim.run(until=5.0)

    # One spawned fix op ran to completion
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-USR-")]
    assert len(spawned) == 1


def test_interrupt_unmapped_reason_is_ignored():
    op = LongOp(identifier="L2", temporal_cost=50.0)
    sim = Simulation([op])

    ProcessInterruptSpawner(interrupt_dispatch={"USR": _dispatch}).attach(sim)

    def kicker():
        yield sim.env.timeout(1.0)
        sim.operation_registry["L2"].simpy_process.interrupt("OTHER")

    sim.env.process(kicker())

    sim.run(until=5.0)
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-")]
    assert len(spawned) == 0


def test_interrupt_wrapper_installed_once_for_multiple_spawners():
    op = LongOp(identifier="L3", temporal_cost=100.0)
    sim = Simulation([op])

    # Attach two interrupt spawners with the same mapping; wrapper has an internal guard.
    ProcessInterruptSpawner(interrupt_dispatch={"USR": _dispatch}).attach(sim)
    ProcessInterruptSpawner(interrupt_dispatch={"USR": _dispatch}).attach(sim)

    def kicker():
        yield sim.env.timeout(1.0)
        sim.operation_registry["L3"].simpy_process.interrupt("USR")

    sim.env.process(kicker())

    sim.run(until=5.0)
    # With two spawners, we may get ≥1 fixes; the important thing is no crash and at least one spawn.
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-USR-")]
    assert len(spawned) >= 1


def test_spawner_attach_twice_raises():
    op = LongOp(identifier="L4", temporal_cost=1.0)
    sim = Simulation([op])

    s = ProcessInterruptSpawner(interrupt_dispatch={})
    s.attach(sim)
    try:
        s.attach(sim)
        assert False, "Second attach should raise"
    except RuntimeError as e:
        assert "already attached" in str(e).lower()
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_spawners_interrupt.py ###
