from __future__ import annotations

import pytest

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_spawn_cycle_length_two_detected() -> None:
    op_a = NoOp(identifier="op-a")
    sim = Simulation([op_a])
    op_b = NoOp(identifier="op-b", required_precedents=["op-a"])
    op_a.required_precedents = ["op-b"]

    with pytest.raises(ValueError, match="Precedent cycle detected"):
        sim.spawn_operation(op_b)


def test_deep_precedent_chain_validates_without_recursion_error() -> None:
    ops: list[Operation] = []
    for idx in range(2000):
        identifier = f"op-{idx}"
        precedents = [f"op-{idx - 1}"] if idx else []
        ops.append(NoOp(identifier=identifier, required_precedents=precedents))

    Simulation(ops)


def test_deferred_precedent_is_rejected() -> None:
    sim = Simulation([])
    deferred = NoOp(identifier="deferred-op")
    sim.spawn_operation(deferred, start_immediately=False)

    dependent = NoOp(identifier="dependent-op", required_precedents=["deferred-op"])
    with pytest.raises(ValueError, match="deferred operation"):
        sim.spawn_operation(dependent)


def test_spawn_same_operation_instance_twice_rejected() -> None:
    sim = Simulation([])
    op = NoOp(identifier="repeat-op")
    sim.spawn_operation(op)

    with pytest.raises(ValueError, match="already exists"):
        sim.spawn_operation(op)


def test_operation_identifier_requires_identifier_form() -> None:
    sim = Simulation([])
    op = NoOp(identifier="https://libsyn-sim/kg/canonical-op")
    with pytest.raises(ValueError, match="identifier form"):
        sim.spawn_operation(op)


def test_interrupt_before_process_start_records_terminal_event() -> None:
    sim = Simulation([])
    op = NoOp(identifier="interrupt-before-start")
    proc = sim.spawn_operation(op, start_immediately=False)

    proc.request_interrupt("deferred-cancel")

    assert any(
        evt.event_type == "OPERATION_INTERRUPT" and evt.operation_id == "interrupt-before-start"
        for evt in sim.history_log
    )
