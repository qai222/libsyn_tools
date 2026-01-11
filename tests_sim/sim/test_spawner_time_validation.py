from __future__ import annotations

import pytest

from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.spawner import KGInspectorSpawner, TimerSpawner, ValidationAuditSpawner


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_timer_spawner_interval_invalid(value: float) -> None:
    with pytest.raises(ValueError, match="interval"):
        TimerSpawner(op_factory=lambda s: _NoOp(), interval=value)


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_timer_spawner_start_offset_invalid(value: float) -> None:
    with pytest.raises(ValueError, match="start_offset"):
        TimerSpawner(op_factory=lambda s: _NoOp(), interval=1.0, start_offset=value)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_kg_inspector_interval_invalid(value: float) -> None:
    with pytest.raises(ValueError, match="inspect_interval"):
        KGInspectorSpawner(shape_dispatch={}, inspect_interval=value)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_validation_audit_interval_invalid(value: float) -> None:
    with pytest.raises(ValueError, match="inspect_interval"):
        ValidationAuditSpawner(inspect_interval=value)
