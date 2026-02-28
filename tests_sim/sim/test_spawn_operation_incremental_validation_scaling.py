from __future__ import annotations

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_spawn_operation_uses_incremental_validation(monkeypatch) -> None:
    sim = Simulation([])

    def _must_not_run(*_args, **_kwargs):
        raise AssertionError("spawn_operation should not call full _validate_precedents")

    monkeypatch.setattr(Simulation, "_validate_precedents", staticmethod(_must_not_run))

    previous_id: str | None = None
    total = 300
    for idx in range(total):
        op_id = f"spawn-incremental-{idx}"
        precedents = [previous_id] if previous_id is not None else []
        sim.spawn_operation(
            _NoOp(identifier=op_id, required_precedents=precedents),
            start_immediately=True,
        )
        previous_id = op_id

    assert len(sim.operation_registry) == total
