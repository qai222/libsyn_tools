from __future__ import annotations

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from pydantic import PrivateAttr

from libsyn_tools.sim.spawner import Spawner


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class _CountingSpawner(Spawner):
    _count: int = PrivateAttr(default=0)

    def _on_attach(self, sim: Simulation) -> None:
        def _on_end(_proc):
            self._count += 1

        sim.callbacks.on_operation_end.append(_on_end)
        self._on_end = _on_end

    def _detach_safe(self, sim: Simulation) -> None:
        if hasattr(self, "_on_end"):
            try:
                sim.callbacks.on_operation_end.remove(self._on_end)
            except ValueError:
                pass


def test_spawner_detached_between_runs() -> None:
    spawner = _CountingSpawner()
    sim_one = Simulation([_NoOp(identifier="first")])
    spawner.attach(sim_one)
    sim_one.run(until=0.1)

    sim_two = Simulation([_NoOp(identifier="second")])
    spawner.attach(sim_two)
    sim_two.run(until=0.1)

    assert spawner._count == 2
