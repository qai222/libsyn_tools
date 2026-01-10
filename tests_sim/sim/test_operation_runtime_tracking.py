from __future__ import annotations

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class DualRoleNoop(Operation):
    participant_left: str
    participant_right: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_recent_operations_dedupes_duplicate_bindings() -> None:
    container = MaterialContainer(identifier="recent-dedupe")
    container.has_pool_type.add("POOL_RECENT_DEDUPE")
    container.is_present = {True}
    Create(instance_1_iri=container.identifier).apply()

    op = DualRoleNoop(
        participant_left=container.identifier,
        participant_right=container.identifier,
    )
    sim = Simulation([op])
    sim.run()

    rs = get_runtime_state(container, sim.env)
    assert rs.recent_operations.count(op) == 1
