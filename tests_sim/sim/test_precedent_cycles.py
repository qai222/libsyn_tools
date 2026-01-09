# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_precedent_cycles.py ###
from __future__ import annotations

import pytest
from pydantic import Field

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class NoOp(Operation):
    participant_obj: str = Field(default_factory=lambda: LabObject().identifier)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_precedent_cycle_fails_fast():
    op_a = NoOp(identifier="A", required_precedents=["B"])
    op_b = NoOp(identifier="B", required_precedents=["A"])

    with pytest.raises(ValueError, match=r"Precedent cycle detected: A -> B -> A"):
        Simulation([op_a, op_b])
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_precedent_cycles.py ###
