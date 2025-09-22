# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_missing_precedent.py ###
from __future__ import annotations

import pytest
from pydantic import Field

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class NoOp(Operation):
    participant_obj: str = Field(default_factory=lambda: LabObject().identifier)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_missing_precedent_id_raises_keyerror():
    # One op that claims a nonexistent precedent
    op = NoOp(identifier="B", required_precedents=["does-not-exist"])
    sim = Simulation([op])
    with pytest.raises(KeyError):
        sim.run()
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_missing_precedent.py ###
