# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_pre_act_interrupt.py ###
from __future__ import annotations

import simpy
import pytest
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector, _OpState
from libsyn_tools.sim.operation.selector import AttributeSelector


class _WaitForTool(Operation):
    """Operation that blocks in pre_act waiting for a pooled tool."""

    participant_tool: StrOrSelector = Field(...)

    def get_operation_effects(self):
        return []


def test_pre_act_interrupt_propagates_and_does_not_run() -> None:
    env = simpy.Environment()

    tool = MaterialContainer()
    tool.has_pool_type.add("POOL_PREACT_DIRECT")
    tool.is_present = {False}
    KnowledgeGraph.get_object_from_lookup(tool.identifier)

    op = _WaitForTool(
        identifier="wait-direct",
        participant_tool=AttributeSelector(pool_type="POOL_PREACT_DIRECT", predicate=lambda _obj: True),
        temporal_cost=0.0,
    )

    pre_act_proc = op.pre_act(env)

    def _interrupt_later():
        yield env.timeout(0.1)
        pre_act_proc.interrupt("direct-cancel")

    env.process(_interrupt_later())

    with pytest.raises(simpy.Interrupt):
        env.run()

    assert op.sim_state is _OpState.PREPARED
    assert op.operation_effects == []
    assert op.resolved_resources == {}
    assert op.resources == []
    assert op.locks == []

    op.cleanup(env)
    assert op.sim_state is _OpState.FINISHED
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_pre_act_interrupt.py ###
