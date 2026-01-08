from __future__ import annotations

from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.selector import AttributeSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry


class _WaitForTool(Operation):
    """Operation that blocks in pre_act waiting for a pooled tool."""

    participant_tool: StrOrSelector = Field(...)

    def get_operation_effects(self):
        return []


def test_interrupt_while_waiting_in_pre_act_does_not_strand_locks():
    """If an operation is interrupted while waiting for a selector in pre_act,
    the simulation should:
      - not crash
      - cancel the in-flight pre_act process
      - release any locks acquired so far
      - leave future pooled objects available
    """

    # Tool exists but is not present yet => it won't be put in the pool store.
    tool = MaterialContainer()
    tool.has_pool_type.add("TOOLS")
    tool.is_present = {False}
    KnowledgeGraph.get_object_from_lookup(tool.identifier)

    def _any(_obj) -> bool:
        return True

    op = _WaitForTool(
        identifier="wait",
        participant_tool=AttributeSelector(pool_type="TOOLS", predicate=_any),
        temporal_cost=0.0,
    )

    sim = Simulation([op])
    env = sim.env
    proc = sim.operation_registry[op.identifier]
    proc.simpy_process = env.process(proc.run())

    # Interrupt while pre_act is blocked on the empty pool store.
    def _interrupt_later():
        yield env.timeout(0.1)
        proc.simpy_process.interrupt("cancel")

    env.process(_interrupt_later())

    # After the interrupt, make the tool available. If the pre_act child process
    # was not cancelled, it would wake up, acquire the tool, and strand its lock.
    def _make_tool_available():
        yield env.timeout(0.2)
        tool.is_present = {True}
        FilterStoreRegistry.put_obj_into_filter_store(tool, env)

    env.process(_make_tool_available())

    env.run()

    # Tool remains available and unlocked.
    rs = get_runtime_state(tool, env)
    assert rs.lock.count == 0

    store = FilterStoreRegistry.get_filter_store("TOOLS", env)
    assert tool in store.items

    # The operation finished its interrupt path.
    assert proc.done_event.triggered
    assert any(r.event_type == "OPERATION_INTERRUPT" for r in sim.history_log)
