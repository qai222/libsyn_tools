from __future__ import annotations

import simpy

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class SlowPreActOp(Operation):
    participant_resource: StrOrSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []

    def _pre_act_implementation(self, env: simpy.Environment) -> None:
        yield from super()._pre_act_implementation(env)
        yield env.timeout(1.0)


def test_interrupt_pre_act_directly_aborts() -> None:
    container = MaterialContainer()
    container.has_pool_type.add("POOL_PRE_ACT_INTERRUPT")
    container.is_present = {True}

    op = SlowPreActOp(participant_resource=container.identifier)
    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]
    proc.simpy_process = sim.env.process(proc.run())

    def interrupter(env: simpy.Environment) -> simpy.events.Event:
        while proc._pre_act_process is None:
            yield env.timeout(0)
        yield env.timeout(0.1)
        proc._pre_act_process.interrupt("direct-pre-act")

    sim.env.process(interrupter(sim.env))
    sim.env.run()

    assert any(event.event_type == "OPERATION_INTERRUPT" for event in sim.history_log)
    assert not any(event.event_type == "OPERATION_START" for event in sim.history_log)


def test_pre_act_interrupt_requeues_pool_item() -> None:
    container = MaterialContainer()
    container.has_pool_type.add("POOL_PRE_ACT_REQUEUE")
    container.is_present = {True}

    op = SlowPreActOp(participant_resource=container.identifier)
    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]
    proc.simpy_process = sim.env.process(proc.run())

    def interrupter(env: simpy.Environment) -> simpy.events.Event:
        while proc._pre_act_process is None:
            yield env.timeout(0)
        yield env.timeout(0.2)
        proc._pre_act_process.interrupt("requeue")

    sim.env.process(interrupter(sim.env))
    sim.env.run()

    rs = get_runtime_state(container, sim.env)
    assert rs.lock.count == 0

    store = FilterStoreRegistry.get_filter_store("POOL_PRE_ACT_REQUEUE", sim.env)
    assert container in store.items
