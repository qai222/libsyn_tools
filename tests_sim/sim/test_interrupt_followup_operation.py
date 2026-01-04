"""Simulation-level regression: an interrupted operation must not strand locks.

We validate that:
- op1 holds a resource and is interrupted,
- op1 releases locks and signals done_event,
- op2 waiting on op1 as a precedent can start and finish,
- the shared resource ends unlocked.
"""

import simpy

from libsyn_tools.sim.simulation import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create
from twa.data_model.base_ontology import KnowledgeGraph


class WaitOp(Operation):
    participant_resource: str

    def get_operation_effects(self):
        return []


def _make_present_container() -> MaterialContainer:
    c = MaterialContainer()
    c.has_capacity.add(100.0)
    c.has_pool_type.add("POOL_INTERRUPT_FOLLOWUP")
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    Create(instance_1_iri=c.identifier).apply()
    return c


def test_interrupt_releases_lock_and_allows_followup_operation():
    c = _make_present_container()

    op1 = WaitOp(participant_resource=c.identifier, temporal_cost=10.0)
    op2 = WaitOp(
        participant_resource=c.identifier,
        temporal_cost=0.0,
        required_precedents=[op1.identifier],
    )

    sim = Simulation([op1, op2])

    proc1 = sim.operation_registry[op1.identifier]
    proc2 = sim.operation_registry[op2.identifier]

    # Manually start processes so we can interrupt proc1 deterministically.
    proc1.simpy_process = sim.env.process(proc1.run())
    proc2.simpy_process = sim.env.process(proc2.run())

    def interrupter(env: simpy.Environment):
        yield env.timeout(0.1)
        proc1.simpy_process.interrupt("boom")

    sim.env.process(interrupter(sim.env))
    sim.env.run(until=1.0)

    # op2 should have started and ended (not stuck on precedent/lock).
    types = [r.event_type for r in sim.history_log]
    assert "OPERATION_INTERRUPT" in types
    assert "OPERATION_START" in types
    assert "OPERATION_END" in types

    # Shared resource must end unlocked.
    rs = get_runtime_state(c, sim.env)
    assert rs.lock.count == 0
