from __future__ import annotations

import simpy

from rdflib import Graph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.policy import PolicyBundle, PolicyRule
from libsyn_tools.sim.knowledge_graph import MaterialContainer


ttl_no_interrupt = """
@prefix sh: <http://www.w3.org/ns/shacl#> .
@prefix lib: <https://libsyn-sim/kg/> .

lib:NoInterruptEventsShape a sh:NodeShape ;
    sh:targetSubjectsOf lib:has_interrupt_events ;
    sh:sparql [
        a sh:SPARQLConstraint ;
        sh:select "SELECT ?this WHERE { ?this lib:has_interrupt_events ?e . }" ;
    ] .
"""


class LongWaitOp(Operation):
    participant_resource: StrOrSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_interrupt_does_not_strand_resource_locks() -> None:
    # Resource to lock
    c = MaterialContainer()
    c.has_capacity.add(100.0)
    c.is_present = {True}

    op = LongWaitOp(participant_resource=c.identifier, temporal_cost=10.0)

    shapes = Graph()
    shapes.parse(data=ttl_no_interrupt, format="turtle")

    sim = Simulation([op], shacl_shapes=shapes)
    sim.effect_engine.policy = PolicyBundle(
        per_shape={
            "https://libsyn-sim/kg/NoInterruptEventsShape": PolicyRule(severity="hard", disposition="aborted")
        }
    )

    proc = sim.operation_registry[op.identifier]

    # Start the operation process manually so we can schedule an interrupt.
    proc.simpy_process = sim.env.process(proc.run())

    def interrupter(env: simpy.Environment) -> simpy.events.Event:
        yield env.timeout(0.1)
        proc.simpy_process.interrupt("boom")

    sim.env.process(interrupter(sim.env))
    sim.env.run()

    # The operation should have finished and released its lock.
    rs = get_runtime_state(c, sim.env)
    assert rs.lock.count == 0
    assert proc.done_event.triggered

    # Should have recorded an interrupt in history.
    assert any(e.event_type == "OPERATION_INTERRUPT" for e in sim.history_log)


def test_interrupt_bookkeeping_exception_does_not_crash(monkeypatch) -> None:
    c = MaterialContainer()
    c.is_present = {True}

    op = LongWaitOp(participant_resource=c.identifier, temporal_cost=1.0)
    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]

    def boom_apply(*args, **kwargs):
        raise RuntimeError("bookkeeping boom")

    monkeypatch.setattr(sim.effect_engine, "apply", boom_apply)

    def interrupter(env: simpy.Environment) -> simpy.events.Event:
        yield env.timeout(0.1)
        proc.simpy_process.interrupt("boom")

    sim.env.process(interrupter(sim.env))
    sim.run()

    assert proc.done_event.triggered
    assert any(e.event_type == "OPERATION_INTERRUPT" for e in sim.history_log)
