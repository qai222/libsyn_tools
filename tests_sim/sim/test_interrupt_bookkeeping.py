from __future__ import annotations

import simpy
from pydantic import Field
import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import AttributeSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class _InterruptNoLockOp(Operation):
    participant_target: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class _InterruptHoldOneOp(Operation):
    participant_a: str = Field(...)
    participant_b: AttributeSelector = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_interrupt_bookkeeping_skips_when_no_locks_held() -> None:
    obj = LabObject()
    obj.is_present = {True}
    obj.has_pool_type.add("POOL_INTERRUPT_NOLOCK")
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()

    op = _InterruptNoLockOp(
        participant_target=obj.identifier,
        scheduled_start_time=10.0,
    )
    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]

    apply_calls: list[object] = []

    def _apply_fail(*_args, **_kwargs) -> None:
        apply_calls.append("called")
        raise RuntimeError("apply should not be called")

    sim.effect_engine.apply = _apply_fail  # type: ignore[assignment]

    proc.request_interrupt("test-no-locks")
    sim.run(until=1.0)

    assert not apply_calls
    assert any(event.event_type == "OPERATION_INTERRUPT" for event in sim.history_log)


def test_interrupt_bookkeeping_only_targets_held_locks() -> None:
    obj_a = LabObject()
    obj_a.is_present = {True}
    obj_a.has_pool_type.add("POOL_INTERRUPT_HELD")
    KnowledgeGraph.get_object_from_lookup(obj_a.identifier)
    Create(instance_1_iri=obj_a.identifier).apply()

    selector_b = AttributeSelector(
        pool_type="zz_POOL_INTERRUPT_WAIT",
        predicate=lambda _candidate: True,
    )
    op = _InterruptHoldOneOp(
        participant_a=obj_a.identifier,
        participant_b=selector_b,
    )
    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]
    get_runtime_state(obj_a, sim.env)

    captured: dict[str, object] = {}

    def _apply_capture(edits, _env, *, locked_iris, **_kwargs) -> None:
        captured["locked_iris"] = set(locked_iris or [])
        captured["edits"] = list(edits)

    sim.effect_engine.apply = _apply_capture  # type: ignore[assignment]

    def _interrupter(env: simpy.Environment):
        while proc._pre_act_process is None:
            yield env.timeout(0)
        for _ in range(100):
            if proc.operation.locks:
                break
            yield env.timeout(0)
        proc.request_interrupt("test-held")

    sim.env.process(_interrupter(sim.env))
    sim.run(until=1.0)

    assert captured.get("locked_iris") == {obj_a.identifier}
    edits = captured.get("edits", [])
    assert edits
    assert all(getattr(edit, "instance_1_iri", None) == obj_a.identifier for edit in edits)
    assert any(event.event_type == "OPERATION_INTERRUPT" for event in sim.history_log)
