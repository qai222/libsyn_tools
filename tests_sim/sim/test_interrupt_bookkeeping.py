from __future__ import annotations

from typing import Generator

from pydantic import Field
import simpy

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer, identifier_from_iri
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector, _OpState
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.runtime import get_object_for_resource
from libsyn_tools.sim.operation.unitary_edit import Create


class _InterruptWaiter(Operation):
    participant_tool: StrOrSelector = Field(...)

    def get_operation_effects(self):
        return []


class _DualParticipantWait(Operation):
    participant_a: StrOrSelector = Field(...)
    participant_b: StrOrSelector = Field(...)

    def get_operation_effects(self):
        return []


def test_interrupt_bookkeeping_skips_when_no_locks_held() -> None:
    tool = MaterialContainer()
    tool.has_pool_type.add("POOL_INTERRUPT_NO_LOCK")
    tool.is_present = {False}
    KnowledgeGraph.get_object_from_lookup(tool.identifier)

    def _any(_obj) -> bool:
        return True

    op = _InterruptWaiter(
        identifier="interrupt-no-locks",
        participant_tool=AttributeSelector(pool_type="POOL_INTERRUPT_NO_LOCK", predicate=_any),
        temporal_cost=0.0,
    )
    sim = Simulation([op])
    env = sim.env
    proc = sim.operation_registry[op.identifier]

    called = {"apply": False}

    def _apply_stub(*_args, **_kwargs):
        called["apply"] = True
        raise RuntimeError("apply should not be called")

    sim.effect_engine.apply = _apply_stub

    def _interrupt_later() -> Generator[simpy.events.Event, None, None]:
        yield env.timeout(0.1)
        proc.request_interrupt("no-locks")

    env.process(_interrupt_later())
    sim.run()

    assert proc.done_event.triggered
    assert called["apply"] is False


def test_interrupt_bookkeeping_only_targets_locked_objects() -> None:
    obj_a = MaterialContainer(identifier="https://libsyn-sim/kg/interrupt-a")
    obj_a.has_pool_type.add("POOL_INTERRUPT_A")
    KnowledgeGraph.get_object_from_lookup(obj_a.identifier)
    Create(instance_1_iri=obj_a.identifier).apply()

    obj_b = MaterialContainer(identifier="https://libsyn-sim/kg/interrupt-b")
    obj_b.has_pool_type.add("POOL_INTERRUPT_B")
    KnowledgeGraph.get_object_from_lookup(obj_b.identifier)
    Create(instance_1_iri=obj_b.identifier).apply()

    op = _DualParticipantWait(
        identifier="interrupt-held-lock",
        participant_a=obj_a.identifier,
        participant_b=obj_b.identifier,
        temporal_cost=1.0,
    )

    sim = Simulation([op])
    env = sim.env
    proc = sim.operation_registry[op.identifier]

    captured: dict[str, list[str]] = {"locked_iris": [], "edit_iris": []}

    def _apply_stub(edits, *_args, **_kwargs):
        captured["locked_iris"] = list(_kwargs.get("locked_iris", []))
        captured["edit_iris"] = [edit.instance_1_iri for edit in edits]
        return None

    sim.effect_engine.apply = _apply_stub

    def _interrupt_later() -> Generator[simpy.events.Event, None, None]:
        for _ in range(200):
            if op.sim_state is _OpState.RUNNING:
                break
            if op.sim_state is _OpState.FINISHED:
                return
            yield env.timeout(0.01)
        if op.sim_state is not _OpState.RUNNING:
            return
        for req in list(op.locks):
            try:
                obj = get_object_for_resource(req.resource)
            except Exception:
                continue
            if obj.identifier == obj_b.identifier:
                req.resource.release(req)
                op.locks.remove(req)
                FilterStoreRegistry.put_obj_into_filter_store(obj_b, env)
                break
        yield env.timeout(0.05)
        proc.request_interrupt("held-lock-only")

    env.process(_interrupt_later())
    sim.run()

    assert proc.done_event.triggered
    assert identifier_from_iri(obj_a.identifier) in captured["locked_iris"]
    assert identifier_from_iri(obj_b.identifier) not in captured["locked_iris"]
    assert captured["edit_iris"] == [obj_a.identifier]
    store = FilterStoreRegistry.get_filter_store("POOL_INTERRUPT_A", env)
    assert any(item.identifier == obj_a.identifier for item in store.items)
