from __future__ import annotations

import pytest
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit
from libsyn_tools.sim.spawner import KGInspectorSpawner, TimerSpawner


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class _SelectorRaisesOp(Operation):
    participant_target: AttributeSelector = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_periodic_spawner_requires_until() -> None:
    sim = Simulation([])

    TimerSpawner(op_factory=lambda _sim: _NoOp(), interval=1.0).attach(sim)
    with pytest.raises(ValueError, match="requires an explicit until"):
        sim.run()


def test_periodic_inspector_requires_until() -> None:
    sim = Simulation([])

    KGInspectorSpawner(shape_dispatch={}, inspect_interval=1.0).attach(sim)
    with pytest.raises(ValueError, match="requires an explicit until"):
        sim.run()


def test_selector_predicate_exception_does_not_drain_pool() -> None:
    obj = LabObject()
    obj.is_present = {True}
    obj.has_pool_type.add("POOL_SELECTOR_EXCEPTION")
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()

    predicate_state = {"calls": 0}

    def _predicate(_candidate: LabObject) -> bool:
        predicate_state["calls"] += 1
        if predicate_state["calls"] > 1:
            raise RuntimeError("boom")
        return True

    selector = AttributeSelector(
        pool_type="POOL_SELECTOR_EXCEPTION",
        predicate=_predicate,
    )
    op = _SelectorRaisesOp(participant_target=selector)
    sim = Simulation([op])
    sim.run()

    rs = get_runtime_state(obj, sim.env)
    store = FilterStoreRegistry.get_filter_store("POOL_SELECTOR_EXCEPTION", sim.env)
    assert rs.lock.count == 0
    assert not rs.lock.queue
    assert obj in store.items
    assert any(event.event_type == "OPERATION_ABORT" for event in sim.history_log)


def test_spawn_operation_does_not_double_start() -> None:
    counter: dict[str, int] = {"count": 0}

    class _CountOp(_NoOp):
        def post_act(self, env) -> None:
            super().post_act(env)
            counter["count"] += 1

    sim = Simulation([])
    op = _CountOp(identifier="COUNT_ONCE")
    sim.spawn_operation(op, start_immediately=True)
    sim.run(until=1.0)

    starts = [
        event for event in sim.history_log
        if event.event_type == "OPERATION_START" and event.operation_id == "COUNT_ONCE"
    ]
    assert len(starts) == 1
    assert counter["count"] == 1
