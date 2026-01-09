from __future__ import annotations

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class RuntimeErrorOp(Operation):
    participant_resource: StrOrSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        raise RuntimeError("boom")


class NoopOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_runtime_error_aborts_and_cleans_up() -> None:
    container = MaterialContainer()
    container.has_pool_type.add("POOL_RUNTIME_ERROR")
    container.is_present = {True}

    op = RuntimeErrorOp(participant_resource=container.identifier)
    sim = Simulation([op])
    sim.run()

    rs = get_runtime_state(container, sim.env)
    assert rs.lock.count == 0

    store = FilterStoreRegistry.get_filter_store("POOL_RUNTIME_ERROR", sim.env)
    assert container in store.items

    abort_events = [event for event in sim.history_log if event.event_type == "OPERATION_ABORT"]
    assert abort_events
    assert abort_events[0].operation_data["error_type"] == "RuntimeError"

    report = sim.build_report()
    assert report.summary["operation_counts"]["abort"] == 1
    assert report.summary["operation_counts"]["in_progress"] == 0


def test_callback_exception_is_logged() -> None:
    op = NoopOp()
    sim = Simulation([op])

    def bad_callback(proc) -> None:
        raise RuntimeError("callback boom")

    sim.callbacks.on_operation_end.append(bad_callback)
    sim.run()

    assert any(event.event_type == "CALLBACK_ERROR" for event in sim.history_log)
    assert any(event.event_type == "OPERATION_END" for event in sim.history_log)


def test_build_report_empty_history() -> None:
    sim = Simulation([])
    report = sim.build_report()

    assert report.event_log.empty
    assert report.summary["operation_counts"]["start"] == 0
    assert report.summary["operation_counts"]["in_progress"] == 0
