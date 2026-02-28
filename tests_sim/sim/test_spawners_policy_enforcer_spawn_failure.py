from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_shacl import SHACLViolationRecord
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit
from libsyn_tools.sim.spawner import PolicyEnforcerSpawner


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _violation(*, shape_iri: str, focus_iri: str, op_id: str = "source-op") -> SHACLViolationRecord:
    return SHACLViolationRecord(
        sim_time=0.0,
        operation_id=op_id,
        origin="SHACL",
        severity="soft",
        disposition="committed",
        shape_iri=shape_iri,
        focus_iri=focus_iri,
    )


def test_policy_enforcer_spawn_failure_consumes_failed_quota_only() -> None:
    shape = "urn:shape:spawn-failure"
    focus = MaterialContainer(identifier="https://libsyn-sim/kg/focus-spawn-failure")
    KnowledgeGraph.get_object_from_lookup(focus.identifier)
    Create(instance_1_iri=focus.identifier).apply()

    calls = {"n": 0}

    def _factory(_record: SHACLViolationRecord):
        calls["n"] += 1
        if calls["n"] == 1:
            return _NoOp(identifier="existing-op")
        return _NoOp(identifier="spawn-ok-1")

    sim = Simulation([_NoOp(identifier="existing-op")])
    spawner = PolicyEnforcerSpawner(
        shape_dispatch={shape: _factory},
        dedupe=False,
        max_spawned_remediations_per_focus=1,
        max_failed_remediations_per_focus=2,
    )
    spawner.attach(sim)

    sim.callbacks.emit_violation(_violation(shape_iri=shape, focus_iri=focus.identifier))
    sim.callbacks.emit_violation(_violation(shape_iri=shape, focus_iri=focus.identifier))
    sim.callbacks.emit_violation(_violation(shape_iri=shape, focus_iri=focus.identifier))

    assert calls["n"] == 2
    assert "spawn-ok-1" in sim.operation_registry
    assert "spawn-ok-2" not in sim.operation_registry


def test_policy_enforcer_failed_quota_caps_spawn_failures() -> None:
    shape = "urn:shape:spawn-failure-cap"
    focus = MaterialContainer(identifier="https://libsyn-sim/kg/focus-failure-cap")
    KnowledgeGraph.get_object_from_lookup(focus.identifier)
    Create(instance_1_iri=focus.identifier).apply()

    calls = {"n": 0}

    def _factory(_record: SHACLViolationRecord):
        calls["n"] += 1
        return _NoOp(identifier="existing-op")

    sim = Simulation([_NoOp(identifier="existing-op")])
    spawner = PolicyEnforcerSpawner(
        shape_dispatch={shape: _factory},
        dedupe=False,
        max_spawned_remediations_per_focus=5,
        max_failed_remediations_per_focus=1,
    )
    spawner.attach(sim)

    for _ in range(3):
        sim.callbacks.emit_violation(_violation(shape_iri=shape, focus_iri=focus.identifier))

    assert calls["n"] == 1
    assert "existing-op" in sim.operation_registry
    assert len([op_id for op_id in sim.operation_registry if op_id.startswith("spawn-ok-")]) == 0
