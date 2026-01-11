from __future__ import annotations

import pytest

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from libsyn_tools.sim.knowledge_graph import (
    Has_capacity,
    Has_interrupt_events,
    Has_pool_type,
    MaterialContainer,
    canonical_iri,
)
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import (
    AddDataProperty,
    ChangeDataProperty,
    Create,
)


def _make_container(identifier: str, *, present: bool = True) -> MaterialContainer:
    container = MaterialContainer(identifier=identifier)
    container.is_present = {present}
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    return container


def test_canonical_iri_edits_apply_with_normalized_locks(env) -> None:
    eng = EffectEngine()
    container = _make_container("canon-container")
    edit = AddDataProperty(
        instance_1_iri=str(canonical_iri(container.identifier)),
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="ok",
    )

    eng.apply(
        [edit],
        env,
        operation_id="canon-edit",
        locked_iris=[str(canonical_iri(container.identifier))],
    )

    assert "ok" in container.has_interrupt_events


def test_rollback_truncates_recent_edits(env, monkeypatch) -> None:
    eng = EffectEngine()
    container = _make_container("rollback-container")
    eng._register_if_new(container, env)
    runtime_state = get_runtime_state(container, env)

    edit_ok = AddDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="first",
    )
    edit_fail = AddDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="second",
    )

    def boom() -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(edit_fail, "apply", boom)

    with pytest.raises(EngineMechanicalError):
        eng.apply_tx(
            [edit_ok, edit_fail],
            env,
            operation_id="rollback",
            locked_iris=[container.identifier],
        )

    assert len(runtime_state.recent_edits) == 0


def test_unhashable_data_value_rejected(env) -> None:
    eng = EffectEngine()
    container = _make_container("unhashable-container")
    edit = AddDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value=["nope"],
    )

    with pytest.raises(EngineMechanicalError, match="unhashable data_value"):
        eng.apply_tx(
            [edit],
            env,
            operation_id="unhashable",
            locked_iris=[container.identifier],
        )


def test_functional_data_property_add_is_rejected(env) -> None:
    eng = EffectEngine()
    container = _make_container("capacity-container")
    container.has_capacity = {1.0}
    edit = AddDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_capacity.predicate_iri,
        data_value=2.0,
    )

    with pytest.raises(EngineMechanicalError, match="functional data property"):
        eng.apply_tx(
            [edit],
            env,
            operation_id="capacity-add",
            locked_iris=[container.identifier],
        )


def test_create_requires_lock_for_runtime_objects(env) -> None:
    eng = EffectEngine()
    container = _make_container("create-lock-container", present=False)
    edit = Create(instance_1_iri=container.identifier)

    eng.apply_tx([edit], env, operation_id="create-no-lock", locked_iris=[])

    container_locked = _make_container("create-lock-container-locked", present=False)
    container_locked.has_pool_type.add("POOL_TX_CREATE")
    store = FilterStoreRegistry.get_filter_store("POOL_TX_CREATE", env)
    store.items.append(container_locked)
    edit_locked = Create(instance_1_iri=container_locked.identifier)

    with pytest.raises(EngineMechanicalError, match="requires lock"):
        eng.apply_tx([edit_locked], env, operation_id="create-needs-lock", locked_iris=[])

    eng.apply_tx(
        [edit_locked],
        env,
        operation_id="create-with-lock",
        locked_iris=[container_locked.identifier],
    )

    assert container_locked.is_present == {True}


def test_pool_type_change_removes_from_old_store(env) -> None:
    eng = EffectEngine()
    container = _make_container("pool-switch-container")
    container.has_pool_type.add("POOL_TX_A")
    eng._register_if_new(container, env)

    edit = ChangeDataProperty(
        instance_1_iri=container.identifier,
        property_iri=Has_pool_type.predicate_iri,
        data_value="POOL_TX_B",
    )

    eng.apply_tx(
        [edit],
        env,
        operation_id="pool-switch",
        locked_iris=[container.identifier],
    )

    store_a = FilterStoreRegistry.get_filter_store("POOL_TX_A", env)
    store_b = FilterStoreRegistry.get_filter_store("POOL_TX_B", env)
    assert container.identifier not in {obj.identifier for obj in store_a.items}
    assert container.identifier in {obj.identifier for obj in store_b.items}
