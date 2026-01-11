from __future__ import annotations

import pytest

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer, PortionOfMaterial
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.remediation import make_drain_to_capacity


class _NoOp(Operation):
    participant_container: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _pom_with_volume(volume: float) -> PortionOfMaterial:
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=volume, density=1.0, molecular_weight=18.0))
    return pom


def test_pom_zero_volume_portioning_raises() -> None:
    pom = _pom_with_volume(0.0)
    with pytest.raises(ValueError):
        pom.get_portion(0.5)
    with pytest.raises(ValueError):
        pom.get_portion_by_volume(0.1)


def test_pom_invalid_portion_size_raises() -> None:
    pom = _pom_with_volume(1.0)
    with pytest.raises(ValueError):
        pom.get_portion(0.0)
    with pytest.raises(ValueError):
        pom.get_portion(1.5)


def test_make_drain_to_capacity_respects_capacity() -> None:
    focus = MaterialContainer(identifier="focus-cap")
    waste = MaterialContainer(identifier="waste-cap")
    for obj in (focus, waste):
        obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)

    focus.has_capacity.add(1.0)
    op = make_drain_to_capacity(focus.identifier, waste.identifier)
    assert op.target_volume <= focus.capacity

    focus_small = MaterialContainer(identifier="focus-small")
    focus_small.is_present = {True}
    focus_small.has_capacity.add(1e-9)
    KnowledgeGraph.get_object_from_lookup(focus_small.identifier)
    op_small = make_drain_to_capacity(focus_small.identifier, waste.identifier)
    assert op_small.target_volume <= focus_small.capacity


def test_container_capacity_missing_raises() -> None:
    container = MaterialContainer(identifier="no-capacity")
    with pytest.raises(ValueError, match="has_capacity for no-capacity is invalid"):
        _ = container.capacity


def test_container_capacity_multiple_values_raises() -> None:
    container = MaterialContainer(identifier="cap-multi")
    container.has_capacity.update({1.0, 2.0})
    with pytest.raises(ValueError, match="has_capacity for cap-multi is invalid"):
        _ = container.capacity


def test_negative_temporal_cost_raises_early() -> None:
    container = MaterialContainer(identifier="neg-temp")
    container.has_pool_type.add("POOL_NEG_TEMP")
    container.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(container.identifier)

    op = _NoOp(participant_container=container.identifier)
    op.temporal_cost = -1.0

    sim = Simulation([op])
    sim.run()

    rs = get_runtime_state(container, sim.env)
    store = FilterStoreRegistry.get_filter_store("POOL_NEG_TEMP", sim.env)
    assert rs.lock.count == 0
    assert container in store.items
    assert any(event.event_type == "OPERATION_ABORT" for event in sim.history_log)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_temporal_cost_nonfinite_assignment_aborts(value: float) -> None:
    pool_type = f"POOL_TEMP_INVALID_{value}"
    container = MaterialContainer(identifier=f"temp-invalid-{value}")
    container.has_pool_type.add(pool_type)
    container.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(container.identifier)

    op = _NoOp(participant_container=container.identifier)
    op.temporal_cost = value

    sim = Simulation([op])
    sim.run()

    rs = get_runtime_state(container, sim.env)
    store = FilterStoreRegistry.get_filter_store(pool_type, sim.env)
    assert rs.lock.count == 0
    assert container in store.items
    assert any(event.event_type == "OPERATION_ABORT" for event in sim.history_log)


@pytest.mark.parametrize("value", [-0.5, float("nan"), float("inf"), float("-inf")])
def test_temporal_cost_invalid_init_raises(value: float) -> None:
    with pytest.raises(ValueError):
        _NoOp(participant_container="x", temporal_cost=value)
