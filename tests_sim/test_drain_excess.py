from __future__ import annotations

from typing import Iterable

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject, MaterialContainer, PortionOfMaterial
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create
from libsyn_tools.sim.operation_preset.drain import DrainExcess


def _make_container(*, capacity: float = 100.0, present: bool = True) -> MaterialContainer:
    c = MaterialContainer()
    c.has_capacity.add(capacity)
    c.is_present = {present}
    return c


def _make_pom(*, volume: float, container: MaterialContainer) -> PortionOfMaterial:
    # Use density=1 g/mL so mass == volume in mL.
    chem = Chemical(
        smiles="O",
        molecular_weight=18.015,
        density=1.0,
        mass=volume,
    )
    pom = PortionOfMaterial()
    pom.add_chemical(chem)
    pom.is_present = {True}
    pom.is_directly_contained_by.add(container)
    return pom


def _direct_poms(container: MaterialContainer) -> list[PortionOfMaterial]:
    return LabObject.get_directly_contained_individuals(container, PortionOfMaterial, only_present=True)


def test_drain_excess_handles_full_pom_drain_and_optional_device() -> None:
    """Regression: draining a full POM previously attempted to create a 0-volume residual POM."""

    src = _make_container(capacity=100.0)
    dst = _make_container(capacity=100.0)

    # Two POMs so that the first one is fully drained (residual volume == 0).
    _make_pom(volume=5.0, container=src)
    _make_pom(volume=5.0, container=src)

    op = DrainExcess(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        # participant_device intentionally omitted (None)
        target_volume=3.0,
    )

    sim = Simulation([op])
    sim.run()

    assert pytest.approx(src.directly_contained_pom_volume, rel=1e-6, abs=1e-6) == 3.0
    assert pytest.approx(dst.directly_contained_pom_volume, rel=1e-6, abs=1e-6) == 7.0

    # Ensure we didn't create any present 0-volume POMs as a "residual".
    for container in (src, dst):
        for pom in _direct_poms(container):
            assert pom.volume > 1e-6


def test_drain_excess_rejects_non_container_destination() -> None:
    src = _make_container(capacity=100.0)
    dst = LabObject()
    dst.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(dst.identifier)
    Create(instance_1_iri=dst.identifier).apply()

    op = DrainExcess(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        target_volume=1.0,
    )

    sim = Simulation([op])
    sim.run()
    aborts = [record for record in sim.history_log if record.event_type == "OPERATION_ABORT"]
    assert aborts
    assert "destination must be a MaterialContainer" in aborts[0].operation_data.get("error", "")

    for obj in (src, dst):
        rs = get_runtime_state(obj, sim.env)
        assert rs.lock.count == 0
        assert not rs.lock.users
        assert not rs.lock.queue
