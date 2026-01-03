from __future__ import annotations

from typing import Iterable

import pytest

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject, MaterialContainer, PortionOfMaterial
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
