from __future__ import annotations

import pytest

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph import PortionOfMaterial


def _pom_with_volume(volume: float) -> PortionOfMaterial:
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=volume, density=1.0))
    return pom


def test_get_portion_by_volume_clamps_to_total_volume() -> None:
    pom = _pom_with_volume(10.0)
    portion = pom.get_portion_by_volume(10.0 + 5e-10)
    assert portion.volume == pytest.approx(pom.volume)


def test_pom_math_rejects_missing_mass_or_density() -> None:
    pom = PortionOfMaterial()
    with pytest.raises(ValueError, match="mass"):
        pom.add_chemical(Chemical(mass=None, density=1.0))

    with pytest.raises(ValueError, match="density"):
        pom.add_chemical(Chemical(mass=1.0, density=None))


def test_split_merge_cycles_conserve_volume() -> None:
    original = _pom_with_volume(12.0)
    pom = original
    for _ in range(5):
        part = pom.get_portion(0.37)
        rest = pom.get_portion(0.63)
        pom = part.mix_with(rest)
    assert pom.volume == pytest.approx(original.volume, rel=1e-9, abs=1e-9)
