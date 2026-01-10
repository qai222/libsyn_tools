from __future__ import annotations

import pytest

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph import PortionOfMaterial


def test_get_portion_by_volume_clamps_within_tolerance() -> None:
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=2.0, density=1.0))

    portion = pom.get_portion_by_volume(pom.volume + 1e-10)
    assert portion.volume == pytest.approx(pom.volume)


def test_get_portion_by_volume_rejects_out_of_tolerance() -> None:
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=2.0, density=1.0))

    with pytest.raises(ValueError):
        pom.get_portion_by_volume(pom.volume + 1e-3)
