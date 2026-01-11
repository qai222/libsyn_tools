from __future__ import annotations

import json
import math

import pytest

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.knowledge_graph import MaterialContainer, PortionOfMaterial, canonical_iri
from libsyn_tools.sim.overlay.chemistry_overlay import ChemistryOverlayProvider, LIB


def test_chemical_density_zero_raises() -> None:
    chem = Chemical(mass=1.0, density=0.0)
    with pytest.raises(ValueError, match="density"):
        _ = chem.volume


@pytest.mark.parametrize(
    ("mass", "density", "pattern"),
    [
        (None, 1.0, "mass"),
        (1.0, None, "density"),
        (1.0, math.inf, "density"),
        (math.inf, 1.0, "mass"),
    ],
)
def test_chemical_mass_density_invalid_raises(mass, density, pattern: str) -> None:
    chem = Chemical(mass=mass, density=density)
    with pytest.raises(ValueError, match=pattern):
        _ = chem.volume


def test_material_container_capacity_requires_finite_positive() -> None:
    container = MaterialContainer(identifier="cap-invalid")
    container.has_capacity.add(0.0)
    with pytest.raises(ValueError, match="finite and > 0"):
        _ = container.capacity

    container.has_capacity = {math.inf}
    with pytest.raises(ValueError, match="finite and > 0"):
        _ = container.capacity


def test_split_and_merge_conserves_volume_within_tolerance() -> None:
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=10.0, density=2.0))
    pom.add_chemical(Chemical(mass=5.0, density=1.0))
    original_volume = pom.volume

    for _ in range(5):
        part_a = pom.get_portion(0.33)
        part_b = pom.get_portion(0.67)
        pom = part_a.mix_with(part_b)

    assert pom.volume == pytest.approx(original_volume, abs=1e-6)


def test_overlay_skips_malformed_ingredient_blob() -> None:
    pom = PortionOfMaterial(identifier="pom-bad-ingredient")
    pom.is_present = {True}
    valid_blob = json.dumps(
        PortionOfMaterial._normalize_chemical_payload(Chemical(mass=1.0, density=1.0)),
        sort_keys=True,
    )
    pom.has_ingredient = {"not-json", valid_blob}

    graph = ChemistryOverlayProvider().snapshot()
    pom_iri = canonical_iri(pom.identifier)
    assert (pom_iri, LIB.hasIngredient, None) in graph
