from __future__ import annotations

import pytest
from rdflib import Graph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_engine import EngineMechanicalError


def test_strict_overlays_fail_validation_by_default_with_shapes():
    shapes_graph = Graph()
    sim = Simulation([], shacl_shapes=shapes_graph)

    def _broken_overlay():
        raise RuntimeError("overlay broke")

    sim.effect_engine.register_overlay_provider(_broken_overlay)

    with pytest.raises(EngineMechanicalError):
        sim.effect_engine.validate_now()


def test_tolerant_overlays_continue_when_disabled():
    shapes_graph = Graph()
    sim = Simulation([], shacl_shapes=shapes_graph, strict_overlays=False)

    def _broken_overlay():
        raise RuntimeError("overlay broke")

    sim.effect_engine.register_overlay_provider(_broken_overlay)

    conforms, _, _ = sim.effect_engine.validate_now()
    assert conforms is True
