from __future__ import annotations

from libsyn_tools.sim.simulation import Simulation
from libsyn_tools.sim.overlay.chemistry_overlay import ChemistryOverlayProvider


def install_chemistry_pack(sim: Simulation) -> None:
    provider = ChemistryOverlayProvider()
    sim.effect_engine.register_overlay_provider(provider.snapshot)
