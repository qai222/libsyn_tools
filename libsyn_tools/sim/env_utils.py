from __future__ import annotations

from typing import TYPE_CHECKING, cast

import simpy

if TYPE_CHECKING:
    from libsyn_tools.sim.effect_engine import EffectEngine


def get_effect_engine(env: simpy.Environment) -> "EffectEngine":
    engine = getattr(env, "_libsyn_effect_engine", None)
    if engine is None:
        raise RuntimeError(
            "EffectEngine is not attached to this environment. "
            "Construct the environment via Simulation so it is registered."
        )
    return cast("EffectEngine", engine)
