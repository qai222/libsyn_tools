# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_simulation_speed_factor_validation.py ###
from __future__ import annotations

import pytest

from libsyn_tools.sim import Simulation


@pytest.mark.parametrize("value", [0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_simulation_speed_factor_invalid_rejected(value: float) -> None:
    with pytest.raises(ValueError, match="simulation_speed_factor"):
        Simulation([], simulation_speed_factor=value)
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_simulation_speed_factor_validation.py ###
