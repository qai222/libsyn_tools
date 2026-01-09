# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_simulation_speed_factor_validation.py ###
from __future__ import annotations

import pytest

from libsyn_tools.sim import Simulation


def test_simulation_speed_factor_zero_rejected():
    with pytest.raises(ValueError, match="simulation_speed_factor"):
        Simulation([], simulation_speed_factor=0)


def test_simulation_speed_factor_negative_rejected():
    with pytest.raises(ValueError, match="simulation_speed_factor"):
        Simulation([], simulation_speed_factor=-1.0)
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_simulation_speed_factor_validation.py ###
