# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_literal_selector_unknown_iri.py ###
from __future__ import annotations

import simpy
import pytest

from libsyn_tools.sim.operation.selector import LiteralSelector


def test_literal_selector_unknown_iri_raises_value_error():
    env = simpy.Environment()
    selector = LiteralSelector("https://libsyn-sim/kg/does-not-exist")
    with pytest.raises(ValueError, match="does-not-exist"):
        env.process(selector.resolve(env))
        env.run()
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_literal_selector_unknown_iri.py ###
