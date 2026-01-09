# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_literal_selector_presence.py ###
from __future__ import annotations

import simpy
import pytest

from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.selector import LiteralSelector


def test_literal_selector_rejects_non_present_object():
    env = simpy.Environment()
    container = MaterialContainer(identifier="https://libsyn-sim/kg/not-present")
    container.is_present = {False}

    selector = LiteralSelector(container.identifier)
    with pytest.raises(ValueError, match="non-present"):
        env.process(selector.resolve(env))
        env.run()
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_literal_selector_presence.py ###
