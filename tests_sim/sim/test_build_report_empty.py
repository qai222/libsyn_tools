# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_build_report_empty.py ###
from __future__ import annotations

from libsyn_tools.sim import Simulation


def test_build_report_empty_history():
    sim = Simulation([])
    report = sim.build_report()

    assert report.summary["operation_counts"]["start"] == 0
    assert report.summary["operation_counts"]["end"] == 0
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_build_report_empty.py ###
