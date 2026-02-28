from __future__ import annotations

from pathlib import Path

from libsyn_tools.chem_schema import Operation, OperationType
from libsyn_tools.opt import SchedulerOutput
from libsyn_tools.sim.adapters.schedule_bridge import compile_schedule_to_simulation


def run(out_dir: str | Path | None = None):
    planned = [Operation(identifier="op1", type=OperationType.TransferLiquid)]
    schedule = SchedulerOutput(
        start_times={"op1": 0.0},
        end_times={"op1": 1.0},
        assignments={"op1": "module_1"},
    )
    sim = compile_schedule_to_simulation(planned, schedule)
    if out_dir:
        sim.run_and_report(out_dir)
    else:
        sim.run()
    return sim


if __name__ == "__main__":
    run()
