from __future__ import annotations

import seaborn as sns

sns.set_theme()

from makefloat_opt.plt_gap import plot_gaps
from makefloat_opt.plt_status import plot_status
from makefloat_opt.plt_gantt import plot_gantt
from makefloat_opt.plt_runtime import plot_runtime
from makefloat_opt.schema import load_runs, SchedulerRun
from makefloat_opt.plt_schedule_diff import plot_diffs_from_runs

from loguru import logger


def report_validation(runs: list[SchedulerRun]):
    for r in runs:
        if not r.is_valid_baseline:
            logger.critical(f"BASELINE INVALID: {r.name}")
        if r.validation_milp is not None and not r.is_valid_milp:
            logger.critical(f"MILP INVALID: {r.name}")


if __name__ == '__main__':
    import sys
    logger.remove(0)
    logger.add(sys.stderr, level="CRITICAL")

    RUNS_FOLDER_MATCH = "../../workplace/RUNS/*"
    RUNS = load_runs(RUNS_FOLDER_MATCH)

    report_validation(RUNS)

    OptimalRUNS = [R for R in RUNS if R.gurobi_status == "Optimal"]

    logger.critical(">> OPTIMAL GAPS <<")
    plot_gaps(OptimalRUNS, figname="float/gaps_optimal.pdf", ms_hue=True).to_csv("gaps_table_optimal.csv", index=False)

    logger.critical(">> ALL GAPS <<")
    plot_gaps(RUNS, figname="float/gaps.pdf", ms_hue=True).to_csv("gaps_table.csv", index=False)

    logger.critical(">> DIFF <<")
    plot_diffs_from_runs(OptimalRUNS, figname="float/scheduler_diff.pdf").to_csv("schedule_diffs_optimal.csv", index=False)

    plot_status(RUNS, figname="float/status.pdf")

    plot_runtime(RUNS, figname="float/runtime.pdf")

    plot_gantt(runs_foler="../../workplace/RUNS", run_name="FDA-03-09-0-0", save_float_folder="./float",
               anno_reaction_index=True, multi_capacity=False)

    plot_gantt(runs_foler="../../workplace/RUNS", run_name="VS-04-06-1-1", save_float_folder="./float",
               anno_reaction_index=False, multi_capacity=True, figsize=(8, 4))
