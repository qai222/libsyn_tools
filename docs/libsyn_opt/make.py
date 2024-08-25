from __future__ import annotations

import seaborn as sns
sns.set_theme()

from makefloat_opt.plt_gap import plot_gaps
from makefloat_opt.plt_status import plot_status
from makefloat_opt.plt_gantt import plot_gantt
from makefloat_opt.plt_runtime import plot_runtime
from makefloat_opt.plt_assign import plot_assign_diff_all
from makefloat_opt.schema import load_runs, SchedulerRun


def report_validation(runs: list[SchedulerRun]):

    try:
        assert all([r.validation_baseline["lmin and lmax precedence valid"] is False for r in runs])
    except AssertionError:
        print("not all baseline violate lmax!")

    try:
        assert all([r.is_valid_milp for r in runs if r.validation_milp is not None])
    except AssertionError:
        print("not all milp valid!")



if __name__ == '__main__':
    RUNS_FOLDER_MATCH = "../../workplace/RUNS/*"
    RUNS = load_runs(RUNS_FOLDER_MATCH)
    report_validation(RUNS)

    OptimalRUNS = [R for R in RUNS if R.gurobi_status == "Optimal"]

    # plot_assign_diff_all(OptimalRUNS, save_float_folder="float")

    # plot_gaps(RUNS, figname="float/gaps.pdf", ms_hue=True)
    # plot_gaps(OptimalRUNS, figname="float/gaps_optimal.pdf", ms_hue=True)

    plot_status(RUNS, figname="float/status.pdf")
    #
    # plot_runtime(RUNS, figname="float/runtime.pdf")

    # plot_gantt(runs_foler="../../workplace/RUNS", run_name="FDA-03-09-0-0", save_float_folder="./float", anno_reaction_index=True, multi_capacity=False)
    #
    # plot_gantt(runs_foler="../../workplace/RUNS", run_name="VS-04-06-1-1", save_float_folder="./float", anno_reaction_index=False, multi_capacity=True)
