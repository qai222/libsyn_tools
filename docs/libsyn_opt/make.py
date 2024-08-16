from __future__ import annotations

from makefloat_opt.plt_gap import plot_gaps
from makefloat_opt.plt_status import plot_status
from makefloat_opt.plt_gantt import plot_gantt
from makefloat_opt.schema import load_runs

if __name__ == '__main__':
    RUNS_FOLDER_MATCH = "../../workplace/RUNS/*"
    RUNS = load_runs(RUNS_FOLDER_MATCH)
    OptimalRUNS = [R for R in RUNS if R.gurobi_status == "Optimal"]

    # # plot_gaps(RUNS, figname="float/gaps.pdf", ms_hue=True)
    # # print("OPTIMAL RUNS")
    # plot_gaps(OptimalRUNS, figname="float/gaps_optimal.pdf", ms_hue=True)

    # plot_status(RUNS, figname="float/status.pdf")

    plot_gantt(runs_foler="../../workplace/RUNS", run_name="FDA-03-09-0-0", save_float_folder="./float", anno_reaction_index=True)

    # plot_gantt(runs_foler="../../workplace/RUNS", run_name="VS-04-06-0-1", save_float_folder="./float", anno_reaction_index=False)
