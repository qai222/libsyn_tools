from __future__ import annotations

import glob
import itertools
import math
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import pandas as pd
import seaborn as sns
from pydantic import BaseModel
from tqdm import tqdm

from libsyn_tools.opt.formulation_baseline import SolverBaseline
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.utils import json_load, FilePath, json_dump

sns.set_theme()

# https://www.gurobi.com/documentation/current/refman/optimization_status_codes.html
gurobi_status_dict = {
    9: "Feasible",
    3: "Infeasible",
    2: "Optimal",
}


class SchedulerRun(BaseModel):
    folder: str
    ifm: int
    name: str
    prefix: str
    network_index: int
    n_target: int
    workshift: bool

    gurobi_status: str | None

    validation_milp: dict | None
    validation_baseline: dict

    solution_baseline: float
    solution_milp: float | None

    runtime_milp: float | None

    @property
    def percentage_gap(self) -> float | None:
        if self.solution_milp is None:
            return None
        else:
            return round(100 * (self.solution_baseline - self.solution_milp) / self.solution_baseline, 1)

    @classmethod
    def from_run_folder(cls, run_folder: FilePath):

        assert os.path.isdir(run_folder)

        name = os.path.basename(run_folder)
        assert name
        prefix, x, y, ifm, shift = name.split('-')

        solver_baseline = os.path.join(run_folder, "solver_baseline.json")
        solver_baseline = SolverBaseline(**json_load(solver_baseline))
        # revalidate just to be safe
        solver_baseline.output.notes['validation'] = solver_baseline.output.validate_schedule(solver_baseline.input)

        solver_milp = os.path.join(run_folder, "solver_milp.json")
        if os.path.exists(solver_milp) and json_load(solver_milp)['output']['end_times']:
            solver_milp = SolverMILP(**json_load(solver_milp))
            solver_milp.output.notes['validation'] = solver_milp.output.validate_schedule(solver_milp.input)
            # revalidate just to be safe
            validation_milp = solver_baseline.output.notes['validation']
            solution_milp = max(solver_milp.output.end_times.values())
            gurobi_status = gurobi_status_dict[solver_milp.opt_log['gurobi status']]
            runtime_milp = solver_milp.opt_log['time solved']
        else:
            validation_milp = None
            solution_milp = None
            gurobi_status = None
            runtime_milp = None

        return cls(
            folder=run_folder,
            ifm=int(ifm),
            name=name,
            prefix=prefix,
            network_index=int(y),
            n_target=int(x),
            workshift=bool(int(shift)),
            validation_milp=validation_milp,
            validation_baseline=solver_baseline.output.notes['validation'],
            solution_baseline=max(solver_baseline.output.end_times.values()),
            solution_milp=solution_milp,
            gurobi_status=gurobi_status,
            runtime_milp=runtime_milp,
        )


def get_gaps(
        runs: list[SchedulerRun],
        min_cap=-math.inf,
        min_fill=None,
        # min_cap=-20,
        # min_fill = -20,
):
    gaps = []
    for r in runs:
        pgap = r.percentage_gap
        if pgap is None:
            gaps.append(pgap)
        elif pgap >= min_cap:
            gaps.append(pgap)
        else:
            gaps.append(min_fill)
    return gaps


def plot_gaps_ax(runs: list[SchedulerRun], subfig_title: str, ax: plt.Axes, remove_legend, ms_hue):
    df = pd.DataFrame.from_records([r.model_dump() for r in runs])
    df['percentage gap'] = get_gaps(runs)
    df['Module set'] = [f"FMS-{ifm}" for ifm in df['ifm'].tolist()]

    marker_shape = itertools.cycle(['D', 'o', 'p', 'd', 'v', 'D', '*', 'o', ])
    color_selector = itertools.cycle(
        ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
         'tab:olive', 'tab:cyan', ])
    markers = {i: next(marker_shape) for i in df["Module set"].unique()}
    colors = {i: next(color_selector) for i in df["Module set"].unique()}

    if ms_hue:
        for ms in df["Module set"].unique():
            sns.stripplot(
                df[df["Module set"] == ms], x="percentage gap", y="n_target", ax=ax,
                dodge=True,
                hue="Module set", orient="h", alpha=.5, marker=markers[ms], s=10, jitter=False, linewidth=1,
                palette=[colors[ms]]
            )
            gap_series = df[df["Module set"] == ms]['percentage gap']
            print(subfig_title, ms, gap_series.min(), gap_series.max(), gap_series.mean())
            df[['name', 'percentage gap']].to_csv(f"gaps--{subfig_title}.csv")
        # p = sns.catplot(df, x="percentage gap", y="n_target", ax=ax, hue="Module set", orient="h", alpha=.5, markers=markers, s=10, jitter=False)
    else:
        p = sns.stripplot(df, x="percentage gap", y="n_target", ax=ax, orient="h", alpha=.5, s=10, marker="D",
                          jitter=False)

    gaps = [v for v in df['percentage gap'].tolist() if not pd.isna(v)]
    print(subfig_title, sum(gaps) / len(gaps))
    ax.set_xlabel("Gap (%)")
    ax.set_ylabel("Library size")
    if ms_hue:
        if remove_legend:
            try:
                ax.get_legend().remove()
            except AttributeError:
                pass
    # ax.set_title(subfig_title, fontsize='large', loc='center')
    ax.set_title(subfig_title, loc='left')


def split_runs(runs: list[SchedulerRun]):
    runs_fda_ws0 = []
    runs_fda_ws1 = []
    runs_vs_ws0 = []
    runs_vs_ws1 = []

    for r in runs:
        assert r.validation_baseline['lmin and lmax precedence valid'] is False
        if r.workshift:
            if r.prefix == "FDA":
                runs_fda_ws1.append(r)
            else:
                runs_vs_ws1.append(r)
        else:
            if r.prefix == "FDA":
                runs_fda_ws0.append(r)
            else:
                runs_vs_ws0.append(r)
    return runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1


def plot_gaps(runs: list[SchedulerRun], figname: str, ms_hue: bool):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 10), sharey="row", sharex=True)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey="row", sharex=True)
    plot_gaps_ax(runs_fda_ws0, subfig_title="(A) FDA", ax=ax1, remove_legend=False, ms_hue=ms_hue)
    plot_gaps_ax(runs_fda_ws1, subfig_title="(B) FDA + work shifts", ax=ax2, remove_legend=True, ms_hue=ms_hue)
    plot_gaps_ax(runs_vs_ws0, subfig_title="(C) VS", ax=ax3, remove_legend=True, ms_hue=ms_hue)
    plot_gaps_ax(runs_vs_ws1, subfig_title="(D) VS + work shifts", ax=ax4, remove_legend=True, ms_hue=ms_hue)

    fig.tight_layout()
    fig.savefig(figname, dpi=600)


def plot_status_ax(runs: list[SchedulerRun], subfig_title: str, ax: plt.Axes, remove_legend):
    df = pd.DataFrame.from_records([r.model_dump() for r in runs])
    df['Solver status'] = df['gurobi_status'].fillna(value='No solution')
    p = sns.histplot(df, y="n_target", ax=ax, hue="Solver status", multiple="stack", discrete=True, )
    ax.set_ylabel("Library size")
    ax.set_xlabel("Scheduler instance")
    if remove_legend:
        ax.get_legend().remove()
    ax.xaxis.set_major_locator(tck.MultipleLocator())
    ax.yaxis.set_major_locator(tck.MultipleLocator(base=1))
    ax.set_xlim([0, 20])
    # ax.set_title(subfig_title, fontsize='large', loc='center')
    ax.set_title(subfig_title, loc='left')


def plot_status(runs: list[SchedulerRun], figname: str):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 10), sharey=True, sharex=True)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey="row", sharex=True)
    plot_status_ax(runs_fda_ws0, subfig_title="(A) FDA", ax=ax1, remove_legend=False)
    plot_status_ax(runs_fda_ws1, subfig_title="(B) FDA + work shifts", ax=ax2, remove_legend=True)
    plot_status_ax(runs_vs_ws0, subfig_title="(C) VS", ax=ax3, remove_legend=True)
    plot_status_ax(runs_vs_ws1, subfig_title="(D) VS + work shifts", ax=ax4, remove_legend=True)

    fig.tight_layout()
    fig.savefig(figname, dpi=600)


def plot_runtime_ax(runs: list[SchedulerRun], subfig_title: str, ax: plt.Axes, remove_legend):
    df = pd.DataFrame.from_records([r.model_dump() for r in runs])
    p = sns.stripplot(
        df, x="n_target", y="runtime_milp", ax=ax, log_scale=True,
        # dodge = True,
        # orient = "h",
        alpha=.5, s=10, jitter=False, linewidth=1, marker="x"

    )
    ax.axhline(y=3600, label="1 hr", color="red", ls=":")
    # ax.axhline(y=10800, label="3 hr", color="red")
    ax.set_ylabel("MILP runtime")
    ax.set_xlabel("Library size")
    ax.legend()
    if remove_legend:
        try:
            ax.get_legend().remove()
        except AttributeError:
            pass
    # ax.xaxis.set_major_locator(tck.MultipleLocator())
    # ax.yaxis.set_major_locator(tck.MultipleLocator(base=1))
    # ax.set_xlim([0, 20])
    # ax.set_title(subfig_title, fontsize='large', loc='center')
    ax.set_title(subfig_title, loc='left')


def plot_runtime(runs: list[SchedulerRun], figname: str):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # runs_fda = runs_fda_ws0 + runs_fda_ws1
    # runs_vs = runs_vs_ws0 + runs_vs_ws1
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey=True, sharex=True)
    plot_runtime_ax(runs_fda_ws0, subfig_title="(A) FDA", ax=ax1, remove_legend=False)
    plot_runtime_ax(runs_fda_ws1, subfig_title="(B) FDA + work shifts", ax=ax2, remove_legend=True)
    plot_runtime_ax(runs_vs_ws0, subfig_title="(C) VS", ax=ax3, remove_legend=True)
    plot_runtime_ax(runs_vs_ws1, subfig_title="(D) VS + work shifts", ax=ax4, remove_legend=True)
    fig.tight_layout()
    fig.savefig(figname, dpi=600)


if __name__ == '__main__':
    RUNS_FOLDER_MATCH = "../workplace/RUNS/*"
    try:
        RUNS = [SchedulerRun(**R) for R in json_load("RUNS.json")]
    except FileNotFoundError:
        RUNS = []
        for FOLDER in tqdm(sorted(glob.glob(RUNS_FOLDER_MATCH))):
            if not os.path.isdir(FOLDER):
                continue
            if not os.path.exists(os.path.join(FOLDER, "solver_baseline.json")):
                continue
            RUN = SchedulerRun.from_run_folder(FOLDER)
            if RUN.ifm > 1:
                continue
            RUNS.append(RUN)
        json_dump([R.model_dump() for R in RUNS], "RUNS.json")

    plot_status(RUNS, figname="float/status.pdf")
    plot_gaps(RUNS, figname="float/gaps.pdf", ms_hue=True)
    plot_runtime(RUNS, figname="float/runtimes.pdf")
    print("OPTIMAL RUNS")
    OptimalRUNS = [R for R in RUNS if R.gurobi_status == "Optimal"]
    plot_gaps(OptimalRUNS, figname="float/gaps_optimal.pdf", ms_hue=True)
    plot_runtime(OptimalRUNS, figname="float/runtimes_optimal.pdf")
