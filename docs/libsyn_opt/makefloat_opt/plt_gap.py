from __future__ import annotations

import itertools
import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from loguru import logger

from makefloat_opt.schema import SchedulerRun, split_runs


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


def plot_gaps_ax(runs: list[SchedulerRun], subfig_title: str, prefix: str, ax: plt.Axes, remove_legend, ms_hue):
    df = pd.DataFrame.from_records([r.model_dump() for r in runs])
    df['percentage gap'] = get_gaps(runs)
    fms_translator = {
        "FMS-0": "LAB-1",
        "FMS-1": "LAB-2",
    }
    df['Module set'] = [fms_translator[f"FMS-{ifm}"] for ifm in df['ifm'].tolist()]

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
            logger.critical(
                f"GAPS {subfig_title} {ms} min:{gap_series.min()} max:{gap_series.max()}"
            )
            # df[['name', 'percentage gap']].to_csv(f"gaps--{subfig_title}.csv")
            # df[['name', 'gurobi_status', 'percentage gap']].to_csv(f"gaps--{subfig_title}.csv")
        # p = sns.catplot(df, x="percentage gap", y="n_target", ax=ax, hue="Module set", orient="h", alpha=.5, markers=markers, s=10, jitter=False)
    else:
        p = sns.stripplot(df, x="percentage gap", y="n_target", ax=ax, orient="h", alpha=.5, s=10, marker="D",
                          jitter=False)

    gaps = [v for v in df['percentage gap'].tolist() if not pd.isna(v)]
    logger.critical(f"mean: {sum(gaps) / len(gaps)}")
    ax.set_xlabel("Makespan gap (%)", fontsize=13)
    ax.set_ylabel("")
    labels = [item.get_text() for item in ax.get_yticklabels()]
    yticks = ax.get_yticks()
    if "FDA" == prefix:
        labels = [f"FDA.{int(l):02}" for l in labels]
    else:
        labels = [f"VS.{int(l):02}" for l in labels]
    ax.set_yticks(yticks, labels)
    if ms_hue:
        if remove_legend:
            try:
                ax.get_legend().remove()
            except AttributeError:
                pass
    # ax.set_title(subfig_title, fontsize='large', loc='center')
    ax.set_title(subfig_title, y=1.0, x=-0.085)
    return df[['name', 'gurobi_status', 'percentage gap']]


def plot_gaps(runs: list[SchedulerRun], figname: str, ms_hue: bool):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 10), sharey="row", sharex=True)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey="row", sharex=True)
    df_a = plot_gaps_ax(runs_fda_ws0, subfig_title="(A)", ax=ax1, remove_legend=True, ms_hue=ms_hue, prefix="FDA")
    df_b = plot_gaps_ax(runs_fda_ws1, subfig_title="(B)", ax=ax2, remove_legend=True, ms_hue=ms_hue, prefix="FDA")
    df_c = plot_gaps_ax(runs_vs_ws0, subfig_title="(C)", ax=ax3, remove_legend=True, ms_hue=ms_hue, prefix="VS")
    df_d = plot_gaps_ax(runs_vs_ws1, subfig_title="(D)", ax=ax4, remove_legend=True, ms_hue=ms_hue, prefix="VS")

    ax4.legend(loc='lower right', bbox_to_anchor=(1.0, -.45),
               fancybox=True, shadow=False, ncol=2)

    fig.supylabel("Chemical library", fontsize=13)
    fig.tight_layout()
    fig.savefig(figname, dpi=600)

    df = pd.concat([df_a, df_b, df_c, df_d], axis=0)
    fms_translator = {
        "0": "LAB-1",
        "1": "LAB-2",
    }
    records = []
    for r in df.to_dict("records"):
        name = r["name"]
        percentage_gap = r["percentage gap"]
        solver_status = r["gurobi_status"]
        lib, ntarget, nindex, fms, workshift = name.split("-")
        record = {
            "Library name": ".".join([lib, ntarget, nindex]),
            "Module set": fms_translator[fms],
            "Work shift": bool(int(workshift)),
            "Percentage gap": "N/A" if pd.isna(percentage_gap) else percentage_gap,
            "Optimal": True if solver_status == "Optimal" else False,
        }
        records.append(record)
    df = pd.DataFrame.from_records(records)
    return df
