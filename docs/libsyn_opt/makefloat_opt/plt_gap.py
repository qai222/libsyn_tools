from __future__ import annotations

import itertools
import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from makefloat_opt.schema import SchedulerRun, split_runs

sns.set_theme()


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
            # df[['name', 'percentage gap']].to_csv(f"gaps--{subfig_title}.csv")
        # p = sns.catplot(df, x="percentage gap", y="n_target", ax=ax, hue="Module set", orient="h", alpha=.5, markers=markers, s=10, jitter=False)
    else:
        p = sns.stripplot(df, x="percentage gap", y="n_target", ax=ax, orient="h", alpha=.5, s=10, marker="D",
                          jitter=False)

    gaps = [v for v in df['percentage gap'].tolist() if not pd.isna(v)]
    print(subfig_title, sum(gaps) / len(gaps))
    ax.set_xlabel("Gap (%)")
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
    ax.set_title(subfig_title, loc='left')


def plot_gaps(runs: list[SchedulerRun], figname: str, ms_hue: bool):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 10), sharey="row", sharex=True)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey="row", sharex=True)
    plot_gaps_ax(runs_fda_ws0, subfig_title="(A)", ax=ax1, remove_legend=True, ms_hue=ms_hue, prefix="FDA")
    plot_gaps_ax(runs_fda_ws1, subfig_title="(B)", ax=ax2, remove_legend=True, ms_hue=ms_hue, prefix="FDA")
    plot_gaps_ax(runs_vs_ws0, subfig_title="(C)", ax=ax3, remove_legend=True, ms_hue=ms_hue, prefix="VS")
    plot_gaps_ax(runs_vs_ws1, subfig_title="(D)", ax=ax4, remove_legend=True, ms_hue=ms_hue, prefix="VS")

    ax4.legend(loc='lower right', bbox_to_anchor=(1.0, -.45),
               fancybox=True, shadow=False, ncol=2)

    # fig.supylabel("Number of target chemicals")
    fig.tight_layout()
    fig.savefig(figname, dpi=600)
