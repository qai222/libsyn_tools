from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import pandas as pd
import seaborn as sns

from makefloat_opt.schema import SchedulerRun, split_runs

sns.set_theme()


def plot_status_ax(runs: list[SchedulerRun], subfig_title: str, ax: plt.Axes, prefix:str):
    df = pd.DataFrame.from_records([r.model_dump() for r in runs])
    df['Solver status'] = df['gurobi_status'].fillna(value='No solution')
    p = sns.histplot(df, y="n_target", ax=ax, hue="Solver status", multiple="stack", discrete=True, )
    # ax.set_ylabel("Number of target chemicals")
    ax.set_ylabel("")
    ax.set_xlabel("Scheduler instance count")

    ax.xaxis.set_major_locator(tck.MultipleLocator())
    ax.yaxis.set_major_locator(tck.MultipleLocator(base=1))
    ax.set_xlim([0, 20])

    labels = [item.get_text() for item in ax.get_yticklabels()]
    yticks = ax.get_yticks()
    if "FDA" == prefix:
        labels = [f"FDA.{int(l):02}" for l in labels]
    else:
        labels = [f"VS.{int(l):02}" for l in labels]
    ax.set_yticks(yticks, labels)
    ax.set_ylim([9.5, 0.5])

    # ax.set_title(subfig_title, fontsize='large', loc='center')
    ax.set_title(subfig_title, loc='left')
    return ax


def plot_status(runs: list[SchedulerRun], figname: str):
    runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 10), sharey=True, sharex=True)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey="row", sharex=True)
    plot_status_ax(runs_fda_ws0, subfig_title="(A)", ax=ax1, prefix="FDA")
    plot_status_ax(runs_fda_ws1, subfig_title="(B)", ax=ax2, prefix="FDA")
    plot_status_ax(runs_vs_ws0, subfig_title="(C)", ax=ax3, prefix="VS")
    plot_status_ax(runs_vs_ws1, subfig_title="(D)", ax=ax4, prefix="VS")

    sns.move_legend(ax4, loc='upper center', bbox_to_anchor=(0.5, -.45),
               fancybox=True, shadow=False, ncol=3)
    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()
    fig.tight_layout()
    fig.savefig(figname, dpi=600)
