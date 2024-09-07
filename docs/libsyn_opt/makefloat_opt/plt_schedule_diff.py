from __future__ import annotations

import os
import string
from collections import defaultdict

import Levenshtein
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from scipy import stats

from libsyn_tools.chem_schema import OperationNetwork
from libsyn_tools.opt.formulation_baseline import SolverBaseline, SchedulerOutput
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.opt.schema import OperationType
from libsyn_tools.utils import json_load
from makefloat_opt.schema import SchedulerRun

sns.set_theme()


def get_sequences(output: SchedulerOutput):
    seqs = defaultdict(list)
    for oid, mid in output.assignments.items():
        seqs[mid].append(oid)
    for mid in seqs.keys():
        seqs[mid] = sorted(seqs[mid], key=lambda x: output.start_times[x])
    return seqs


def calculate_diffs(runs: list[SchedulerRun]):
    records = []

    for run in runs:
        run_folder = run.folder
        # if run.n_target < 4:
        #     continue
        # if run.percentage_gap > 3:
        #     continue
        n_diff_assign = 0
        n_diff_assign_selected_modules = 0
        n_assign = 0
        n_assign_selected_modules = 0
        solver_baseline = SolverBaseline(**json_load(os.path.join(run_folder, "solver_baseline.json")))
        solver_milp = SolverMILP(**json_load(os.path.join(run_folder, "solver_milp.json")))
        operation_network = OperationNetwork(**json_load(os.path.join(run_folder, "operation_network.json")))

        for oid in solver_baseline.output.assignments:
            baseline_assign = solver_baseline.output.assignments[oid]
            milp_assign = solver_milp.output.assignments[oid]
            n_assign += 1

            if operation_network.operation_dictionary[oid].type in [OperationType.ConcentrationAndPurification,
                                                                    OperationType.Heating]:
                n_assign_selected_modules += 1

            if baseline_assign != milp_assign:
                n_diff_assign += 1
                if operation_network.operation_dictionary[oid].type in [OperationType.ConcentrationAndPurification,
                                                                        OperationType.Heating]:
                    n_diff_assign_selected_modules += 1

        seqs_baseline = get_sequences(solver_baseline.output)
        seqs_milp = get_sequences(solver_milp.output)

        d_levs = []
        # d_hammings = []
        for mid in seqs_baseline:
            seq_baseline = seqs_baseline[mid]
            seq_milp = seqs_milp[mid]
            # d_hamming = hamming(seq_baseline, seq_milp)
            d_lev = Levenshtein.ratio(seq_baseline, seq_milp)
            d_levs.append(d_lev)
            # d_hammings.append(d_hamming)
        prefix, n_target, n_index, ifm, ws = run.name.split("-")
        record = {
            "Library name": ".".join([prefix, n_target, n_index]),
            "Module set": "LAB-1" if int(ifm) == 0 else "LAB-2",
            "Work shift": True if int(ws) == 1 else False,
            "n_diff_assign": n_diff_assign,
            "n_assign": n_assign,
            "n_diff_assign_selected_modules": n_diff_assign_selected_modules,
            "n_assign_selected_modules": n_assign_selected_modules,
            "percentage_gap": run.percentage_gap,
            "targets": run.n_target,
            "percentage_n_assign_selected_modules": round(
                n_diff_assign_selected_modules / n_assign_selected_modules * 100, 1),
            "percentage_n_assign": round(n_diff_assign / n_assign * 100, 1),
            # "d_hamming_mean": np.array(d_hamming).mean(),
            "d_lev_mean": round(np.array(d_levs).mean(), 2),
        }
        records.append(record)
    df = pd.DataFrame.from_records(records)
    # print("n assignment percentage diff", df['percentage_n_assign'].mean())
    # # print(df['d_hamming_mean'].mean())
    # print("lev distance", df['d_lev_mean'].mean(), df['d_lev_mean'].min(), df['d_lev_mean'].max())
    return df


def plot_diffs(df, figname):
    df["Sequence difference"] = [1 - v for v in df["d_lev_mean"]]

    all_percentage_n_assign = df["percentage_n_assign"]
    all_sequence_difference = df["Sequence difference"]

    logger.critical(f"""
    global DIFF
    percentage_n_assign: {all_percentage_n_assign.min()} - {all_percentage_n_assign.max()}, mean: {all_percentage_n_assign.mean()}
    Sequence_difference: {all_sequence_difference.min()} - {all_sequence_difference.max()}, mean: {all_sequence_difference.mean()}
    """)

    # # df = df[df["targets"] > 2]
    # print(.min(), df["percentage_n_assign"].max(), df["percentage_n_assign"].mean())
    # print(df["Sequence difference"].min(), df["Sequence difference"].max(), df["Sequence difference"].mean())
    target_max = 7

    fig, axs = plt.subplots(target_max // 2 + 1, 2, figsize=(5 * 2, 3 * (target_max // 2 + 1)), sharex=True,
                            sharey=True)
    ax_inv = axs[-1][-1]
    ax_inv.xaxis.set_visible(False)
    ax_inv.yaxis.set_visible(False)
    ax_inv.set_axis_off()

    for target in range(1, target_max + 1):
        df_sub = df[df["targets"] == target]
        ax = axs[(target - 1) // 2, (target - 1) % 2]
        sp = sns.scatterplot(
            x="percentage_n_assign",
            y="Sequence difference",
            hue="percentage_gap",
            # size="percentage_gap",
            data=df_sub, ax=ax,
            legend="brief"
        )
        logger.critical(
            f"""
            >>>> DIFF STATS Target # {target} <<<<
            {get_corr_string(df_sub)}
            """
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        # ax.set_title(f"({string.ascii_lowercase[target - 1].upper()}) Target number = {target}", loc='left')
        ax.set_title(f"({string.ascii_lowercase[target - 1].upper()})", loc='left')
        sp.legend_.set_title("Gap (%)")
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1.0, 1.0), ncol=1)
    # sns.scatterplot(x="percentage_n_assign_selected_modules", y="percentage_gap", data=df, ax=ax1)
    # sns.scatterplot(x="d_lev_mean", y="percentage_gap", data=df, ax=ax2)
    # ax2.set_xlabel("Processing sequence similarity")
    fig.supxlabel("Assignment difference (%)", x=0.5)
    fig.supylabel("Sequence difference")
    fig.tight_layout()
    fig.savefig(figname)
    return df


def get_corr_string(df):
    s = f"""
    percentage_n_assign: {df["percentage_n_assign"].min()} - {df["percentage_n_assign"].max()}, mean: {df["percentage_n_assign"].mean()}
    Sequence_difference: {df["Sequence difference"].min()} - {df["Sequence difference"].max()}, mean: {df["Sequence difference"].mean()}
    """
    for x, y in [
        ["percentage_n_assign", "percentage_gap"],
        ["Sequence difference", "percentage_gap"],
        ["percentage_n_assign", "Sequence difference"],
    ]:
        correlation, p_value = stats.pearsonr(df[x], df[y])
        s += f"""between {x} and {y} -- correlation: {correlation} p_value: {p_value}"""
    return s


def plot_diffs_from_runs(runs: list[SchedulerRun], figname):
    return plot_diffs(calculate_diffs(runs), figname)
