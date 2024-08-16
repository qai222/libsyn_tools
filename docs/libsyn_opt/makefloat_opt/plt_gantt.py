from __future__ import annotations

import itertools
import os
from collections import defaultdict

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from libsyn_tools.chem_schema import OperationNetwork
from libsyn_tools.opt.formulation_baseline import SolverBaseline, SchedulerOutput
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.opt.schema import SchedulerInput
from libsyn_tools.utils import json_load, FilePath

sns.set_theme()


# # RUN_NAME = "FDA-03-09-0-0"
# # RUN_NAME = "FDA-03-09-1-0"
# # RUN_NAME = "VS-07-10-1-1"
# RUN_NAME = "VS-04-06-0-1"
# RUNS_FOLDER = os.path.abspath("../workplace/RUNS")


def plot_gantt(runs_foler: FilePath, run_name: str, save_float_folder: FilePath, anno_reaction_index: bool):
    runs_folder = os.path.abspath(runs_foler)
    run_folder = os.path.join(runs_folder, run_name)
    save_float_folder = os.path.abspath(save_float_folder)
    save_folder = os.path.join(save_float_folder, run_name)
    os.makedirs(save_folder, exist_ok=True)
    solver_baseline = SolverBaseline(**json_load(os.path.join(run_folder, "solver_baseline.json")))
    solver_milp = SolverMILP(**json_load(os.path.join(run_folder, "solver_milp.json")))
    operation_network = OperationNetwork(**json_load(os.path.join(run_folder, "operation_network.json")))

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8, 4), sharex=True, sharey=True)
    ax2.set_xlabel("Time (min)")
    plot_gantt_ax(solver_milp, operation_network, ax=ax1, trans=False, anno_reaction_index=anno_reaction_index, subfig_title="(A) Optimal schedule")
    plot_gantt_ax(solver_baseline, operation_network, ax=ax2, trans=False, anno_reaction_index=anno_reaction_index, subfig_title="(B) Baseline schedule")
    fig.tight_layout()
    fig.savefig(os.path.join(save_folder, f"gantt-{run_name}.pdf"))

    return solver_baseline, solver_milp, operation_network


def get_schedule_df(scheduler_output: SchedulerOutput, scheduler_input: SchedulerInput,
                    operation_network: OperationNetwork):
    records = scheduler_output.operation_view()
    df = pd.DataFrame.from_records(records)
    fm_groups = defaultdict(list)
    for fm in scheduler_input.functional_modules:
        fm_groups[fm.can_process[0]].append(fm.identifier)

    fm_rename = dict()
    for fm in scheduler_input.functional_modules:
        fm_new_name = f"{fm.can_process[0]}-{fm_groups[fm.can_process[0]].index(fm.identifier)}"
        print(fm.identifier, fm_new_name)
        if "Heating" in fm_new_name:
            fm_new_name = fm_new_name.replace("Heating", "Heater")
        elif "ConcentrationAndPurification" in fm_new_name:
            fm_new_name = fm_new_name.replace("ConcentrationAndPurification", "Workup station")
        fm_rename[fm.identifier] = fm_new_name.split(".")[-1]

    df["assigned_to_new_name"] = [fm_rename[n] for n in df['assigned_to'].tolist()]
    df["reaction_identifier"] = [operation_network.operation_dictionary[n].from_reaction for n in
                                 df['operation_identifier'].tolist()]
    df["duration"] = df["end_time"] - df["start_time"]
    df.to_csv("gantt_chart.csv", index=False)
    return df


def plot_gantt_ax(solver: SolverMILP | SolverBaseline, operation_network: OperationNetwork, ax: plt.Axes, trans: bool,
                  anno_reaction_index: bool,
                  subfig_title: str):
    height = 0.5
    df = get_schedule_df(solver.output, solver.input, operation_network)
    fms = sorted(df['assigned_to_new_name'].unique())
    reaction_ids = sorted(df["reaction_identifier"].unique().tolist())
    reaction_indexer = {
        rid: ii for ii, rid in enumerate(reaction_ids)
    }
    for rid, ii in reaction_indexer.items():
        print(ii, rid)

    color_list = mcolors.TABLEAU_COLORS.copy()
    pattern_list = ('', '/', '-', '.', 'x', '+', '\\', '*', 'o', 'O',)
    pattern_list = [p * 3 for p in pattern_list]
    pc_tuples = itertools.product(pattern_list, color_list)
    pc_dict = dict()
    for rid in reaction_ids:
        pc_dict[rid] = next(pc_tuples)

    fms = [fm for fm in fms if fm.startswith("H") or fm.startswith("Workup")]

    for reaction_id, operations in operation_network.operations_by_reaction.items():
        reaction_df = df[df['reaction_identifier'] == reaction_id].copy()
        reaction_df = reaction_df[reaction_df["assigned_to_new_name"].isin(fms)]
        reaction_df["gantt_y"] = [fms.index(fm) for fm in reaction_df['assigned_to_new_name'].tolist()]
        reaction_pattern, reaction_color = pc_dict[reaction_id]
        reaction_color_rgba = mcolors.to_rgba(reaction_color)
        if trans:
            reaction_color_rgba = list(reaction_color_rgba)
            reaction_color_rgba[-1] = 0.4
            reaction_color_rgba = tuple(reaction_color_rgba)
        bars = ax.barh(
            y=reaction_df["gantt_y"],
            width=reaction_df['duration'],
            height=height,
            left=reaction_df['start_time'],
            linewidth=0.0,
            color=reaction_color_rgba,
        )
        if reaction_pattern != "":
            for b in bars:
                b.set_hatch(reaction_pattern)
                b._hatch_color = mcolors.to_rgb("k")
                b.stale = True

        if anno_reaction_index:
            xcenters = reaction_df['start_time'].values + reaction_df['duration'].values / 2
            # text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
            # text_ys = reaction_df["gantt_y"]
            # for text_y, text_x in zip(text_ys, xcenters):
            #     ax.text(text_x, text_y, reaction_indexer[reaction_id], ha='center', va='center', color=text_color)
            text_color = 'black'
            text_ys = reaction_df["gantt_y"]
            for text_y, text_x in zip(text_ys, xcenters):
                ax.text(text_x, text_y + 0.5 * height, str(reaction_indexer[reaction_id]), ha='center', va='bottom',
                        color=text_color, fontsize=7)

    if solver.include_shift_constraints:
        for ws in solver.input.frak_W:
            ax.vlines(ws.start_time, ymin=0 - height * 0.5, ymax=len(fms) - height * 1.5, ls=":", colors="gray")
            ax.vlines(ws.end_time, ymin=0 - height * 0.5, ymax=len(fms) - height * 1.5, ls=":", colors="gray")
            ax.hlines(y=0 - height * 0.5, xmin=ws.start_time, xmax=ws.end_time, ls=":", colors="gray")
            ax.hlines(y=len(fms) - height * 1.5, xmin=ws.start_time, xmax=ws.end_time, ls=":", colors="gray")
            if ws.end_time > df["end_time"].max():
                break

    ax.set_yticks(range(len(fms)))
    ax.set_ylim([0-0.5, len(fms)-0.5+0.2])
    ax.set_yticklabels(fms)
    ax.set_title(subfig_title, loc='left')
