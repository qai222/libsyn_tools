from __future__ import annotations

import itertools
import os
from collections import defaultdict

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from libsyn_tools.chem_schema import OperationNetwork, FunctionalModule, OperationType
from libsyn_tools.opt.formulation_baseline import SolverBaseline, SchedulerOutput
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.opt.schema import SchedulerInput
from libsyn_tools.utils import json_load, FilePath

plt.rcParams['svg.fonttype'] = 'none'


def plot_gantt(runs_foler: FilePath, run_name: str, save_float_folder: FilePath, anno_reaction_index: bool,
               multi_capacity: bool, figsize=None):
    runs_folder = os.path.abspath(runs_foler)
    run_folder = os.path.join(runs_folder, run_name)
    save_float_folder = os.path.abspath(save_float_folder)
    save_folder = os.path.join(save_float_folder, run_name)
    os.makedirs(save_folder, exist_ok=True)
    solver_baseline = SolverBaseline(**json_load(os.path.join(run_folder, "solver_baseline.json")))
    solver_milp = SolverMILP(**json_load(os.path.join(run_folder, "solver_milp.json")))
    operation_network = OperationNetwork(**json_load(os.path.join(run_folder, "operation_network.json")))
    functional_modules = [FunctionalModule(**fm) for fm in
                          json_load(os.path.join(run_folder, "functional_modules.json"))]

    if figsize is None:
        figsize = (8, 4)
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True, sharey=True)
    ax2.set_xlabel("Time (min)")
    if multi_capacity:
        plot_gantt_ax_multi_capacity(solver_milp, operation_network, functional_modules, ax1,
                                     subfig_title="(A) Optimal schedule")
        plot_gantt_ax_multi_capacity(solver_baseline, operation_network, functional_modules, ax2,
                                     subfig_title="(B) Baseline schedule")
    else:
        plot_gantt_ax(solver_milp, operation_network, ax=ax1, trans=False, anno_reaction_index=anno_reaction_index,
                      subfig_title="(A) Optimal schedule")
        plot_gantt_ax(solver_baseline, operation_network, ax=ax2, trans=False, anno_reaction_index=anno_reaction_index,
                      subfig_title="(B) Baseline schedule")
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
        # print(fm.identifier, fm_new_name)
        if "Heating" in fm_new_name:
            fm_new_name = fm_new_name.replace("Heating", "H")
        elif "ConcentrationAndPurification" in fm_new_name:
            fm_new_name = fm_new_name.replace("ConcentrationAndPurification", "W")
        fm_rename[fm.identifier] = fm_new_name.split(".")[-1]

    df["assigned_to_new_name"] = [fm_rename[n] for n in df['assigned_to'].tolist()]
    df["reaction_identifier"] = [operation_network.operation_dictionary[n].from_reaction for n in
                                 df['operation_identifier'].tolist()]
    df["duration"] = df["end_time"] - df["start_time"]
    # df.to_csv("gantt_chart.csv", index=False)
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

    color_list = mcolors.TABLEAU_COLORS.copy()
    pattern_list = ('', '/', '-', '.', 'x', '+', '\\', '*', 'o', 'O',)
    pattern_list = [p * 3 for p in pattern_list]
    pc_tuples = itertools.product(pattern_list, color_list)
    pc_dict = dict()
    for rid in reaction_ids:
        pc_dict[rid] = next(pc_tuples)

    fms = [fm for fm in fms if fm.startswith("H") or fm.startswith("W")]

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
    ax.set_ylim([0 - 0.5, len(fms) - 0.5 + 0.2])
    ax.set_yticklabels(fms)
    ax.set_title(subfig_title, loc='left')


def overlaps_with(scheduled_operation1: dict, scheduled_operation2: dict) -> bool:
    if scheduled_operation1["assigned_to_module"] != scheduled_operation2["assigned_to_module"]:
        return False
    s1, e1 = scheduled_operation1["start_time"], scheduled_operation1["end_time"]
    s2, e2 = scheduled_operation2["start_time"], scheduled_operation2["end_time"]
    return bool(max(0, min(e1, e2) - max(s1, s2)))


def plot_gantt_ax_multi_capacity(
        solver: SolverMILP | SolverBaseline, operation_network: OperationNetwork,
        functional_modules: list[FunctionalModule],
        ax: plt.Axes, subfig_title: str
):
    bar_height = 0.5
    inter_bar_spacing = bar_height * 2
    inter_fm_spacing = bar_height * 5
    df = get_schedule_df(solver.output, solver.input, operation_network)
    only_show_types = (OperationType.Heating, OperationType.ConcentrationAndPurification)
    fms_shortnames = {
        OperationType.Heating: "H",
        OperationType.ConcentrationAndPurification: "W",
    }

    # heater and workup stations only
    fms = [fm for fm in functional_modules if fm.can_process[0] in only_show_types]
    fms = sorted(fms, key=lambda x: x.can_process[0], reverse=True)
    fms_dict = {fm.identifier: fm for fm in fms}
    fms_y_dict = dict()
    for ifm, fm in enumerate(fms):
        fm_height = fm.capacity * bar_height + inter_bar_spacing * (fm.capacity - 1) + inter_fm_spacing
        if ifm == 0:
            fm_y = 0
        else:
            last_fm = fms[ifm - 1]
            last_fm_y = fms_y_dict[last_fm.identifier]["fm_y"]
            last_fm_height = fms_y_dict[last_fm.identifier]["fm_height"]
            fm_y = last_fm_y + last_fm_height * 0.5 + fm_height * 0.5
        fms_y_dict[fm.identifier] = {
            "fm_y": fm_y, "fm_height": fm_height, "capacity": fm.capacity
        }

    reaction_ids = sorted(df["reaction_identifier"].unique().tolist())
    reaction_indexer = {
        rid: ii for ii, rid in enumerate(reaction_ids)
    }

    color_list = mcolors.TABLEAU_COLORS.copy()
    pattern_list = ('', '/', '-', '.', 'x', '+', '\\', '*', 'o', 'O',)
    pattern_list = [p * 3 for p in pattern_list]
    pc_tuples = itertools.product(pattern_list, color_list)
    pc_dict = dict()
    for rid in reaction_ids:
        pc_dict[rid] = next(pc_tuples)

    scheduled_operations = []
    for record in df.to_dict('records'):
        so = {
            "assigned_to_module": record["assigned_to"],
            "identifier": record["operation_identifier"],
            "start_time": record["start_time"],
            "end_time": record["end_time"],
            "duration": record["duration"],
            "reaction": record["reaction_identifier"],
        }
        reaction_pattern, reaction_color = pc_dict[so["reaction"]]
        so["bar_pattern"] = reaction_pattern
        so["bar_color"] = reaction_color
        scheduled_operations.append(so)

    so_overlap_dict = defaultdict(list)
    for so1 in scheduled_operations:
        for so2 in scheduled_operations:
            if so1["identifier"] == so2["identifier"]:
                continue
            if overlaps_with(so1, so2):
                so_overlap_dict[so1["identifier"]].append(so2)

    already_plotted = []
    ax_y_max = 0
    ax_y_min = 1e3
    for so in scheduled_operations:
        if so["identifier"] in already_plotted:
            continue
        try:
            fm_record = fms_y_dict[so["assigned_to_module"]]
        except KeyError:
            continue
        fm_y = fm_record["fm_y"]
        sos_to_plot = [so, ] + so_overlap_dict[so["identifier"]]
        sos_to_plot = [s for s in sos_to_plot if s["identifier"] not in already_plotted]
        sos_y = []
        for iso, so_to_plot in enumerate(sos_to_plot):
            so_y = bar_height * 0.5 + iso * inter_bar_spacing + iso * bar_height
            sos_y.append(so_y)  # the lower boundary of the first bar is 0
        sos_y_center = sum(sos_y) / len(sos_y)
        sos_y = [so_y_ + fm_y - sos_y_center for so_y_ in sos_y]  # translate center of bars to fm_y

        for iso, so_to_plot in enumerate(sos_to_plot):
            bars = ax.barh(
                y=sos_y[iso],
                width=so_to_plot["duration"],
                height=bar_height,
                left=so_to_plot['start_time'],
                color=so_to_plot["bar_color"],
                linewidth=0.0,
            )
            logger.info(
                f"plotting gantt: reaction={reaction_indexer[so_to_plot['reaction']]} assigned_to={so_to_plot['assigned_to_module']}")
            if sos_y[iso] > ax_y_max:
                ax_y_max = sos_y[iso] + bar_height * 0.5
            if sos_y[iso] < ax_y_min:
                ax_y_min = sos_y[iso] - bar_height * 0.5
            already_plotted.append(so_to_plot["identifier"])
            if so_to_plot["bar_pattern"] != "":
                mpl.rc('hatch', color='k', linewidth=0.5)
                for b in bars:
                    b.set_hatch(so_to_plot["bar_pattern"] * 2)
                    b._hatch_color = mcolors.to_rgb("k")
                    b.stale = True

            text_color = 'black'
            ax.text(
                so_to_plot["start_time"] + so_to_plot["duration"] * 0.5,
                sos_y[iso] + 0.5 * bar_height,
                str(reaction_indexer[so_to_plot["reaction"]]),
                ha='center', va='bottom',
                color=text_color, fontsize=7
            )

    if solver.include_shift_constraints:
        lw = 0.8
        ws_linestyles = (0, (lw, 2))
        font_height = 0.2
        ax_y_min = ax_y_min - font_height
        ax_y_max = ax_y_max + font_height
        for ws in solver.input.frak_W:
            ax.vlines(ws.start_time, ymin=ax_y_min, ymax=ax_y_max, linestyles=ws_linestyles,
                      colors="red", lw=lw)
            ax.vlines(ws.end_time, ymin=ax_y_min, ymax=ax_y_max, linestyles=ws_linestyles,
                      colors="red", lw=lw)
            ax.hlines(y=ax_y_min, xmin=ws.start_time, xmax=ws.end_time, linestyles=ws_linestyles, colors="red",
                      lw=lw)
            ax.hlines(y=ax_y_max, xmin=ws.start_time, xmax=ws.end_time, linestyles=ws_linestyles,
                      colors="red", lw=lw)
            if ws.end_time > df["end_time"].max():
                break

    ax.set_yticks(sorted([r["fm_y"] for r in fms_y_dict.values()]))
    # ax.set_ylim([0 - 0.5, len(fms) - 0.5 + 0.2])
    yticklabels = [fms_shortnames[fm.can_process[0]] for fm in fms]
    yticklabels = [f"{ytl}-{ytlidx}" for ytl, ytlidx in zip(yticklabels, [0, 1, 0, 1])]
    ax.set_yticklabels(yticklabels)
    ax.set_title(subfig_title, loc='left')
