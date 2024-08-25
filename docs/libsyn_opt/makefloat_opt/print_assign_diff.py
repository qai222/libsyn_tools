from __future__ import annotations

import os
from collections import defaultdict

import Levenshtein
import numpy as np
import pandas as pd

from libsyn_tools.chem_schema import OperationNetwork
from libsyn_tools.opt.formulation_baseline import SolverBaseline, SchedulerOutput
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.opt.schema import OperationType
from libsyn_tools.utils import json_load, FilePath
from makefloat_opt.schema import SchedulerRun


def get_sequences(output: SchedulerOutput):
    seqs = defaultdict(list)
    for oid, mid in output.assignments.items():
        seqs[mid].append(oid)
    for mid in seqs.keys():
        seqs[mid] = sorted(seqs[mid], key=lambda x: output.start_times[x])
    return seqs


def print_assign_diff_all(runs: list[SchedulerRun], save_float_folder: FilePath):
    records = []

    for run in runs:
        run_folder = run.folder
        if run.n_target < 4:
            continue
        if run.percentage_gap > 3:
            continue
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
    df.to_csv("assign_diff.csv", index=False)
    print(df['percentage_n_assign'].mean())
    # print(df['d_hamming_mean'].mean())
    print(df['d_lev_mean'].mean(), df['d_lev_mean'].min(), df['d_lev_mean'].max())
