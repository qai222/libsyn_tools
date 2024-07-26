import glob
import os.path

import pandas as pd

from libsyn_tools.opt.formulation_baseline import SolverBaseline
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.utils import json_load, FilePath


def collect(runs_folder: FilePath):
    records = []
    for run_folder in sorted(glob.glob(f"{runs_folder}/*/")):
        solver_baseline = os.path.join(run_folder, "solver_baseline.json")
        solver_milp = os.path.join(run_folder, "solver_milp.json")
        try:
            solver_baseline = SolverBaseline(**json_load(solver_baseline))
            solver_milp = SolverMILP(**json_load(solver_milp))
        except FileNotFoundError:
            continue

        record = {
            "name": run_folder,
            "baseline solution": max(solver_baseline.output.end_times.values()),
            "baseline valid": all(v for v in solver_baseline.output.notes['validation'].values() if v is not None),
            "milp valid": all(v for v in solver_milp.output.notes['validation'].values() if v is not None),
            "milp runtime": solver_milp.opt_log['time solved'],
            "milp feasible": solver_milp.opt_log['gurobi status'] != 3,
        }

        for k, v in solver_baseline.output.notes['validation'].items():
            record[f"baseline validation: {k}"] = v
        for k, v in solver_milp.output.notes['validation'].items():
            record[f"milp validation: {k}"] = v

        if solver_milp.output.end_times:
            record["milp solution"] = max(solver_milp.output.end_times.values())
            gap = (- record["milp solution"] + record["baseline solution"]) / record["baseline solution"]
            record["gap"] = f"{gap: .1%}"
        else:
            record["milp solution"] = None
            record["gap"] = None

        records.append(record)
    return pd.DataFrame.from_records(records)


if __name__ == '__main__':
    RUNS_FOLDER = os.path.join(os.getcwd(), "RUNS")
    DF = collect(RUNS_FOLDER)
    DF.to_csv("collect_libs.csv", index=False)
