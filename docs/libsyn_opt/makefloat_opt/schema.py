from __future__ import annotations

import glob
import os

from loguru import logger
from pydantic import BaseModel
from tqdm import tqdm

from libsyn_tools.opt.formulation_baseline import SolverBaseline
from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.utils import FilePath
from libsyn_tools.utils import json_load, json_dump

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
    def is_valid_baseline(self):
        valid = True
        for v in self.validation_baseline.values():
            if v is None:
                continue
            if v is False:
                valid = False
                break
        return valid

    @property
    def is_valid_milp(self):
        valid = True
        for v in self.validation_milp.values():
            if v is None:
                continue
            if v is False:
                valid = False
                break
        return valid

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

            # revalidate just to be safe
            solver_milp.output.notes['validation'] = solver_milp.output.validate_schedule(solver_milp.input)
            validation_milp = solver_milp.output.notes['validation']

            solution_milp = max(solver_milp.output.end_times.values())
            gurobi_status = gurobi_status_dict[solver_milp.opt_log['gurobi status']]
            runtime_milp = solver_milp.opt_log['time solved']
            if int(x) < 8 and runtime_milp > 3600:
                runtime_milp = 3600
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


def split_runs(runs: list[SchedulerRun]):
    runs_fda_ws0 = []
    runs_fda_ws1 = []
    runs_vs_ws0 = []
    runs_vs_ws1 = []

    for r in runs:
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


def load_runs(runs_folder_match: FilePath):
    try:
        runs = [SchedulerRun(**R) for R in json_load("RUNS.json")]
        logger.critical("RUNS.json found in cwd, loaded runs from it")
    except FileNotFoundError:
        logger.critical("RUNS.json not found in cwd, loading from RUNS FOLDER")
        runs = []
        for FOLDER in tqdm(sorted(glob.glob(runs_folder_match))):
            if not os.path.isdir(FOLDER):
                continue
            if not os.path.exists(os.path.join(FOLDER, "solver_baseline.json")):
                continue
            run = SchedulerRun.from_run_folder(FOLDER)
            if run.ifm > 1:
                continue
            runs.append(run)
        json_dump([R.model_dump() for R in runs], "RUNS.json")
    return runs
