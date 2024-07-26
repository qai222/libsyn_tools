from __future__ import annotations

import itertools
import math
import random
from typing import Any

import numpy as np
from pydantic import BaseModel

from libsyn_tools.chem_schema import OperationNetwork, FunctionalModule, OperationType
from libsyn_tools.chem_schema.network_operation import default_process_time
from loguru import logger

class WorkShift(BaseModel):
    start_time: float

    end_time: float

    @property
    def duration(self) -> float:
        return self.end_time - self.start_time

    @staticmethod
    def get_regular_shifts(duration: float, interval: float, horizon: float) -> list[WorkShift]:
        shifts = []
        clock = 0
        while clock <= horizon:
            shift = WorkShift(start_time=clock, end_time=clock + duration)
            shifts.append(shift)
            clock += duration + interval
        return shifts


class SchedulerInput(BaseModel):
    """ input for FJSS """

    lmin: list[list[float]]
    """ integer-indexed min lag constraints, indexed by i, j """

    lmax: list[list[float]]
    """ integer-indexed max lag constraints, indexed by i, j """

    frak_O: list[str]
    """ operation identifiers, indexed by i """

    frak_P: list[str]
    """ subset of `frak_O` that require human supervision (thus subject to work shift constraints) """

    frak_M: list[str]
    """ functional module identifiers, indexed by m """

    K: list[int]
    """ capacity, indexed by m """

    p: list[list[float]]
    """ process time, indexed by i, m """

    C: list[list[int]]
    """ compatability on the same machine, indexed by i, m """

    frak_W: list[WorkShift] | None = None
    """ work shifts, indexed by n  """

    functional_modules: list[FunctionalModule]

    def get_dummy_work_shifts(self):
        p_finite = np.array(self.p, dtype=float)
        p_finite = p_finite[p_finite < math.inf]
        p_max = p_finite.max()
        duration = p_max * 1.5
        interval = duration * 0.5
        horizon = p_max * len(self.p)
        work_shifts = WorkShift.get_regular_shifts(duration, interval, horizon)
        return work_shifts

    @property
    def S(self) -> list[float] | None:
        """ work shift start time, indexed by n """
        if not self.frak_W:
            return
        return [w.start_time for w in self.frak_W]

    @property
    def E(self) -> list[float] | None:
        """ work shift end time, indexed by n """
        if not self.frak_W:
            return
        return [w.end_time for w in self.frak_W]

    @property
    def summary(self):
        d = dict()
        d['# of operations'] = len(self.frak_O)
        d['# of functional modules'] = len(self.frak_M)
        d['# of work shifts'] = None if self.frak_W is None else len(self.frak_W)

        n_precedence = 0
        for i, j in itertools.combinations(range(len(self.frak_O)), 2):
            if self.lmin[i][j] > - math.inf or self.lmax[i][j] < math.inf:
                n_precedence += 1
        d['# of precedence relations'] = n_precedence
        return d

    @classmethod
    def build_input(
            cls,
            rng: random.Random,
            operation_network: OperationNetwork, functional_modules: list[FunctionalModule],
            # temperature_threshold: float,
            human_supervision_types: list[OperationType] = None, work_shifts: list[WorkShift] = None,
    ):
        # check if operable
        operation_types = set([o.type for o in operation_network.operations])
        operable_operation_types = []
        for m in functional_modules:
            operable_operation_types += m.can_process
        assert operation_types.issubset(set(operable_operation_types))

        # calculate process time
        default_process_time(operation_network.operations, functional_modules, rng)

        # operation helper
        n_operations = len(operation_network.operations)
        g = operation_network.nx_digraph

        frak_O = []
        frak_P = []
        for operation in operation_network.operations:
            frak_O.append(operation.identifier)
            if human_supervision_types is None or operation.type in human_supervision_types:
                frak_P.append(operation.identifier)

        # get lmin and lmax by integer index
        lmin = np.zeros((n_operations, n_operations))
        lmin[:] = - math.inf  # no precedence

        lmax = np.zeros((n_operations, n_operations))
        lmax[:] = math.inf
        for u, v, d in g.edges(data=True):
            try:
                _lmin = d['lmin']
            except KeyError:
                _lmin = 0
            if lmin is None:
                _lmin = 0

            try:
                _lmax = d['lmax']
            except KeyError:
                _lmax = math.inf
            if lmax is None:
                _lmax = math.inf
            uid = frak_O.index(u)
            vid = frak_O.index(v)

            lmax[uid][vid] = _lmax
            lmin[uid][vid] = _lmin

        # functional module params
        frak_M = []
        K = []
        for m in functional_modules:
            frak_M.append(m.identifier)
            K.append(m.capacity)

        # process time
        p = np.zeros((n_operations, len(frak_M)), dtype=float)
        p[:] = math.inf
        for operation in operation_network.operations:
            i = frak_O.index(operation.identifier)
            for module_identifier, pt in operation.process_times.items():
                m = frak_M.index(module_identifier)
                p[i][m] = pt

        # compatability
        C = np.zeros((n_operations, n_operations), dtype=int)
        # compatability = operation_network.get_default_compatability(temperature_threshold=temperature_threshold)
        compatability = operation_network.get_default_compatability_ad_hoc()
        for i, oid in enumerate(frak_O):
            for j, ojd in enumerate(frak_O):
                C[i][j] = compatability[oid][ojd]

        scheduler_input = cls(
            lmin=lmin.tolist(),
            lmax=lmax.tolist(),
            frak_O=frak_O,
            frak_P=frak_P,
            frak_M=frak_M,
            K=K,
            p=p.tolist(),
            frak_W=work_shifts,
            functional_modules=functional_modules,
            C=C.tolist(),
        )
        return scheduler_input


class SchedulerOutput(BaseModel):
    """ output for FJSS """

    start_times: dict[str, float] = dict()
    """ operation identifier to start time """

    end_times: dict[str, float] = dict()
    """ operation identifier to end time """

    assignments: dict[str, str] = dict()
    """ operation identifier to functional module identifier """

    functional_modules: list[str] = []
    """ functional module identifiers, consistent with input, note this is a subset of the range of `assignments` """

    notes: dict = dict()
    """ additional notes """

    def operation_view(self):
        view = []
        for i, o in enumerate(self.start_times):
            d = dict(
                operation_identifier=o,
                milp_index=i,
                assigned_to=self.assignments[o],
                start_time=self.start_times[o],
                end_time=self.end_times[o],
            )
            view.append(d)
        return view

    @classmethod
    def from_MILP(cls, scheduler_input: SchedulerInput, var_s: np.ndarray, var_e: np.ndarray, var_a: np.ndarray):
        output = cls()

        for i, oid in enumerate(scheduler_input.frak_O):
            output.start_times[oid] = var_s[i]
            output.end_times[oid] = var_e[i]
            for m, mid in enumerate(scheduler_input.frak_M):
                if var_a[i][m]:
                    output.assignments[oid] = mid

        output.functional_modules = [*scheduler_input.frak_M]  # make a copy
        return output

    def validate_schedule(self, scheduler_input: SchedulerInput, eps=1e-3, consider_work_shifts: bool = True):
        report = {
            "ordinary precedence valid": [],
            "lmax only precedence valid": [],
            "lmin only precedence valid": [],
            "lmin and lmax precedence valid": [],
            "work shift valid": None,
        }
        size_i = len(scheduler_input.frak_O)

        if scheduler_input.frak_W and consider_work_shifts:
            report["work shift valid"] = []
            for i in range(size_i):
                if scheduler_input.frak_O[i] not in scheduler_input.frak_P:
                    continue
                assigned_n = None
                for n in range(len(scheduler_input.frak_W)):
                    if scheduler_input.S[n] <= self.start_times[scheduler_input.frak_O[i]] <= self.end_times[
                        scheduler_input.frak_O[i]] <= scheduler_input.E[n]:
                        assigned_n = n
                if assigned_n is None:
                    report["work shift valid"].append(False)
                else:
                    report["work shift valid"].append(True)

        for i in range(size_i):
            for j in range(size_i):
                lmax = scheduler_input.lmax[i][j]
                lmin = scheduler_input.lmin[i][j]
                s_j = self.start_times[scheduler_input.frak_O[j]]
                e_i = self.end_times[scheduler_input.frak_O[i]]
                if lmax == math.inf and lmin == - math.inf:
                    continue
                elif lmin == 0 and lmax == math.inf:
                    valid = bool(s_j > e_i - eps)
                    report['ordinary precedence valid'].append(valid)
                elif lmin > - math.inf and lmax == math.inf:
                    valid = bool(s_j > e_i + lmin - eps)
                    report['lmin only precedence valid'].append(valid)
                elif lmin == - math.inf and lmax < math.inf:
                    valid = bool(s_j < e_i + lmax + eps)
                    report['lmax only precedence valid'].append(valid)
                elif lmin > - math.inf and lmax < math.inf:
                    valid = bool(e_i + lmax + eps > s_j > e_i + lmin - eps)
                    report['lmin and lmax precedence valid'].append(valid)
                else:
                    raise ValueError
        simple_report = dict()
        for k, v in report.items():
            if isinstance(v, list):
                simple_report[k] = all(v)
            else:
                simple_report[k] = v
        return simple_report


class Solver(BaseModel):
    """ actual solver of a specific formulation """

    input: SchedulerInput | None = None

    output: SchedulerOutput | None = None

    time_cost: float | None = None

    time_limit: float | None = None

    opt_settings: dict[str, Any] = dict()

    opt_log: dict = dict()

    def solve(self): pass


def random_functional_modules(rng: random.Random, max_capacity: int = 3, max_module_number_per_type: int = 3) -> \
        list[FunctionalModule]:
    fms = []
    capacity_range = [1, max_capacity]
    module_number_range = [1, max_module_number_per_type]
    for t in OperationType:
        module_number = rng.randint(*module_number_range)
        for i in range(module_number):
            m = FunctionalModule(
                name=f"{t}-{i}",
                can_process=[t, ],
                capacity=rng.randint(*capacity_range),
            )
            fms.append(m)
    return fms
