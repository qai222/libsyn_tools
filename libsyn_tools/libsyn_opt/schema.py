import itertools
import math
import random
from typing import Any

import numpy as np
from pydantic import BaseModel

from libsyn_tools.chem_schema import OperationNetwork, FunctionalModule, OperationType
from libsyn_tools.chem_schema.network_operation import default_process_time


class WorkShift(BaseModel):
    start_time: float

    end_time: float

    @property
    def duration(self) -> float:
        return self.end_time - self.start_time


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

    frak_W: list[WorkShift] | None = None
    """ work shifts, indexed by n  """

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
            cls, operation_network: OperationNetwork, functional_modules: list[FunctionalModule],
            human_supervision_types: list[OperationType] = None, work_shifts: list[WorkShift] = None,
    ):
        # check if operable
        operation_types = set([o.type for o in operation_network.operations])
        operable_operation_types = []
        for m in functional_modules:
            operable_operation_types += m.can_process
        assert operation_types.issubset(set(operable_operation_types))

        # calculate process time
        default_process_time(operation_network.operations, functional_modules)

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
        p = np.zeros((n_operations, len(frak_M)))
        p[:] = math.inf
        for operation in operation_network.operations:
            i = frak_O.index(operation.identifier)
            for module_identifier, pt in operation.process_times.items():
                m = frak_M.index(module_identifier)
                p[i][m] = pt

        scheduler_input = cls(
            lmin=lmin.tolist(),
            lmax=lmax.tolist(),
            frak_O=frak_O,
            frak_P=frak_P,
            frak_M=frak_M,
            K=K,
            p=p,
            frak_W=work_shifts,
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

    @classmethod
    def from_MILP(cls, scheduler_input: SchedulerInput, var_s, var_e, var_a):
        output = cls()

        var_s = np.array(var_s)
        var_e = np.array(var_e)
        var_a = np.array(var_a)

        for i, oid in enumerate(scheduler_input.frak_O):
            output.start_times[oid] = var_s[i]
            output.end_times[oid] = var_e[i]
            for m, mid in enumerate(scheduler_input.frak_M):
                if var_a[i][m]:
                    output.assignments[oid] = mid

        output.functional_modules = [*scheduler_input.frak_M]  # make a copy
        return output


class Solver(BaseModel):
    """ actual solver of a specific formulation """

    input: SchedulerInput | None = None

    output: SchedulerOutput | None = None

    time_cost: float | None = None

    opt_settings: dict[str, Any] = dict()

    opt_log: str | None = None

    def solve(self): pass


def random_functional_modules(random_seed: int = 42, max_capacity: int = 3, max_module_number_per_type: int = 3) -> \
list[FunctionalModule]:
    fms = []
    random.seed(random_seed)
    capacity_range = [1, 3]
    module_number_range = [1, 3]
    for t in OperationType:
        module_number = random.randint(*module_number_range)
        for i in range(module_number):
            m = FunctionalModule(
                name=f"{t}-{i}",
                can_process=[t, ],
                capacity=random.randint(*capacity_range),
            )
            fms.append(m)
    return fms
