import itertools
import math
import pprint
from collections import defaultdict
from typing import Any

import numpy as np
from loguru import logger

from .schema import Solver, SchedulerOutput


class SolverBaseline(Solver):
    consider_shifts: bool = True

    @property
    def include_shift_constraints(self):
        # TODO DRY
        return self.consider_shifts and self.input.frak_W

    def model_post_init(self, __context: Any) -> None:
        logger.info("\n" + pprint.pformat(self.input.summary))

    def solve(self):
        size_i = len(self.input.frak_O)
        size_m = len(self.input.frak_M)

        i_to_precedents = defaultdict(list)  # include both ordinary and specific
        for i, j in itertools.product(range(size_i), range(size_i)):
            if self.input.lmin[i][j] >= 0 or self.input.lmax[i][j] < math.inf:
                i_to_precedents[j].append(i)

        i_to_possible_ms = defaultdict(list)
        for i in range(size_i):
            i_to_possible_ms[i] = [m for m in range(size_m) if self.input.p[i][m] < math.inf]

        occupied_until = defaultdict(dict)
        for m in range(size_m):
            for slot in range(self.input.K[m]):
                occupied_until[m][slot] = 0.0

        assigned = []
        var_a_im = np.zeros((size_i, size_m))
        var_s_i = np.zeros(size_i)
        var_e_i = np.zeros(size_i)

        var_a_msi = defaultdict(lambda: defaultdict(dict))
        var_a_ims = defaultdict(lambda: defaultdict(dict))
        for i in range(size_i):
            for m in range(size_m):
                for k in range(self.input.K[m]):
                    var_a_msi[m][k][i] = 0
                    var_a_ims[i][m][k] = 0

        while len(assigned) < size_i:
            # get unassigned operations
            unassigned = set(range(size_i)).difference(set(assigned))

            # do not assign i if its precedents are not already assigned
            can_be_assigned = [i for i in unassigned if set(i_to_precedents[i]).issubset(set(assigned))]

            # prioritize "can be assigned" operations based on possible end times
            assignment_candidates = []
            for i in can_be_assigned:
                # find latest precedent + lmin, note lmax is not considered in this algorithm
                latest_precedent_end_time = 0
                if len(i_to_precedents[i]):
                    for j in i_to_precedents[i]:
                        e_j = var_e_i[j]
                        if e_j + self.input.lmin[j][i] > latest_precedent_end_time:
                            latest_precedent_end_time = e_j + self.input.lmin[j][i]

                # greedy assignment for this candidate
                possible_assignment_for_i = []
                for m in i_to_possible_ms[i]:
                    processing_time = self.input.p[i][m]
                    # find the latest incompatible job that is assigned to the same m
                    incompatible_i = None
                    for _incompatible_i in range(size_i):
                        if var_a_im[_incompatible_i][m] == 1 and not self.input.C[_incompatible_i][i]:
                            incompatible_i = _incompatible_i
                    if incompatible_i is None:
                        incompatible_i_end_time = 0
                    else:
                        incompatible_i_end_time = var_e_i[incompatible_i]

                    for slot in range(self.input.K[m]):
                        # find possible start time
                        possible_start_times = [
                            occupied_until[m][slot],  # cannot start if the slot is occupied
                            latest_precedent_end_time,  # cannot start if precedents are not finished
                            incompatible_i_end_time,  # cannot start if incompatible are not finished
                        ]
                        possible_start_time = max(possible_start_times)

                        if self.include_shift_constraints:
                            # TODO need to validate the work shift list is ordered
                            size_n = len(self.input.frak_W)
                            if self.input.frak_O[i] in self.input.frak_P:
                                n_assigned = None
                                for n in range(size_n):
                                    S_n = self.input.S[n]
                                    E_n = self.input.E[n]
                                    if S_n <= possible_start_time <= E_n and possible_start_time + processing_time <= E_n:
                                        n_assigned = n
                                        break
                                    elif S_n >= possible_start_time:
                                        n_assigned = n
                                        possible_start_time = S_n
                                        assert E_n >= S_n + processing_time, f"the processing time of {i} on {m} is longer than a shift!"
                                        break
                                assert n_assigned is not None, "shift not found!"

                        possible_end_time = possible_start_time + processing_time
                        possible_assignment_for_i.append(
                            [i, possible_start_time, possible_end_time, m, slot]
                        )
                possible_assignment_for_i = sorted(possible_assignment_for_i, key=lambda x: x[2])[0]
                assignment_candidates.append(possible_assignment_for_i)

            assignment_candidate = sorted(assignment_candidates, key=lambda x: x[2])[0]
            assign_i, assign_s_i, assign_e_i, assign_m, assign_slot = assignment_candidate

            # assignment
            var_s_i[assign_i] = assign_s_i
            var_e_i[assign_i] = assign_e_i
            occupied_until[assign_m][assign_slot] = assign_e_i
            var_a_im[assign_i][assign_m] = 1
            assigned.append(assign_i)

        self.output = SchedulerOutput.from_MILP(self.input, var_s_i, var_e_i, var_a_im)
