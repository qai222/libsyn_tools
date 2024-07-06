import itertools
import math
import pprint
from collections import defaultdict
from typing import Any

import numpy as np
from loguru import logger

from .schema import Solver, SchedulerOutput


class SolverBaseline(Solver):

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

        occupied_until = {m: 0.0 for m in range(size_m)}
        assigned = []
        var_a_im = np.zeros((size_i, size_m))
        var_s_i = np.zeros(size_i)
        var_e_i = np.zeros(size_i)

        while len(assigned) < size_i:
            # get unassigned operations
            unassigned = set(range(size_i)).difference(set(assigned))

            # do not assign i if its precedents are not already assigned
            can_be_assigned = [i for i in unassigned if set(i_to_precedents[i]).issubset(set(assigned))]

            # prioritize "can be assigned" operations based on possible end times
            possible_end_times = []
            for i in can_be_assigned:
                precedent_end_time = 0 if len(i_to_precedents[i]) == 0 else max(var_e_i[j] for j in i_to_precedents[i])
                selected_m = None
                selected_possible_end_time = math.inf
                selected_possible_start_time = None
                for m in i_to_possible_ms[i]:
                    m_is_occupied_until = occupied_until[m]
                    processing_time = self.input.p[i][m]
                    min_start_time = max(m_is_occupied_until, precedent_end_time)
                    possible_end_time = min_start_time + processing_time
                    if possible_end_time < selected_possible_end_time:
                        selected_m = m
                        selected_possible_end_time = possible_end_time
                        selected_possible_start_time = min_start_time
                possible_end_times.append((i, selected_m, selected_possible_end_time, selected_possible_start_time))
            to_be_assigned = sorted(possible_end_times, key=lambda x: x[2])[0]
            to_be_assigned_i, to_be_assigned_m, to_be_assigned_end_time, to_be_assigned_start_time = to_be_assigned

            # assignment
            var_s_i[to_be_assigned_i] = to_be_assigned_start_time
            var_e_i[to_be_assigned_i] = to_be_assigned_end_time
            occupied_until[to_be_assigned_m] = to_be_assigned_end_time
            var_a_im[to_be_assigned_i][to_be_assigned_m] = 1
            assigned.append(to_be_assigned_i)

        self.output = SchedulerOutput.from_MILP(self.input, var_s_i, var_e_i, var_a_im)
