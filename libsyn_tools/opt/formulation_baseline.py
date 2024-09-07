from __future__ import annotations

import copy
import itertools
import math
import pprint
from collections import defaultdict

import numpy as np
from loguru import logger
from pydantic import BaseModel

from .schema import Solver, SchedulerOutput, SchedulerInput
from ..chem_schema import ReactionNetwork, OperationNetwork, OperationType


class GreedyAssignment(BaseModel):
    """ an assignment in one iteration of the greedy algorithm """

    operation_index: int
    """ which operation is this assignment about? indexed by i """

    start_time: float
    """ the start time of the assignment """

    end_time: float
    """ the end time of the assignment """

    assigned_to_module: int
    """ to which module is this assignment, indexed by m """

    assigned_to_slot: int
    """ which slot on M_m is this assignment, indexed by the slot index of a specific m """

    has_lmax_subsequent: bool
    """ does this assignment have lmax subsequent """

    possible_start_times: list[float]
    """ possible start times, used for tracking how the assignment is generated """

    latest_precedent: int | None
    """ what is the latest precedent of the assigned operation """

    considers_work_shift: bool
    """ is work shift considered """

    def add_clock_zero(self, clock_zero, scheduler_input: SchedulerInput):
        """
        add global time shift to this assignment
        if work shift is considered, the resulting start time would be shifted further down
        such that it falls in one of the work shifts

        :param clock_zero: time shift value
        :param scheduler_input:
        :return:
        """
        self.start_time += clock_zero
        self.end_time += clock_zero
        if self.considers_work_shift:
            in_a_work_shift = False
            for ws_start, ws_end in zip(scheduler_input.S, scheduler_input.E):
                if ws_start <= self.start_time < self.end_time <= ws_end:
                    in_a_work_shift = True
            if not in_a_work_shift:
                for ws_start in scheduler_input.S:
                    if ws_start > self.start_time:
                        duration = self.end_time - self.start_time
                        self.start_time = ws_start
                        self.end_time = self.start_time + duration
                        break

    @property
    def priority_score(self):
        """ always prioritize the assignment that has a lmax subsequent """
        return self.end_time + int(self.has_lmax_subsequent) * 1e5

    def __gt__(self, other: GreedyAssignment):
        return self.priority_score > other.priority_score

    def __lt__(self, other: GreedyAssignment):
        return self.priority_score < other.priority_score

    def __eq__(self, other: GreedyAssignment):
        return self.priority_score == other.priority_score


def get_greedy_assignment(
        i: int, occupied_until: dict, assigned, var_e_i, var_a_im,
        i_to_precedents: dict, i_to_possible_ms: dict,
        i_to_subsequents_lmax: dict,
        scheduler_input: SchedulerInput,
        consider_work_shift: bool,
) -> GreedyAssignment:
    """
    greedy iteration

    :param i: the operation for which the assignment is to be generated
    :param occupied_until: d[m][slot_index] -> until when the module slot is occupied
    :param assigned: already assigned operations
    :param var_e_i: e_i from the formulation
    :param var_a_im: a_im from the formulation
    :param i_to_precedents: dict[i] -> list of precedents
    :param i_to_possible_ms: dict[i] -> list of possible modules
    :param i_to_subsequents_lmax: dict[i] -> list of subsequents of i that have a lmax
    :param scheduler_input:
    :param consider_work_shift:
    :return:
    """
    # start time should be later than the last precedent plus lmin
    precedents = i_to_precedents[i]
    latest_precedent_end_time = 0
    latest_precedent = None
    if precedents:
        for j in precedents:
            e_j = var_e_i[j]
            if e_j + scheduler_input.lmin[j][i] > latest_precedent_end_time:
                latest_precedent_end_time = e_j + scheduler_input.lmin[j][i]
                latest_precedent = j

    assignments = []
    for m in i_to_possible_ms[i]:
        processing_time = scheduler_input.p[i][m]
        # find the latest incompatible job that is assigned to the same m
        incompatible_i_end_time = 0
        for assigned_operation in assigned:
            if var_a_im[assigned_operation][m] == 1 and not scheduler_input.C[assigned_operation][i]:
                incompatible_i = assigned_operation
                if var_e_i[incompatible_i] > incompatible_i_end_time:
                    incompatible_i_end_time = var_e_i[incompatible_i]

        for slot in range(scheduler_input.K[m]):
            # find possible start time
            possible_start_times = [
                occupied_until[m][slot],  # cannot start if the slot is occupied
                latest_precedent_end_time,  # cannot start if precedents are not finished
                incompatible_i_end_time,  # cannot start if incompatible are not finished
            ]
            possible_start_time = max(possible_start_times)

            if consider_work_shift:
                # if self.include_shift_constraints:
                # TODO need to validate the work shift list is ordered
                size_n = len(scheduler_input.frak_W)
                if scheduler_input.frak_O[i] in scheduler_input.frak_P:
                    n_assigned = None
                    for n in range(size_n):
                        S_n = scheduler_input.S[n]
                        E_n = scheduler_input.E[n]
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
            assignments.append(
                GreedyAssignment(
                    operation_index=i,
                    start_time=possible_start_time,
                    end_time=possible_end_time,
                    assigned_to_module=m,
                    assigned_to_slot=slot,
                    has_lmax_subsequent=i in i_to_subsequents_lmax,
                    possible_start_times=possible_start_times,
                    latest_precedent=latest_precedent,
                    considers_work_shift=bool(consider_work_shift),
                )
            )
    return min(assignments)


def iterate_by_reaction(operation_network: OperationNetwork, reaction_network: ReactionNetwork,
                        scheduler_input: SchedulerInput) -> list[list[int]]:
    """
    get a list of operation list such that the operation lists are sorted in a way that complies with the
    precedence relations in the reaction network

    :param operation_network:
    :param reaction_network:
    :param scheduler_input:
    :return:
    """

    reaction_to_operation_precedents = defaultdict(list)
    for operation in operation_network.operations:
        reaction_to_operation_precedents[operation.from_reaction] += operation.precedents

    reaction_to_reaction_precedents = defaultdict(list)
    for reaction_1, operation_precedents in reaction_to_operation_precedents.items():
        for reaction_2, operations in operation_network.operations_by_reaction.items():
            if set(operation_precedents).intersection(set([o.identifier for o in operations])):
                if reaction_1 != reaction_2:
                    reaction_to_reaction_precedents[reaction_1].append(reaction_2)
    # rid_ = "6ee813a2-7937-4b6f-b051-564bb59c6132"
    # logger.warning(reaction_to_reaction_precedents[rid_])
    # logger.warning([scheduler_input.frak_O.index(oid) for oid in reaction_to_operation_precedents[rid_]])

    visited = []
    n_reaction = len(reaction_network.chemical_reactions)
    niter = 0
    while len(visited) < n_reaction:
        unvisited = [r.identifier for r in reaction_network.chemical_reactions if r.identifier not in visited]
        reactions_can_visit = [rid for rid in unvisited if
                               set(reaction_to_reaction_precedents[rid]).issubset(set(visited))]
        visited.append(reactions_can_visit[0])
        assert len(visited) <= n_reaction
        niter += 1
        if niter > 1000:
            raise RuntimeError
    list_of_operations = []
    for rid in visited:
        operations_this_reaction = operation_network.operations_by_reaction[rid]
        operations_this_reaction = [scheduler_input.frak_O.index(o.identifier) for o in operations_this_reaction]
        # logger.critical(f"{rid}: {operations_this_reaction}")
        list_of_operations.append(operations_this_reaction)
    return list_of_operations


def register_assign(assign, var_s_i, var_e_i, occupied_until, var_a_im, assigned_global):
    """
    when an assignment is generated, it has to be registered so the greedy algo state is updated

    :param assign:
    :param var_s_i:
    :param var_e_i:
    :param occupied_until:
    :param var_a_im:
    :param assigned_global:
    :return:
    """
    var_s_i[assign.operation_index] = assign.start_time
    var_e_i[assign.operation_index] = assign.end_time
    occupied_until[assign.assigned_to_module][assign.assigned_to_slot] = assign.end_time
    var_a_im[assign.operation_index][assign.assigned_to_module] = 1
    assigned_global.add(assign.operation_index)


def lmax_validate(operations: list[int], var_s_i, var_e_i, scheduler_input: SchedulerInput, eps=0.1):
    """
    validate operations that have a lmax subsequent

    :param operations:
    :param var_s_i:
    :param var_e_i:
    :param scheduler_input:
    :param eps:
    :return:
    """
    assert len(operations)
    valids = []
    for i in operations:
        for j in operations:
            lmax = scheduler_input.lmax[i][j]
            lmin = scheduler_input.lmin[i][j]
            s_j = var_s_i[j]
            e_i = var_e_i[i]
            if lmin > - math.inf and lmax < math.inf:
                valid = bool(e_i + lmax + eps > s_j > e_i + lmin - eps)
                # if not valid:
                #     logger.critical(f"invalid lmax: {i} {j} {s_j - e_i} {lmax}")
                valids.append(valid)
    if len(valids) == 0:
        return True
    return all(valids)


def get_linearized_scheduler_input(scheduler_input: SchedulerInput, operations_grouped_by_reaction) -> SchedulerInput:
    """
    add artificial precedence relations to make a linearized input

    :param scheduler_input:
    :param operations_grouped_by_reaction:
    :return:
    """
    scheduler_input_linearized = copy.deepcopy(scheduler_input)
    for operations in operations_grouped_by_reaction:
        for ii in range(len(operations) - 1):
            i = operations[ii]
            j = operations[ii + 1]
            if scheduler_input.lmin[i][j] == - math.inf:
                scheduler_input_linearized.lmin[i][j] = 0
    return scheduler_input_linearized


class SolverBaseline(Solver):
    """ baseline solver, greedy + linearized precedence relations in scheduler input """

    consider_shifts: bool = True
    """ should work shifts be considered """

    operation_network: OperationNetwork
    """ operation network, required to linearize input """

    reaction_network: ReactionNetwork
    """ reaction network, required to linearize input """

    @property
    def include_shift_constraints(self):
        # TODO DRY
        return self.consider_shifts and self.input.frak_W

    def model_post_init(self, __context: Any) -> None:
        logger.info("\n" + pprint.pformat(self.input.summary))

    def solve(self):
        size_i = len(self.input.frak_O)
        size_m = len(self.input.frak_M)

        # get list of operation lists from a linearized reaction network
        operations_grouped_by_reaction = iterate_by_reaction(self.operation_network, self.reaction_network,
                                                             self.input, )
        order_dict = {
            OperationType.MakeSolution: 0,
            OperationType.TransferLiquid: 1,
            OperationType.Heating: 2,
            OperationType.ConcentrationAndPurification: 3,
        }
        # sort each operation list based on operation type
        for i in range(len(operations_grouped_by_reaction)):
            operations_grouped_by_reaction[i] = sorted(operations_grouped_by_reaction[i], key=lambda x: order_dict[
                self.operation_network.operations[x].type])

        # get scheduler input where operations of the same reaction are linearized
        input_linearized = get_linearized_scheduler_input(scheduler_input=self.input,
                                                          operations_grouped_by_reaction=operations_grouped_by_reaction)

        # record relations among operations
        i_to_precedents = defaultdict(list)  # include both ordinary and specific
        i_to_subsequents = defaultdict(list)
        i_to_precedents_lmax = defaultdict(list)
        i_to_subsequents_lmax = defaultdict(list)
        for i, j in itertools.product(range(size_i), range(size_i)):
            if input_linearized.lmin[i][j] >= 0 or input_linearized.lmax[i][j] < math.inf:
                i_to_precedents[j].append(i)
                i_to_subsequents[i].append(j)
            if input_linearized.lmax[i][j] < math.inf:
                i_to_precedents_lmax[j].append(i)
                i_to_subsequents_lmax[i].append(j)

        # record possible operation-to-module assignment
        i_to_possible_ms = defaultdict(list)
        for i in range(size_i):
            i_to_possible_ms[i] = [m for m in range(size_m) if self.input.p[i][m] < math.inf]

        # record status of slots in modules
        occupied_until = defaultdict(dict)
        for m in range(size_m):
            for slot in range(self.input.K[m]):
                occupied_until[m][slot] = 0.0

        # milp params, make sure baseline output schema is consistent with milp output
        var_a_im = np.zeros((size_i, size_m))
        var_s_i = np.zeros(size_i)
        var_e_i = np.zeros(size_i)

        assigned_global = set()
        for operations_this_reaction in operations_grouped_by_reaction:

            clock_zero = 0  # ad-hoc shift for operation start time
            reaction_assign_valid = False
            # logger.critical(f"already assigned: {assigned_global}")
            # logger.critical(f"trying to assign {operations_this_reaction}")

            while not reaction_assign_valid:
                # use copies for assigning operations of this reaction
                assigned_global_ = copy.deepcopy(assigned_global)
                var_s_i_ = copy.deepcopy(var_s_i)
                var_e_i_ = copy.deepcopy(var_e_i)
                occupied_until_ = copy.deepcopy(occupied_until)
                var_a_im_ = copy.deepcopy(var_a_im)
                # logger.critical(f"now clock {clock_zero}")

                for ii, i in enumerate(operations_this_reaction):
                    assign = get_greedy_assignment(
                        i, occupied_until=occupied_until_, assigned=assigned_global_,
                        var_e_i=var_e_i_,
                        var_a_im=var_a_im_,
                        i_to_precedents=i_to_precedents,
                        i_to_possible_ms=i_to_possible_ms,
                        scheduler_input=input_linearized,
                        consider_work_shift=self.include_shift_constraints,
                        i_to_subsequents_lmax=i_to_subsequents_lmax,
                    )
                    if ii == 0:
                        assign.add_clock_zero(clock_zero, self.input)

                    # if assign.operation_index in [16, 131]:
                    #     logger.warning(
                    #         f"assignment greedy: {assign.operation_index} {self.operation_network.operations[assign.operation_index].type} {assign.start_time} {assign.end_time - assign.start_time}")
                    #     logger.warning(f"assign details: {assign.possible_start_times} {assign.latest_precedent}")
                    #     if assign.operation_index in i_to_subsequents_lmax:
                    #         logger.warning(f"^assigning i that has subsequent lmax!!")
                    #     elif assign.operation_index in i_to_precedents_lmax:
                    #         logger.warning("^assigning i that has precedent lmax!!")
                    #     logger.warning(">>>")

                    register_assign(assign, var_s_i_, var_e_i_, occupied_until_, var_a_im_, assigned_global_)
                reaction_assign_valid = lmax_validate(operations_this_reaction, var_s_i_, var_e_i_, self.input)
                clock_zero += 30
                if clock_zero > 1e4:
                    raise RuntimeError("clock zero too large")

            assigned_global = assigned_global_
            var_s_i = var_s_i_
            var_e_i = var_e_i_
            occupied_until = occupied_until_
            var_a_im = var_a_im_

        self.output = SchedulerOutput.from_MILP(self.input, var_s_i, var_e_i, var_a_im)
        logger.info(f"baseline makespan: {max(self.output.end_times.values())}")
