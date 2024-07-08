import itertools
import math
import pprint
import time
from typing import Any

import gurobipy as gp
import numpy as np
from gurobipy import GRB
from loguru import logger

from .schema import Solver, SchedulerOutput
from ..utils import FilePath


class SolverMILP(Solver):
    infinity: float = 1e6  # use a small infinity to prevent numeric issues

    eps: float = 0.0

    dummy_big_m: bool = False

    consider_shifts: bool = True

    supress_gp_log: bool = False

    @property
    def include_shift_constraints(self):
        return self.consider_shifts and self.input.frak_W

    def model_post_init(self, __context: Any) -> None:
        logger.info("\n" + pprint.pformat(self.input.summary))

    def prepare_input(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int], np.ndarray, int, int]:
        # get params and cleanup array types
        p = np.array(self.input.p, dtype=float)
        p[p == math.inf] = self.infinity

        lmin = np.array(self.input.lmin, dtype=float)
        lmin[lmin == - math.inf] = - self.infinity

        lmax = np.array(self.input.lmax, dtype=float)
        lmax[lmax == math.inf] = self.infinity

        K = self.input.K

        C = np.array(self.input.C, dtype=int)

        size_i = len(self.input.frak_O)
        size_m = len(self.input.frak_M)
        return p, lmin, lmax, K, C, size_i, size_m

    def model_add_vars(self, model: gp.Model, size_i, size_m):
        # add vars following table 2
        var_e_max = model.addVar(name="var_e_max", vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
        var_s = model.addMVar(size_i, vtype=GRB.CONTINUOUS, name="var_s_i", lb=0.0, ub=GRB.INFINITY)
        var_e = model.addMVar(size_i, vtype=GRB.CONTINUOUS, name="var_e_i", lb=0.0, ub=GRB.INFINITY)
        var_a = model.addMVar((size_i, size_m), vtype=GRB.BINARY, name="var_a_im")
        var_x = model.addMVar((size_i, size_i, size_m), vtype=GRB.BINARY, name="var_x_ijm")
        var_y = model.addMVar((size_i, size_i, size_m), vtype=GRB.BINARY, name="var_y_ijm")
        var_z = model.addMVar((size_i, size_i, size_m), vtype=GRB.BINARY, name="var_z_ijm")
        var_e_list = var_e.tolist()
        var_s_list = var_s.tolist()
        var_a_list = var_a.tolist()
        var_x_list = var_x.tolist()
        var_y_list = var_y.tolist()
        var_z_list = var_z.tolist()
        return var_e_max, var_s, var_e, var_a, var_x, var_y, var_z, var_e_list, var_s_list, var_a_list, var_x_list, var_y_list, var_z_list

    def model_add_main_constraints(
            self, model: gp.Model, big_m, p, lmin, lmax, K, C, size_i, size_m,
            var_e_max,
            var_s, var_e, var_a,
            var_e_list, var_s_list, var_a_list, var_x_list, var_y_list, var_z_list
    ):
        """
        main constraints for flexible job shop scheduling with constraints including
        - min and max time lags
        - machine capacity and job compatability
        """

        ts_add_constraints_start = time.time()

        # eq:emax
        # this is allowed in Gurobi 10.0
        # see https://support.gurobi.com/hc/en-us/community/posts/360072196032
        model.addConstr(var_e <= var_e_max, name="eq:emax")
        ts_added_eq_emax = time.time()
        self.opt_log['added eq:emax'] = ts_added_eq_emax - ts_add_constraints_start

        # eq:processtime
        for i in range(size_i):
            model.addConstr(
                var_e[i] == var_s[i] + gp.quicksum(
                    p[i, m] * var_a[i, m] for m in range(size_m)
                ), name="eq:processtime"
            )
        ts_added_eq_processtime = time.time()
        self.opt_log['eq:processtime'] = ts_added_eq_processtime - ts_added_eq_emax

        # eq:assign
        for i in range(size_i):
            model.addConstr(gp.quicksum(var_a[i, m] for m in range(size_m)) == 1, name="eq:assign")
        ts_added_eq_assign = time.time()
        self.opt_log['added eq:assign'] = ts_added_eq_assign - ts_added_eq_processtime

        # somehow `var_e[i]` is a `MVar`, I have to access the list to get `Var` to be used in `addLConstr`
        for i, j in itertools.product(range(size_i), range(size_i)):
            if i != j:
                # eq:precedence1
                model.addLConstr(lhs=var_e_list[i] - var_s_list[j], sense="<=", rhs=-lmin[i, j], name="eq:precedence1")
                # eq:precedence2
                model.addLConstr(lhs=var_s_list[j] - var_e_list[i], sense="<=", rhs=lmax[i, j], name="eq:precedence2")
        ts_added_eq_precedence12 = time.time()
        self.opt_log['added eq:precedence1 eq:precedence2'] = ts_added_eq_precedence12 - ts_added_eq_assign

        for i, j in itertools.combinations(range(size_i), 2):  # i < j holds automatically
            for m in range(size_m):
                # eq:overlap1
                model.addLConstr(
                    lhs=var_e_list[i], sense="<=",
                    rhs=var_s_list[j] + big_m * (3 - var_x_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq:overlap1"
                )
                # eq:overlap2
                model.addLConstr(
                    lhs=var_s_list[j], sense="<=",
                    rhs=var_e_list[i] + big_m * (2 + var_x_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq:overlap2"
                )
                # eq:overlap3
                model.addLConstr(
                    lhs=var_e_list[j], sense="<=",
                    rhs=var_s_list[i] + big_m * (3 - var_y_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq:overlap3"
                )
                # eq:overlap4
                model.addLConstr(
                    lhs=var_s_list[i], sense="<=",
                    rhs=var_e_list[j] + big_m * (2 + var_y_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq:overlap4"
                )

                # # implied, may be good for performance
                # model.addConstr(var_x[i, j, m] + var_y[i, j, m] <= 1, name="implied eq 8-11")

                # eq:overlap5
                model.addLConstr(
                    lhs=var_x_list[i][j][m] + var_y_list[i][j][m] + var_z_list[i][j][m], sense="=",
                    rhs=1, name="eq:overlap5"
                )

        ts_added_eq_overlap12345 = time.time()
        self.opt_log['added eq:overlap12345'] = ts_added_eq_overlap12345 - ts_added_eq_precedence12

        # eq:capacity
        for i, m in itertools.product(range(size_i), range(size_m)):
            model.addLConstr(
                lhs=gp.quicksum(var_z_list[i][j][m] for j in range(size_i) if i != j), sense="<=",
                rhs=(K[m] - 1) * var_a_list[i][m] + self.eps, name="eq:capacity"
            )
        ts_added_eq_capacity = time.time()
        self.opt_log['added eq:capacity'] = ts_added_eq_capacity - ts_added_eq_overlap12345

        # eq:compatibility
        for i, j, m in itertools.product(range(size_i), range(size_i), range(size_m)):
            model.addLConstr(lhs=var_z_list[i][j][m], sense="<=", rhs=C[i, j],
                             name="eq:compatibility")  # C[i, i] is 1 by default
        ts_added_eq_compatibility = time.time()
        self.opt_log['added eq:compatibility'] = ts_added_eq_compatibility - ts_added_eq_compatibility

        ts_add_main_constraints_end = time.time()
        self.opt_log['total adding main constraints'] = ts_add_main_constraints_end - ts_add_constraints_start

    def model_add_work_shift_constraints(
            self, model: gp.Model, big_m, size_i,
            var_e_list, var_s_list,
    ):
        # TODO separate param, var, and constraint addition
        # work shift params
        size_n = len(self.input.frak_W)
        S = self.input.S
        E = self.input.E
        frak_P = [self.input.frak_O.index(oid) for oid in self.input.frak_P]

        var_Y = model.addMVar((size_i, size_n), vtype=GRB.BINARY, name="var_Y_in")
        var_Z = model.addMVar((size_i, size_n), vtype=GRB.BINARY, name="var_Z_in")
        var_A = model.addMVar((size_i, size_n), vtype=GRB.BINARY, name="var_A_in")
        var_Y_list = var_Y.tolist()
        var_Z_list = var_Z.tolist()
        var_A_list = var_A.tolist()

        # eq:shift1
        for i in frak_P:
            model.addConstr(gp.quicksum(var_A[i, n] for n in range(size_n)) == 1, name="eq:shift1")

        for i, n in itertools.product(frak_P, range(size_n)):
            # eq:shift2
            model.addLConstr(
                lhs=S[n],
                sense="<=",
                rhs=var_s_list[i] + big_m * (1 - var_Y_list[i][n]),
                name="eq:shift2"
            )
            # eq:shift3
            model.addLConstr(
                lhs=var_s_list[i],
                sense="<=",
                rhs=S[n] + big_m * var_Y_list[i][n],
                name="eq:shift3"
            )
            # eq:shift4
            model.addLConstr(
                lhs=var_e_list[i],
                sense="<=",
                rhs=E[n] + big_m * (1 - var_Z_list[i][n]),
                name="eq:shift4"
            )
            # eq:shift5
            model.addLConstr(
                lhs=E[n],
                sense="<=",
                rhs=var_e_list[i] + big_m * var_Z_list[i][n],
                name="eq:shift5"
            )
            # eq (19)
            model.addLConstr(
                lhs=var_Y_list[i][n] + var_Z_list[i][n] - 2 * var_A_list[i][n],
                sense=">=",
                rhs=0,
                name="eq:shift6"
            )
            # eq (20)
            model.addLConstr(
                lhs=var_A_list[i][n] - var_Y_list[i][n] - var_Z_list[i][n],
                sense=">=",
                rhs=-1,
                name="eq:shift7"
            )

    def estimate_big_m(self, p, lmin, size_i, size_m, dummy: bool = False) -> float:
        if dummy:
            return self.infinity
        eq21_lhs = []  # i.e. worst assignment
        for i in range(size_i):
            p_i_max = 0
            for m in range(size_m):
                pt = p[i][m]
                if pt < self.infinity and pt > p_i_max:
                    p_i_max = pt
            eq21_lhs.append(p_i_max)
        eq22_lhs = []
        for i in range(size_i):
            lmin_i_max = 0
            for j in range(size_i):
                _lmin_i_j = lmin[i][j]
                if _lmin_i_j > lmin_i_max:
                    lmin_i_max = _lmin_i_j
            eq22_lhs.append(lmin_i_max)
        big_m = sum(eq21_lhs) + sum(eq22_lhs)

        # eq (23)
        if self.include_shift_constraints:
            E_max = 0
            for w in self.input.frak_W:
                if w.end_time > E_max:
                    E_max = w.end_time
            if E_max > big_m:
                big_m = E_max

        if big_m > self.infinity * 0.1:
            raise RuntimeError(f"big_m is estimated ({big_m}) to be rather close to self.infinity ({self.infinity})!")

        return big_m + self.eps

    def solve(self, logfile:FilePath = None):
        # TODO inspect scale up

        ts_start = time.time()

        # prepare input
        p, lmin, lmax, K, C, size_i, size_m = self.prepare_input()

        # estimate big m, if dummy then self.infinity is used
        big_m = self.estimate_big_m(p, lmin, size_i, size_m, self.dummy_big_m)

        # setup gurobi env
        env = gp.Env(empty=self.supress_gp_log)
        if self.supress_gp_log:
            env.setParam("OutputFlag", 0)
            env.setParam("LogToConsole", 0)
        env.start()

        # model init
        model = gp.Model("FJSS-MILP")

        if logfile:
            model.setParam("LogFile", logfile)

        if self.time_limit:
            model.setParam("TimeLimit", self.time_limit)

        # focus on feasible solution
        model.setParam("MIPFocus", 1)

        # based on log, see https://support.gurobi.com/hc/en-us/community/posts/360077323472
        model.setParam("BarHomogeneous", 1)

        # add variables
        (
            var_e_max,
            var_s, var_e, var_a, var_x, var_y, var_z,
            var_e_list, var_s_list, var_a_list, var_x_list, var_y_list, var_z_list
        ) = self.model_add_vars(model, size_i, size_m)
        ts_added_vars = time.time()

        # add main constraints
        self.model_add_main_constraints(
            model, big_m, p, lmin, lmax, K, C, size_i, size_m,
            var_e_max,
            var_s, var_e, var_a,
            var_e_list, var_s_list, var_a_list, var_x_list, var_y_list, var_z_list
        )
        ts_added_main_constraints = time.time()

        # add shift constraints
        if self.include_shift_constraints:
            self.model_add_work_shift_constraints(model, big_m, size_i, var_e_list, var_s_list)
        ts_added_shift_constraints = time.time()

        # actual solve
        model.setObjective(var_e_max, GRB.MINIMIZE)
        model.optimize()
        ts_solved = time.time()

        if model.Status == GRB.OPTIMAL:
            logger.info("optimal solution found!")
            logger.info(f"the solution is: {model.objVal}")
        else:
            logger.warning(model.Status)
        env.close()

        # logging
        self.opt_log["time adding vars"] = ts_added_vars - ts_start
        self.opt_log["time adding main constraints"] = ts_added_main_constraints - ts_added_vars
        self.opt_log["time adding shift constraints"] = ts_added_shift_constraints - ts_added_main_constraints
        self.opt_log["time solved"] = ts_solved - ts_added_shift_constraints
        self.opt_log["big m estimated as"] = big_m
        self.opt_log["gurobi status"] = model.Status
        self.opt_log["gurobi solution count"] = model.SolCount
        logger.info("\n" + pprint.pformat(self.opt_log))

        if model.SolCount:
            # collect output
            var_s_values = np.empty(size_i)
            var_e_values = np.empty(size_i)
            var_a_values = np.empty((size_i, size_m))
            for i in range(size_i):
                var_s_values[i] = var_s[i].X
                var_e_values[i] = var_e[i].X
                for m in range(size_m):
                    var_a_values[i, m] = var_a[i, m].X
            self.output = SchedulerOutput.from_MILP(self.input, var_s_values, var_e_values, var_a_values)
        else:
            self.output = SchedulerOutput()  # empty output when no solution found

        # # print out all vars
        # all_vars = model.getVars()
        # values = model.getAttr("X", all_vars)
        # names = model.getAttr("VarName", all_vars)
        # for name, val in zip(names, values):
        #     print(f"{name} = {val}")

        # other info that can be added to log
        # https://www.gurobi.com/documentation/current/refman/attributes.html
