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


class SolverMILP(Solver):
    infinity: float = 1e5  # use a small infinity to prevent numeric issues

    eps: float = 0.0

    def model_post_init(self, __context: Any) -> None:
        logger.info(pprint.pformat(self.input.summary))

    def solve(self):
        """
        Flexible job shop scheduling with constraints including
        - min and max time lags
        - machine capacity
        - work shifts.

        The formulation is established in eq.2-24 in the main text

        :return:
        """
        # TODO work shift
        # TODO tune up based on gp warnings
        # TODO inspect scale up
        # TODO profile constraints addition & bulk addition implementation

        ts_start = time.time()

        # setup gurobi
        # env = gp.Env(empty=True)
        env = gp.Env()
        env.setParam("OutputFlag", 0)
        env.setParam("LogToConsole", 0)
        env.start()
        model = gp.Model("fjs")

        # get params and cleanup array types
        p = np.array(self.input.p, dtype=float)
        p[p == math.inf] = self.infinity

        lmin = np.array(self.input.lmin, dtype=float)
        lmin[lmin == - math.inf] = - self.infinity

        lmax = np.array(self.input.lmax, dtype=float)
        lmax[lmax == math.inf] = self.infinity

        K = self.input.K

        C = np.array(self.input.C, dtype=int)

        # add vars following table 2
        size_i = len(self.input.frak_O)
        size_m = len(self.input.frak_M)

        # estimate big m
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
        big_m = sum(eq21_lhs) + sum(eq22_lhs) + self.eps
        # big_m = 1e5

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

        ts_added_vars = time.time()

        # add constraints

        # eq. (3)
        # this is allowed in Gurobi 10.0
        # see https://support.gurobi.com/hc/en-us/community/posts/360072196032
        model.addConstr(var_e <= var_e_max, name="eq_3")
        ts_added_eq_3 = time.time()
        self.opt_log['added eq_3'] = ts_added_eq_3 - ts_added_vars

        for i in range(size_i):
            # eq. (4)
            model.addConstr(
                var_e[i] == var_s[i] + gp.quicksum(
                    p[i, m] * var_a[i, m] for m in range(size_m)
                ), name="eq_4"
            )
        ts_added_eq_4 = time.time()
        self.opt_log['added eq_4'] = ts_added_eq_4 - ts_added_eq_3

        for i in range(size_i):
            # eq. (5)
            model.addConstr(gp.quicksum(var_a[i, m] for m in range(size_m)) == 1, name="eq_5")
        ts_added_eq_5 = time.time()
        self.opt_log['added eq_5'] = ts_added_eq_5 - ts_added_eq_4

        # somehow `var_e[i]` is a `MVar`, I have to access the list to get `Var` to be used in `addLConstr`
        for i, j in itertools.product(range(size_i), range(size_i)):
            if i != j:
                # eq. (6)
                model.addLConstr(lhs=var_e_list[i] - var_s_list[j], sense="<=", rhs=-lmin[i, j], name="eq_6")
                # eq. (7)
                model.addLConstr(lhs=var_s_list[j] - var_e_list[i], sense="<=", rhs=lmax[i, j], name="eq_7")

        ts_added_eq_6_eq_7 = time.time()
        self.opt_log['added eq_6 eq_7'] = ts_added_eq_6_eq_7 - ts_added_eq_5

        for i, j in itertools.combinations(range(size_i), 2):  # i < j holds automatically
            for m in range(size_m):
                # eq. (8)
                model.addLConstr(
                    lhs=var_e_list[i], sense="<=",
                    rhs=var_s_list[j] + big_m * (3 - var_x_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq_8"
                )
                # eq. (9)
                model.addLConstr(
                    lhs=var_s_list[j], sense="<=",
                    rhs=var_e_list[i] + big_m * (2 + var_x_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq_9"
                )
                # eq. (10)
                model.addLConstr(
                    lhs=var_e_list[j], sense="<=",
                    rhs=var_s_list[i] + big_m * (3 - var_y_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq_10"
                )
                # eq. (11)
                model.addLConstr(
                    lhs=var_s_list[i], sense="<=",
                    rhs=var_e_list[j] + big_m * (2 + var_y_list[i][j][m] - var_a_list[i][m] - var_a_list[j][m]),
                    name="eq_11"
                )

                # eq. (12)  # TODO maybe faster using `addMConstr`?
                model.addLConstr(
                    lhs=var_x_list[i][j][m] + var_y_list[i][j][m] + var_z_list[i][j][m], sense="=",
                    rhs=1, name="eq_12"
                )

                # # implied, may be good for performance
                # model.addConstr(var_x[i, j, m] + var_y[i, j, m] <= 1, name="implied eq 8-11")

        ts_added_eq_8_9_10_11_12 = time.time()
        self.opt_log['added eq_8 eq_9 eq_10 eq_11 eq_12'] = ts_added_eq_8_9_10_11_12 - ts_added_eq_6_eq_7

        # eq. (13)
        for i, m in itertools.product(range(size_i), range(size_m)):
            model.addLConstr(
                lhs=gp.quicksum(var_z_list[i][j][m] for j in range(size_i) if i != j), sense="<=",
                rhs=(K[m] - 1) * var_a_list[i][m] + self.eps, name="eq_13"
            )
        ts_added_eq_13 = time.time()
        self.opt_log['added eq_13'] = ts_added_eq_13 - ts_added_eq_8_9_10_11_12

        # eq. (14) compatability
        for i, j, m in itertools.product(range(size_i), range(size_i), range(size_m)):
            model.addLConstr(lhs=var_z_list[i][j][m], sense="<=", rhs=C[i, j], name="eq_14")  # C[i, i] is 1 by default
        ts_added_eq_14 = time.time()
        self.opt_log['added eq_14'] = ts_added_eq_14 - ts_added_eq_14

        ts_added_constraints = time.time()
        logger.warning("finish adding constraints")

        model.setObjective(var_e_max, GRB.MINIMIZE)
        model.optimize()

        ts_solved = time.time()

        if model.Status == GRB.OPTIMAL:
            logger.info("optimal solution found!")
            logger.info(f"the solution is: {model.objVal}")
        else:
            logger.warning(model.Status)

        env.close()

        var_s_values = np.empty(size_i)
        var_e_values = np.empty(size_i)
        var_a_values = np.empty((size_i, size_m))
        for i in range(size_i):
            var_s_values[i] = var_s[i].X
            var_e_values[i] = var_e[i].X
            for m in range(size_m):
                var_a_values[i, m] = var_a[i, m].X

        # # print out all vars
        # all_vars = model.getVars()
        # values = model.getAttr("X", all_vars)
        # names = model.getAttr("VarName", all_vars)
        # for name, val in zip(names, values):
        #     print(f"{name} = {val}")

        self.output = SchedulerOutput.from_MILP(self.input, var_s_values, var_e_values, var_a_values)
        self.opt_log["time adding vars"] = ts_added_vars - ts_start
        self.opt_log["time adding constraints"] = ts_added_constraints - ts_added_vars
        self.opt_log["time solved"] = ts_solved - ts_added_constraints
        self.opt_log["big m estimated as"] = big_m

        logger.info("\n" + pprint.pformat(self.opt_log))

        # other info that can be added to log
        # https://www.gurobi.com/documentation/current/refman/attributes.html
