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
    # infinity: float = GRB.INFINITY
    infinity: float = 1e10

    eps: float = 1e-6

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
        p = np.array(self.input.p)
        p[p == math.inf] = self.infinity

        lmin = np.array(self.input.lmin)
        lmin[lmin == - math.inf] = - self.infinity

        lmax = np.array(self.input.lmax)
        lmax[lmax == math.inf] = self.infinity

        K = self.input.K

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

        ts_added_vars = time.time()

        # add constraints
        for i in range(size_i):
            # eq. (3)
            model.addConstr(var_e_max >= var_e[i], name="eq_3")
            # eq. (4)
            model.addConstr(
                var_e[i] == var_s[i] + gp.quicksum(
                    p[i, m] * var_a[i, m] for m in range(size_m)
                ), name="eq_4"
            )
            # eq. (5)
            model.addConstr(gp.quicksum(var_a[i, m] for m in range(size_m)) == 1, name="eq_5")

        # TODO maybe `combinations` is enough?
        for i, j in itertools.product(range(size_i), range(size_i)):
            if i != j:
                # eq. (6)
                model.addConstr(var_s[j] >= var_e[i] + lmin[i, j], name="eq_6")
                # eq. (7)
                model.addConstr(var_s[j] <= var_e[i] + lmax[i, j], name="eq_7")

        for i, j in itertools.combinations(range(size_i), 2):  # i < j holds automatically
            for m in range(size_m):
                # eq. (8)
                model.addConstr(
                    var_e[i] <= var_s[j] + big_m * (3 - var_x[i, j, m] - var_a[i, m] - var_a[j, m]), name="eq_8"
                )
                # eq. (9)
                model.addConstr(
                    var_e[i] >= var_s[j] - big_m * (2 + var_x[i, j, m] - var_a[i, m] - var_a[j, m]), name="eq_9"
                )
                # eq. (10)
                model.addConstr(
                    var_e[j] <= var_s[i] + big_m * (3 - var_y[i, j, m] - var_a[i, m] - var_a[j, m]), name="eq_10"
                )
                # eq. (11)
                model.addConstr(
                    var_e[j] >= var_s[i] - big_m * (2 + var_y[i, j, m] - var_a[i, m] - var_a[j, m]), name="eq_11"
                )

                # implied, may be good for performance
                model.addConstr(var_x[i, j, m] + var_y[i, j, m] <= 1, name="implied eq 8-11")

                # eq. (12)
                model.addConstr(var_x[i, j, m] + var_y[i, j, m] + var_z[i, j, m] == 1, name="eq_12")

        # eq. (13)
        for i, m in itertools.product(range(size_i), range(size_m)):
            model.addConstr(
                gp.quicksum(var_z[i, j, m] for j in range(size_i) if i != j) <= (K[m] - 1) * var_a[i, m] + self.eps,
                name="eq_13",
            )

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

        logger.info(pprint.pformat(self.output.operation_view()))
        logger.info("\n" + pprint.pformat(self.opt_log))

        # other info that can be added to log
        # https://www.gurobi.com/documentation/current/refman/attributes.html
