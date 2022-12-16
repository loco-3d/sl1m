import numpy as np
import quadprog
from enum import Enum
from scipy.optimize import linprog
from time import perf_counter as clock

try:
    import glpk
except ImportError:
    print("Import error: No module GLPK")
    pass

try:
    import gurobipy as grb

    grb.setParam("LogFile", "")
    grb.setParam("OutputFlag", 0)
except ImportError:
    print("Import error: No module Gurobipy")
    pass

try:
    import cvxpy
except ImportError:
    print("Import error: No module cvxpy")
    pass


# -------------------------- RESULT DATA STRUCTURE -------------------------------------------------


class ResultData:
    def __init__(self, success, time, x=None):
        self.success = success
        self.time = time
        self.x = x

    def __str__(self):
        string = "ResultData: "
        string += "\n    Success: " + str(self.success)
        string += "\n    Time: " + str(self.time)
        return string

    def __repr__(self):
        return self.__str__()


# -------------------------- CONSTANTS AND HELPERS -------------------------------------------------

EPSILON = 0.00001


class Solvers(Enum):
    GUROBI = 0
    GLPK = 1
    LINPROG = 2
    QUADPROG = 3
    CVXPY = 4


def call_LP_solver(q, G=None, h=None, C=None, d=None, solver=Solvers.GUROBI):
    """
    Solve the LP problem with a specific solver
    :param: q, G, h, C, d problem data
    :param: SOLVER Solver choice
    :return: None if wrong SOLVER, else ResultData
    """
    if solver == Solvers.GUROBI:
        result = solve_lp_gurobi(q, G, h, C, d)
    elif solver == Solvers.GLPK:
        result = solve_lp_glpk(q, G, h, C, d)
    elif solver == Solvers.LINPROG:
        result = solve_lp_linprog(q, G, h, C, d)
    elif solver == Solvers.CVXPY:
        result = solve_qp_lp_cvxpy(None, None, G, h, C, d)
    else:
        print("Unknown LP solver asked : ", solver)
        return
    return result


def call_QP_solver(P, q, G=None, h=None, C=None, d=None, solver=Solvers.GUROBI):
    """
    Solve the QP problem with a specific solver
    :param: P, q, G, h, C, d problem data
    :param: SOLVER Solver choice
    :return: None if wrong SOLVER, else ResultData
    """
    if solver == Solvers.GUROBI:
        result = solve_qp_gurobi(P, q, G, h, C, d)

    elif solver == Solvers.QUADPROG:
        result = solve_qp_quaprog(P, q, G, h, C, d)
    elif solver == Solvers.CVXPY:
        result = solve_qp_lp_cvxpy(P, q, G, h, C, d)
    else:
        print("Unknown QP solver asked : ", solver)
        return
    return result


def call_MIP_solver(
    slack_selection_vector,
    P=None,
    q=None,
    G=None,
    h=None,
    C=None,
    d=None,
    solver=Solvers.GUROBI,
):
    """
    Solve the MIP problem with a specific solver
    :param: q, G, h, C, d problem data
    :param: SOLVER Solver choice
    :param: slack_selection_vector Slack variables selection matrice of the problem
    :return: None if wrong SOLVER, else ResultData
    """
    hasCost = not (P is None or q is None)
    result = None
    if hasCost:
        if solver == Solvers.GUROBI:
            result = solve_MIP_gurobi_cost(slack_selection_vector, P, q, G, h, C, d)
        elif solver == Solvers.CVXPY:
            result = solve_MIP_cvxpy(slack_selection_vector, G, h, C, d)
        else:
            print("Unknown MIP solver asked : ", solver)
            return
    else:
        if solver == Solvers.GUROBI:
            result = solve_MIP_gurobi(slack_selection_vector, G, h, C, d)

        elif solver == Solvers.CVXPY:
            result = solve_MIP_cvxpy(slack_selection_vector, G, h, C, d)
        else:
            print("Unknown MIP solver asked : ", solver)
            return
    return result


def ms(t):
    """
    Convert a time in microseconds to milliseconds
    """
    return t * 1000.0


def get_nonzero_rows(M):
    """
    Return the Matrix M with only the rows not null
    """
    M_result = {}
    rows, cols = M.nonzero()
    for ij in zip(rows, cols):
        i, j = ij
        if i not in M_result:
            M_result[i] = []
        M_result[i].append(j)
    return M_result


# -------------------------- SOLVE METHODS ---------------------------------------------------------


def solve_least_square(A, b, G=None, h=None, C=None, d=None, solver=Solvers.GUROBI):
    """
    min ||Ax-b||**2
    subject to  G x <= h
    subject to  C x  = d
    """
    P = np.dot(A.T, A)
    q = -np.dot(A.T, b).reshape(A.shape[1])
    return call_QP_solver(P, q, G, h, C, d, solver)


# -------------------------- LINEAR PROGRAM --------------------------------------------------------


def solve_lp_linprog(q, G=None, h=None, C=None, d=None):
    """
    Solve linear programm using scipy.linprog
    min q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    bounds = [(-100000.0, 10000.0) for _ in range(q.shape[0])]
    res = linprog(
        q,
        A_ub=G,
        b_ub=h,
        A_eq=C,
        b_eq=d,
        bounds=bounds,
        method="interior-point",
        callback=None,
        options={"presolve": True},
    )
    t_end = clock()
    return ResultData(res["success"], ms(t_end - t_init), res["x"])


def solve_lp_glpk(q, G=None, h=None, C=None, d=None):
    """
    Solve linear programm using the simplex method with glpk
    min q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    lp = glpk.LPX()
    lp.obj.maximize = False

    n_equality_constraints = 0
    n_inequality_constraints = 0
    if C is not None:
        n_equality_constraints = C.shape[0]
        n_variables = C.shape[1]
    if G is not None:
        n_inequality_constraints = G.shape[0]
        n_variables = G.shape[1]
    n_constraints = n_equality_constraints + n_inequality_constraints
    lp.cols.add(n_variables)

    for c in lp.cols:
        c.bounds = None, None
    if n_constraints > 0:
        lp.rows.add(n_constraints)
        i_start = 0
        mat = []
        for i in range(n_inequality_constraints):
            lp.rows[i_start].bounds = None, h[i]
            for j in range(n_variables):
                if abs(G[i, j]) > EPSILON:
                    mat.append((i_start, j, G[i, j]))
            i_start += 1
        for i in range(n_equality_constraints):
            lp.rows[i_start].bounds = d[i]
            for j in range(n_variables):
                if abs(C[i, j]) > EPSILON:
                    mat.append((i_start, j, C[i, j]))
            i_start += 1
        lp.matrix = mat

    lp.obj[:] = q.tolist()
    lp.simplex()
    t_end = clock()
    return ResultData(
        lp.status == "opt",
        ms(t_end - t_init),
        np.array([c.primal for c in lp.cols]),
    )


def solve_lp_gurobi(q, G=None, h=None, C=None, d=None):
    """
    Solve linear programm using the simplex method with glpk
    min q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    model = grb.Model("lp")

    # add continuous variables
    cVars = []
    for el in q:
        cVars.append(
            model.addVar(
                vtype=grb.GRB.CONTINUOUS,
                lb=-grb.GRB.INFINITY,
                ub=grb.GRB.INFINITY,
            )
        )

    # Update model to integrate new variables
    model.update()
    x = np.array(model.getVars(), copy=False)

    # Equality constraints
    if C.shape[0] > 0:
        for i in range(C.shape[0]):
            idx = [j for j, el in enumerate(C[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = C[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, d[i])

    # Inequality constraints
    if G.shape[0] > 0:
        for i in range(G.shape[0]):
            idx = [j for j, el in enumerate(G[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = G[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, h[i])

    slack_indices = [i for i, el in enumerate(q) if el > 0]

    # set objective
    variables = []
    previousL = 0
    obj = 0
    for i, el in enumerate(slack_indices):
        if i != 0 and el - previousL > 2.0:
            assert len(variables) > 0
            expr = grb.LinExpr(np.ones(len(variables)), variables)
            obj += expr
            variables = [x[el]]
        elif el != 0:
            variables += [x[el]]
        previousL = el
    if len(variables) > 1:
        expr = grb.LinExpr(np.ones(len(variables)), variables)
        obj += expr

    model.setObjective(obj, grb.GRB.MINIMIZE)
    model.optimize()
    t_end = clock()
    try:
        res = [el.x for el in cVars]
        return ResultData(model.Status == grb.GRB.OPTIMAL, ms(t_end - t_init), res)
    except:
        return ResultData(False, ms(t_end - t_init))


# -------------------------- QUADRATIC PROGRAM -----------------------------------------------------


def solve_qp_quaprog(P, q, G=None, h=None, C=None, d=None):
    """
    Solve the QP problem using Quadprog
    min (1/2)x' P x + q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    qp_G = 0.5 * (P + P.T)
    qp_a = -q
    if C is not None and d is not None:
        if G is not None and h is not None:
            qp_C = -np.vstack([C, G]).T
            qp_b = -np.hstack([d, h])
        else:
            qp_C = -C.transpose()
            qp_b = -d
        meq = C.shape[0]
    else:
        qp_C = -G.T
        qp_b = -h
        meq = 0
    try:
        res = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)
        t_end = clock()
        return ResultData(True, ms(t_end - t_init), res[0])
    except:
        t_end = clock()
        return ResultData(False, ms(t_end - t_init))


def solve_qp_lp_cvxpy(P, q, G=None, h=None, C=None, d=None):
    """
    Solve the Mixed-Integer problem using CVXPY
    find x
    min x'P'x + q' x
    subject to  G x <= h
    subject to  C x  = d
    """

    t_init = clock()
    if G is not None:
        n_variables = G.shape[1]
    elif C is not none:
        n_variables = C.shape[1]

    variables = cvxpy.Variable(n_variables)
    obj = 0.0
    if P is not None:
        n_variables = P.shape[1]
        obj = cvxpy.Minimize(cvxpy.quad_form(variables, P) + q.T * variables)

    constraints = []
    if G.shape[0] > 0:
        ineq_constraints = G * variables <= h
        constraints.append(ineq_constraints)

    if C.shape[0] > 0:
        eq_constraints = C * variables == d
        constraints.append(eq_constraints)

    prob = cvxpy.Problem(obj, constraints)
    prob.solve(verbose=False)
    t_end = clock()
    if prob.status not in ["infeasible", "unbounded"]:
        res = np.array([v.value for v in variables])
        return ResultData(True, ms(t_end - t_init), res)
    else:
        return ResultData(False, ms(t_end - t_init))


def solve_qp_gurobi(P, q, G=None, h=None, C=None, d=None, verbose=False):
    """
    Solve the QP problem using Quadprog
    min (1/2)x' P x + q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    grb.setParam("OutputFlag", 1 if verbose else 0)
    n = P.shape[1]
    model = grb.Model("qp")

    cVars = []
    for i in range(n):
        cVars.append(
            model.addVar(
                vtype=grb.GRB.CONTINUOUS,
                lb=-grb.GRB.INFINITY,
                ub=grb.GRB.INFINITY,
            )
        )

    # Update model to integrate new variables
    model.update()
    x = np.array(model.getVars(), copy=False)
    obj = grb.QuadExpr()
    rows, cols = P.nonzero()
    for i, j in zip(rows, cols):
        obj += 0.5 * x[i] * P[i, j] * x[j]
    for i in range(n):
        obj += q[i] * x[i]
    model.setObjective(obj, grb.GRB.MINIMIZE)

    # Inquality constraints
    if G is not None:
        G_nonzero_rows = get_nonzero_rows(G)
        for i, row in G_nonzero_rows.items():
            model.addConstr(grb.quicksum(G[i, j] * x[j] for j in row) <= h[i])

    # Equality constraints
    if C is not None:
        C_nonzero_rows = get_nonzero_rows(C)
        for i, row in C_nonzero_rows.items():
            model.addConstr(grb.quicksum(C[i, j] * x[j] for j in row) == d[i])

    model.optimize()
    t_end = clock()
    try:
        res = [el.x for el in cVars]
        return ResultData(model.Status == grb.GRB.OPTIMAL, ms(t_end - t_init), res)
    except:
        return ResultData(False, ms(t_end - t_init))


# -------------------------- MIXED_INTEGER PROGRAM -------------------------------------------------


def solve_MIP_gurobi(slack_selection_vector, G=None, h=None, C=None, d=None):
    """
    Solve the Mixed-Integer problem using Gurobipy
    Choose one alpha per phase
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    model = grb.Model("mip")

    cVars = []
    for variable in slack_selection_vector:
        if variable > 0:
            cVars.append(model.addVar(lb=0, ub=1, vtype=grb.GRB.BINARY, name="slack"))
        else:
            cVars.append(
                model.addVar(
                    lb=-grb.GRB.INFINITY,
                    ub=grb.GRB.INFINITY,
                    vtype=grb.GRB.CONTINUOUS,
                    name="x",
                )
            )

    # Update model to integrate new variables
    model.update()
    x = np.array(model.getVars(), copy=False)

    # Inequality constraints
    if G.shape[0] > 0:
        for i in range(G.shape[0]):
            idx = [j for j, el in enumerate(G[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = G[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, h[i])
    model.update()

    # Equality constraints
    if C.shape[0] > 0:
        for i in range(C.shape[0]):
            idx = [j for j, el in enumerate(C[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = C[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, d[i])
    model.update()

    slack_indices = [i for i, el in enumerate(slack_selection_vector) if el > 0]

    # equality slack sum
    variables = []
    previousL = 0
    for i, el in enumerate(slack_indices):
        if i != 0 and el - previousL > 2.0:
            assert len(variables) > 0
            expr = grb.LinExpr(np.ones(len(variables)), variables)
            model.addConstr(expr, grb.GRB.EQUAL, len(variables) - 1)
            variables = [x[el]]
        elif el != 0:
            variables += [x[el]]
        previousL = el
    if len(variables) > 1:
        expr = grb.LinExpr(np.ones(len(variables)), variables)
        model.addConstr(expr, grb.GRB.EQUAL, len(variables) - 1)
    model.update()
    model.optimize()
    t_end = clock()
    try:
        res = [el.x for el in cVars]
        return ResultData(model.Status == grb.GRB.OPTIMAL, ms(t_end - t_init), res)
    except:
        return ResultData(False, ms(t_end - t_init))


def solve_MIP_gurobi_cost(slack_selection_vector, P, q, G=None, h=None, C=None, d=None):
    """
    Solve the Mixed-Integer problem using Gurobipy
    min (1/2)x' P x + q' x
    subject to  G x <= h
    subject to  C x  = d
    """
    grb.setParam("LogFile", "")
    grb.setParam("OutputFlag", 0)

    slack_indices = [i for i, el in enumerate(slack_selection_vector) if el > 0]
    n_variables = len(slack_selection_vector)

    model = grb.Model("mip")

    # add continuous variables
    cVars = []
    for variable in slack_selection_vector:
        if variable > 0:
            cVars.append(model.addVar(lb=0, ub=1, vtype=grb.GRB.BINARY, name="slack"))
        else:
            cVars.append(
                model.addVar(
                    lb=-grb.GRB.INFINITY,
                    ub=grb.GRB.INFINITY,
                    vtype=grb.GRB.CONTINUOUS,
                    name="x",
                )
            )

    # Update model to integrate new variables
    model.update()
    x = np.array(model.getVars(), copy=False)

    # Inequality constraints
    if G.shape[0] > 0:
        for i in range(G.shape[0]):
            idx = [j for j, el in enumerate(G[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = G[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, h[i])
    model.update()

    # Equality constraints
    if C.shape[0] > 0:
        for i in range(C.shape[0]):
            idx = [j for j, el in enumerate(C[i].tolist()) if el != 0.0]
            variables = x[idx]
            coeff = C[i, idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, d[i])
    model.update()

    obj = grb.QuadExpr()
    rows, cols = P.nonzero()
    for i, j in zip(rows, cols):
        obj += 0.5 * x[i] * P[i, j] * x[j]
    for i in range(n_variables):
        obj += q[i] * x[i]
    model.setObjective(obj, grb.GRB.MINIMIZE)

    t_init = clock()
    model.optimize()
    t_end = clock()
    try:
        res = [variable.x for variable in cVars]
        return ResultData(model.Status == grb.GRB.OPTIMAL, ms(t_end - t_init), res)
    except:
        print("Failed to solve the MIP")
        return ResultData(False, ms(t_end - t_init))


def solve_MIP_cvxpy(slack_selection_vector, G=None, h=None, C=None, d=None):
    """
    Solve the Mixed-Integer problem using CVXPY
    find x
    subject to  G x <= h
    subject to  C x  = d
    """
    t_init = clock()
    n_variables = G.shape[1]
    variables = cvxpy.Variable(n_variables)
    constraints = []
    if G.shape[0] > 0:
        ineq_constraints = G * variables <= h
        constraints.append(ineq_constraints)

    if C.shape[0] > 0:
        eq_constraints = C * variables == d
        constraints.append(eq_constraints)

    slack_indices = [i for i, el in enumerate(slack_selection_vector) if el > 0]
    n_slack_variables = len([el for el in slack_selection_vector if el > 0])
    obj = cvxpy.Minimize(slack_selection_vector * variables)

    if n_slack_variables > 0:
        boolvars = cvxpy.Variable(n_slack_variables, boolean=True)
        constraints = constraints + [
            variables[el] <= boolvars[i] for i, el in enumerate(slack_indices)
        ]

    currentSum = []
    previousL = 0
    for i, el in enumerate(slack_indices):
        if i != 0 and el - previousL > 2.0:
            assert len(currentSum) > 0
            constraints = constraints + [sum(currentSum) == len(currentSum) - 1]
            currentSum = [boolvars[i]]
        elif el != 0:
            currentSum = currentSum + [boolvars[i]]
        previousL = el
    if len(currentSum) > 1:
        constraints = constraints + [sum(currentSum) == len(currentSum) - 1]
    obj = cvxpy.Minimize(0.0)
    prob = cvxpy.Problem(obj, constraints)
    prob.solve(solver=cvxpy.CBC, verbose=False)
    t_end = clock()
    if prob.status not in ["infeasible", "unbounded"]:
        res = np.array([v.value for v in variables])
        return ResultData(True, ms(t_end - t_init), res)
    else:
        return ResultData(False, ms(t_end - t_init))


# -------------------------- L1-Norm minimisation --------------------------------------------------
