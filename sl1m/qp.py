import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity, ones

from scipy.optimize import linprog

# import gurobi and cvxpy
GUROBI_OK = False
try:
    import gurobipy as grb
    import cvxpy as cp
    GUROBI_OK = True
    grb.setParam('LogFile', '')
    grb.setParam('OutputFlag', 0)

except ImportError:
    pass


############### QUADRATIC PROGRAM SOLVER ###############

#min (1/2)x' P x + q' x
#subject to  G x <= h
#subject to  C x  = d
def quadprog_solve_qp(P, q, G=None, h=None, C=None, d=None, verbose = False):
    #~ qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if C is not None:
        if G is not None:
            qp_C = -vstack([C, G]).T
            qp_b = -hstack([d, h])
        else:
            qp_C = -C.transpose()
            qp_b = -d
        meq = C.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    res = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)
    if verbose:
        return res
    #~ print 'qp status ', res
    return res[0]


#min ||Ax-b||**2
#subject to  G x <= h
#subject to  C x  = d
def solve_least_square(A,b,G=None, h=None, C=None, d=None):
    P = dot(A.T, A)
    #~ q = 2*dot(b, A).reshape(b.shape[0])
    q = -dot(A.T, b).reshape(b.shape[0])
    #~ q = 2*dot(b, A)
    return solve_qp_quadprog(P, q, G, h, C, d)

#min x' C x + c x
#subject to  A x <= b
#subject to  E x  = e
def gurobi_solve_qp(C, c, A=None, b=None, E=None, e=None):
    if not GUROBI_OK:
        print("Requires gurobi")
        raise ImportError
    model = grb.Model("lp")

    #add continuous variables
    cVars = []
    for el in (c):
        cVars.append(model.addVar(vtype=grb.GRB.CONTINUOUS,  lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))

    # Update model to integrate new variables
    model.update()
    x = array(model.getVars(), copy=False)

    # equality constraints
    if E.shape[0] > 0:
        for i in range(E.shape[0]):
            idx = [j for j, el in enumerate(E[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = E[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, e[i])
    model.update()

    # inequality constraints
    if A.shape[0] > 0:
        for i in range(A.shape[0]):
            idx = [j for j, el in enumerate(A[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = A[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, b[i])
    model.update()

    # set objective
    obj = grb.QuadExpr()
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if C[i][j] != 0:
                obj += C[i][j]*cVars[i]*cVars[j]
    for el in c:
        if el != 0:
            obj += el*cVars[j]
    model.setObjective(obj, grb.GRB.MINIMIZE)
    model.optimize()

    try:
        res = [el.x for el in cVars]
        return res
    except:
        return None

############### LINEAR PROGRAM SOLVER ###############

#min q' x
#subject to  G x <= h
#subject to  C x  = d
def solve_lp(q, G=None, h=None, C=None, d=None): 
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=[(-100000.,10000.) for _ in range(q.shape[0])], method='interior-point', callback=None, options={'presolve': True})
    return res

#min sum (alpha)
#subject to  A x <= b + alpha
#subject to  E x  = e + alpha
def gurobi_solve_lp(c, A=None, b=None, E=None, e=None):
    if not GUROBI_OK:
        print("Requires gurobi")
        raise ImportError
    model = grb.Model("lp")

    #add continuous variables
    cVars = []
    for el in (c):
        cVars.append(model.addVar(vtype=grb.GRB.CONTINUOUS,  lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))

    # Update model to integrate new variables
    model.update()
    x = array(model.getVars(), copy=False)

    # equality constraints
    if E.shape[0] > 0:
        for i in range(E.shape[0]):
            idx = [j for j, el in enumerate(E[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = E[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, e[i])
    # model.update()

    # inequality constraints
    if A.shape[0] > 0:
        for i in range(A.shape[0]):
            idx = [j for j, el in enumerate(A[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = A[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, b[i])


    slackIndices = [i for i,el in enumerate (c) if el > 0]

    # set objective
    variables = []
    previousL = 0
    obj = 0
    for i, el in enumerate(slackIndices):
        if i != 0 and el - previousL > 2.:
            assert len(variables) > 0
            expr = grb.LinExpr(ones(len(variables)), variables)
            obj += expr
            variables = [x[el]]
        elif el!=0:
            variables += [x[el]]
        previousL = el
    if len(variables) > 1:
        expr = grb.LinExpr(ones(len(variables)), variables)
        obj += expr

    model.setObjective(obj, grb.GRB.MINIMIZE)
    model.optimize()

    try:
        res = [el.x for el in cVars]
        return res
    except:
        return None

############### MIXED-INTEGER PROGRAM SOLVER ###############

#subject to  A x <= b + a
#subject to  E x  = e + a
#subject to  sum(a) = m -1
def gurobi_solve_mip(c, A=None, b=None, E=None, e=None):
    if not GUROBI_OK:
        print("Requires gurobi")
        raise ImportError

    model = grb.Model("mip")

    cVars = []
    for el in (c):
        if el > 0:
            cVars.append(model.addVar(lb=0, ub=1, vtype=grb.GRB.BINARY, name = 'slack'))
        else:
            cVars.append(model.addVar(lb=-grb.GRB.INFINITY, ub=grb.GRB.INFINITY, vtype=grb.GRB.CONTINUOUS, name = 'x'))

    # Update model to integrate new variables
    model.update()
    x = array(model.getVars(), copy=False)

    # inequality constraints
    if A.shape[0] > 0:
        for i in range(A.shape[0]):
            idx = [j for j, el in enumerate(A[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = A[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.LESS_EQUAL, b[i])
    model.update()

    # equality constraints
    if E.shape[0] > 0:
        for i in range(E.shape[0]):
            idx = [j for j, el in enumerate(E[i].tolist()) if el != 0.]
            variables = x[idx]
            coeff = E[i,idx]
            expr = grb.LinExpr(coeff, variables)
            model.addConstr(expr, grb.GRB.EQUAL, e[i])
    model.update()

    slackIndices = [i for i,el in enumerate (c) if el > 0]

    # equality slack sum
    variables = []
    previousL = 0
    for i, el in enumerate(slackIndices):
        if i != 0 and el - previousL > 2.:
            assert len(variables) > 0
            expr = grb.LinExpr(ones(len(variables)), variables)
            model.addConstr(expr, grb.GRB.EQUAL, len(variables) -1)
            variables = [x[el]]
        elif el!=0:
            variables += [x[el]]
        previousL = el
    if len(variables) > 1:
        expr = grb.LinExpr(ones(len(variables)), variables)
        model.addConstr(expr, grb.GRB.EQUAL, len(variables) -1)
    model.update()
    model.optimize()
    try:
        res = [el.x for el in cVars]
        return res
    except:
        return None

def tovals(variables):
    return array([el.value for el in variables])

#subject to  A x <= b + alpha
#subject to  E x  = e + alpha
#subject to  alpha < 100 * a
#subject to  sum(a) = m -1
def cvxpy_solve_mip (c,A=None, b=None, E=None, e=None):
    if not GUROBI_OK:
        print("Requires gurobi packaged in cvxpy")
        raise ImportError

    #real variables
    varReal = cp.Variable(len(c))
    #boolean vars
    slackIndices = [i for i,el in enumerate (c) if el > 0]
    numSlackVariables = len([el for el in c if el > 0])
    boolvars = cp.Variable(numSlackVariables, boolean=True)

    #constraints
    constraints = []
    constraintNormalIneq = A @ varReal <= b
    constraintNormalEq   = E @ varReal == e
    constraints = [constraintNormalIneq, constraintNormalEq]
    constraints = constraints + [varReal[el] <= 100. * boolvars[i] for i, el in enumerate(slackIndices)]

    currentSum = []
    previousL = 0
    for i, el in enumerate(slackIndices):
        if i!= 0 and el - previousL > 2.:
            assert len(currentSum) > 0
            constraints = constraints + [sum(currentSum) == len(currentSum) -1 ]
            currentSum = [boolvars[i]]
        elif el !=0:
            currentSum = currentSum + [boolvars[i]]
        previousL  = el
    if len(currentSum) > 1:
        constraints = constraints + [sum(currentSum) == len(currentSum) -1 ]

    obj = cp.Minimize(0)
    prob = cp.Problem(obj, constraints)
    res = prob.solve(solver=cp.GUROBI, verbose=False )
    res = tovals(varReal)
    return res


if __name__ == '__main__':

    from numpy.linalg import norm

    A = array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
    b = array([3., 2., 3.])
    P = dot(A.T, A)
    q = 2*dot(b, A).reshape((3,))
    G = array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
    h = array([3., 2., -2.]).reshape((3,))

    res2 = solve_least_square(A, b, G, h)
    res1 =  quadprog_solve_qp(P, q, G, h)
    print(res1)
    print(res2)
