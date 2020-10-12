from numpy import array, dot, vstack, hstack, asmatrix, identity, ones
from time import clock

import quadprog
from scipy.optimize import linprog

def timMs(t1, t2):
    return (t2-t1) * 1000.

__eps = 0.00001

# import glpk
GLPK_OK = False  
try:
    import glpk
    GLPK_OK = True

except ImportError:
    print ("import error")
    pass

# import gurobipy
GUROBI_OK = False
try:
    import gurobipy
    GUROBI_OK = True
    import gurobipy as grb
    grb.setParam('LogFile', '')
    grb.setParam('OutputFlag', 0)

except ImportError:
    pass

class ResultData:
    def __init__(self, x, status, success, cost, time):
        self.x = x
        self.status = status
        self.success = success
        self.cost = cost
        self.time = time
        
    def __str__(self):
        return "ResultData: \n \t solver status: " + str(self.status) + "\n \t success: " + str(self.success) + "\n \t x: " + str(self.x) + "\n \t cost: " + str(self.cost)
   
    def __repr__(self):
        return self.__str__()

def get_nonzero_rows(M):
    nonzero_rows = {}
    rows, cols = M.nonzero()
    for ij in zip(rows, cols):
        i, j = ij
        if i not in nonzero_rows:
            nonzero_rows[i] = []
        nonzero_rows[i].append(j)
    return nonzero_rows
    
#min (1/2)x' P x + q' x  
#subject to  G x <= h
#subject to  C x  = d
def quadprog_solve_qp(P, q, G=None, h=None, C=None, d=None, verbose = False):
    #~ qp_G = .5 * (P + P.T)   # make sure P is symmetric
    t1 = clock()
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
    try:
        res = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)
        t2 = clock()
        return ResultData(res[0], 'opt', True, res[1], timMs(t1, t2))
    except:
        t2 = clock()
        return ResultData(None,  'unfeasible', False, 0., timMs(t1, t2))


#min ||Ax-b||**2 
#subject to  G x <= h
#subject to  C x  = d
def solve_least_square(A,b,G=None, h=None, C=None, d=None):
        P = dot(A.T, A)
        #~ q = 2*dot(b, A).reshape(b.shape[0])
        q = -dot(A.T, b).reshape(b.shape[0])
        #~ q = 2*dot(b, A) 
        return quadprog_solve_qp(P, q, G, h, C, d)


#min q' x  
#subject to  G x <= h
#subject to  C x  = d
def solve_lp(q, G=None, h=None, C=None, d=None): 
    t1 = clock()
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=[(-100000.,10000.) for _ in range(q.shape[0])], method='interior-point', callback=None, options={'presolve': True})
    t2 = cock()
    return ResultData(res['x'],  res['status'], res['success'], res['fun'], timMs(t1, t2))

    
if GLPK_OK:    
    
    #solve linear programm using the simplex method with glpk
    # min q' x  
    #subject to  CI x <= ci0
    #subject to  CE x  = ce0
    def solve_lp_glpk(q, CI=None, ci0=None, CE=None, ce0=None):
        
        t1 = clock()
        lp = glpk.LPX() 
        # ~ lp.name = 'sample'
        lp.obj.maximize = False
        
        numEqConstraints = 0
        numIneqConstraints = 0
        if CE is not None:
            numEqConstraints   = CE.shape[0]
            xsize = CE.shape[1]
        if CI is not None:
            numIneqConstraints = CI.shape[0]
            xsize = CI.shape[1]
        
        numConstraints = numEqConstraints + numIneqConstraints
        
        lp.cols.add(xsize) 
        
        for c in lp.cols:      # Iterate over all columns
            c.bounds = None, None
        
        if numConstraints > 0:        
            lp.rows.add(numConstraints) 
        
            idrow = 0
            idcol = 0
            idConsMat = 1
            
            mat = []
            
            for i in range(numIneqConstraints):
                lp.rows[i].bounds = None, ci0[i]
                for j in range(xsize):
                     if abs(CI[i,j]) > __eps:
                         mat.append((idrow, j, CI[i,j]))
                idrow+=1
                
            for i in range(numEqConstraints):
                lp.rows[idrow].bounds = ce0[i]
                for j in range(xsize):
                     if abs(CE[i,j]) > __eps:
                         mat.append((idrow, j, CE[i,j]))
                idrow+=1
        
            lp.matrix = mat
        lp.obj[:] = q.tolist()
        lp.simplex()
        t2 = clock()
        print (lp.obj.value)
        return ResultData(array([c.primal for c in lp.cols]), lp.status, lp.status == "opt", lp.obj.value, timMs(t1, t2))

if GUROBI_OK:  

    def solve_qp_gurobi(P, q, G=None, h=None, A=None, b=None, verbose=False):
    
        t1 = clock()
        grb.setParam('OutputFlag', 1 if verbose else 0)
        n = P.shape[1]
        model = grb.Model("qp")

        cVars = []
        for i in range(n):
            cVars.append(model.addVar(vtype=grb.GRB.CONTINUOUS,  lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
        
        # Update model to integrate new variables
        model.update()
        x = array(model.getVars(), copy=False)

        # minimize
        #     1/2 x.T * P * x + q * x
        obj = grb.QuadExpr()
        rows, cols = P.nonzero()
        for i, j in zip(rows, cols):
            obj += 0.5 * x[i] * P[i, j] * x[j]
        for i in range(n):
            obj += q[i] * x[i]
        model.setObjective(obj, grb.GRB.MINIMIZE)

        # subject to
        #     G * x <= h
        if G is not None:
            G_nonzero_rows = get_nonzero_rows(G)
            for i, row in G_nonzero_rows.items():
                model.addConstr(grb.quicksum(G[i, j] * x[j] for j in row) <= h[i])

        # subject to
        #     A * x == b
        if A is not None:
            A_nonzero_rows = get_nonzero_rows(A)
            for i, row in A_nonzero_rows.items():
                model.addConstr(grb.quicksum(A[i, j] * x[j] for j in row) == b[i])

        model.optimize()

        t2 = clock()
        try:
            res = [el.x for el in cVars]
            print (model.ObjVal)
            return ResultData(res, model.Status, model.Status == grb.GRB.OPTIMAL, model.ObjVal, timMs(t1, t2))
        except:
            return ResultData(None,  model.Status, False, 0., timMs(t1, t2))

    #solve linear programm using gurobi
    # min q' x  
    #subject to  A x <= b
    #subject to  E x  = e
    def solve_lp_gurobi(c, A=None, b=None, E=None, e=None):
        t1 = clock()
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

        
        # slackIndices = [i for i,el in enumerate (c) if el > 0]
            
        # # set objective
        # variables = []
        # previousL = 0
        # obj = 0
        # for i, el in enumerate(slackIndices):
        #     if i != 0 and el - previousL > 2.:
        #         assert len(variables) > 0
        #         expr = grb.LinExpr(ones(len(variables))/(len(variables)-1), variables)
        #         # expr = grb.LinExpr(ones(len(variables)), variables)
        #         obj += expr
        #         variables = [x[el]]
        #     elif el!=0:
        #         variables += [x[el]]
        #     previousL = el
        # if len(variables) > 1:
        #     expr = grb.LinExpr(ones(len(variables))/(len(variables)-1), variables)
        #     # expr = grb.LinExpr(ones(len(variables)), variables)
        #     obj += expr 
            
        model.setObjective(obj, grb.GRB.MINIMIZE)
        # model.modelSense = grb.GRB.MINIMIZE
        # model.update()            
        model.optimize()
        
        t2 = clock()
        try:
            res = [el.x for el in cVars]
            print (model.ObjVal)
            return ResultData(res, model.Status, model.Status == grb.GRB.OPTIMAL, model.ObjVal, timMs(t1, t2))
        except:
            return ResultData(None,  model.Status, False, 0., timMs(t1, t2))


    def solve_MIP_gurobi(c, A=None, b=None, E=None, e=None):
        t1 = clock()
        model = grb.Model("mip")
        
        cVars = []
        for el in (c):
            if el > 0:
                cVars.append(model.addVar(lb=0, ub=1, vtype=grb.GRB.BINARY, name = 'slack'))
            else:
                cVars.append(model.addVar(lb=-grb.GRB.INFINITY, ub=grb.GRB.INFINITY, vtype=grb.GRB.CONTINUOUS, name = 'x'))
        #     cVars.append(model.addVar(obj = el, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY, name = 'x'))
        
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
        t2 = clock()
        try:
            res = [el.x for el in cVars]
            return ResultData(res, model.Status, model.Status == grb.GRB.OPTIMAL, model.ObjVal, timMs(t1, t2))
        except:
            return ResultData(None,  model.Status, False, 0., timMs(t1, t2))    
        

if __name__ == '__main__':
        
    from numpy.linalg import norm
    
    # A = array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
    # b = array([3., 2., 3.])
    # P = dot(A.T, A)
    # q = 2*dot(b, A).reshape((3,))
    # G = array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
    # h = array([3., 2., -2.]).reshape((3,))

    # res2 = solve_least_square(A, b, G, h)
    # res1 =  quadprog_solve_qp(P, q, G, h)
    # print(res1)
    # print(res2)

    A = array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    b = array([1, 0, 3])
    C = array([[1, 1, 1]])
    d = array([1.])
    
    resglpk = solve_lp_glpk(b, A, b, C, d)
    resgurobi = solve_lp_gurobi(b, A, b, C, d)
    print (resglpk.x  )
    print (resgurobi.x)
