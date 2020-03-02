import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity, abs, ones

from scipy.optimize import linprog


__eps = 0.00001

GLPK_OK = False  
try:
    import glpk
    GLPK_OK = True

except ImportError:
    pass

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
    def __init__(self, x, status, success, cost):
        self.x = x
        self.status = status
        self.success = success
        self.cost = cost
        
    def __str__(self):
        return "ResultData: \n \t solver status: " + str(self.status) + "\n \t success: " + str(self.success) + "\n \t x: " + str(self.x) + "\n \t cost: " + str(self.cost)
   
    def __repr__(self):
        return self.__str__()
            

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
    try:
        res = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)
        return ResultData(res[0],  'opt', True, res[1])
    except:
        return ResultData(-ones(P.shape[1]),  'unfeasible', False, 0.)
        


#min ||Ax-b||**2 
#subject to  G x <= h
#subject to  C x  = d
def solve_least_square(A,b,G=None, h=None, C=None, d=None):
        P = dot(A.T, A)
        #~ q = 2*dot(b, A).reshape(b.shape[0])
        q = -dot(A.T, b).reshape(b.shape[0])
        return quadprog_solve_qp(P, q, G, h, C, d)

#min q' x  
#subject to  G x <= h
#subject to  C x  = d
def solve_lp(q, G=None, h=None, C=None, d=None): 
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=[(-100000.,10000.) for _ in range(q.shape[0])], method='interior-point', callback=None, options={'presolve': True})
    return ResultData(res['x'],  res['status'], res['success'], res['fun'])
        
if GLPK_OK:    
    
    #solve linear programm using the simplex method with glpk
    # min q' x  
    #subject to  CI x <= ci0
    #subject to  CE x  = ce0
    def solve_lp_glpk(q, CI=None, ci0=None, CE=None, ce0=None):
        lp = glpk.LPX() 
        # ~ lp.name = 'sample'
        lp.obj.maximize = False
        
        numEqConstraints = 0
        numIneqConstraints = 0
        if CE is not None:
            numEqConstraints   = CE.shape[0];
            xsize = CE.shape[1];  
        if CI is not None:
            numIneqConstraints = CI.shape[0];
            xsize = CI.shape[1];  
        
        numConstraints = numEqConstraints + numIneqConstraints
        
        lp.cols.add(xsize) 
        
        for c in lp.cols:      # Iterate over all columns
            c.bounds = None, None
        
        if numConstraints > 0:        
            lp.rows.add(numConstraints) 
        
            idrow = 0;
            idcol = 0;
            idConsMat = 1;
            
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
        from time import clock
        t1 = clock()
        lp.simplex()
        t2 = clock()
        return ResultData(array([c.primal for c in lp.cols]), lp.status, lp.status == "opt", lp.obj.value)


if GUROBI_OK:  
    #solve linear programm using gurobi
    # min q' x  
    #subject to  A x <= b
    #subject to  E x  = e
    def solve_lp_gurobi(c, A=None, b=None, E=None, e=None):
        
        model = grb.Model("lp")
        
        rdim = A.shape[1]
        
        #add continuous variables
        cVars = []
        for el in (c):
            cVars.append(model.addVar(obj = el, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
        
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
            
        model.modelSense = grb.GRB.MINIMIZE
        model.optimize()
        res = [el.x for el in cVars]
        return ResultData(res, model.Status, model.Status == grb.GRB.OPTIMAL, model.ObjVal)
        
        
        

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

        A = array([[-1.,0., 0.], [-0., -1., 0.], [0., 0., -1.]])
        b = array([-1., 2., 3.])
        C = array([1.,1.,1.]).reshape([1,-1])
        d = array([1.])
        reslp = solve_lp(b, A, b,C,d)
        resglpk = solve_lp_glpk(b, A, b,C,d)
        resgurobi = solve_lp_gurobi(b, A, b,C,d)

        print("lp", reslp.x, b.dot(reslp.x))
        print("lpk", resglpk.x,  b.dot(resglpk.x))
        print("gurobi", resgurobi.x,  b.dot(resgurobi.x))
