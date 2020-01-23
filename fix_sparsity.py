import numpy as np


from sl1m.constants_and_tools import *
from sl1m import planner_l1 as pl1
from sl1m import planner    as pl

import qp


# try to import mixed integer solver
MIP_OK = False  
try:
    #~ from gurobipy import *
    import gurobipy as grb
    import cvxpy as cp
    from scipy.sparse import csr_matrix
    MIP_OK = True

except ImportError:
    pass


from time import clock

np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})



### This solver is called when the sparsity is fixed. It assumes the first contact surface for each phase
### is the one used for contact creation.
def solve(pb,surfaces, draw_scene = None, plot = True ):  
        
    t1 = clock()
    A, b, E, e = pl.convertProblemToLp(pb)    
    C = identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    t3 = clock()
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = pl.retrieve_points_from_res(pb, res)
    
    plot = plot and draw_scene is not None 
    if plot:
        ax = draw_scene(surfaces)
        pl.plotQPRes(pb, res, ax=ax)
    
    return pb, coms, footpos, allfeetpos, res


### Calls the sl1m solver. Brute-forcedly tries to solve non fixed sparsity by handling the combinatorial.
### Ultimately calls solve which provides the approriate cost function
def solveL1(pb, surfaces, draw_scene = None, plot = True):     
    A, b, E, e = pl1.convertProblemToLp(pb)    
    C = identity(A.shape[1]) * 0.00001
    c = pl1.slackSelectionMatrix(pb)
    
    
    t1 = clock()
    res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    t2 = clock()
    
    print "time to solve initial qp ", timMs(t1,t2)
        
    ok = pl1.isSparsityFixed(pb, res)
    solutionIndices = None
    solutionComb = None
    if not ok:
        pbs = pl1.generateAllFixedScenariosWithFixedSparsity(pb, res)
        
        t3 = clock()
        
        for (pbComb, comb, indices) in pbs:
            A, b, E, e = pl1.convertProblemToLp(pbComb, convertSurfaces = False)
            C = identity(A.shape[1]) * 0.00001
            c = pl1.slackSelectionMatrix(pbComb)
            try:
                res = qp.quadprog_solve_qp(C, c,A,b,E,e)
                if pl1.isSparsityFixed(pbComb, res):       
                    coms, footpos, allfeetpos = pl1.retrieve_points_from_res(pbComb, res)
                    pb = pbComb
                    ok = True
                    solutionIndices = indices[:]
                    solutionComb = comb
                    if plot:
                        ax = draw_scene(surfaces)
                        pl1.plotQPRes(pb, res, ax=ax)
                    break
            except:
                print "unfeasible problem"
                pass
            
        t4 = clock()      
        
        print "time to solve combinatorial ", timMs(t3,t4)
    
    if ok:
        surfacesret, indices = pl1.bestSelectedSurfaces(pb, res)        
        for i, phase in enumerate(pb["phaseData"]): 
            phase["S"] = [surfaces[i][indices[i]]]
        if solutionIndices is not None:
            for i, idx in enumerate(solutionIndices):
                pb["phaseData"][idx]["S"] = [surfaces[idx][solutionComb[i]]]
        
        return solve(pb,surfaces, draw_scene = draw_scene, plot = True )  


############### MIXED-INTEGER SOLVER ###############

def tovals(variables):
    return array([el.value for el in variables])



def solveMIPcvx(pb, surfaces, MIP = True, draw_scene = None, plot = True):  
    if not MIP_OK:
        print "Mixed integer formulation requires gurobi packaged in cvxpy"
        raise ImportError
        
    grb.setParam('LogFile', '')
    grb.setParam('OutputFlag', 0)
       
    A, b, E, e = pl1.convertProblemToLp(pb)   
    slackMatrix = pl1.slackSelectionMatrix(pb)
    
    rdim = A.shape[1]
    varReal = cp.Variable(rdim)
    constraints = []
    constraintNormalIneq = A * varReal <= b
    constraintNormalEq   = E * varReal == e
    
    constraints = [constraintNormalIneq, constraintNormalEq]
    #creating boolean vars
    
    slackIndices = [i for i,el in enumerate (slackMatrix) if el > 0]
    numSlackVariables = len([el for el in slackMatrix if el > 0])
    boolvars = cp.Variable(numSlackVariables, boolean=True)    
    obj = cp.Minimize(slackMatrix * varReal)
    
    if MIP:    
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
            
    obj = cp.Minimize(ones(numSlackVariables) * boolvars)
    prob = cp.Problem(obj, constraints)
    t1 = clock()
    res = prob.solve(solver=cp.GUROBI, verbose=False )
    t2 = clock()
    res = tovals(varReal)
    print "time to solve MIP ", timMs(t1,t2)

    
    plot = plot and draw_scene is not None 
    if plot:
        ax = draw_scene(surfaces)
        pl1.plotQPRes(pb, res, ax=ax)
    
    return timMs(t1,t2)
        

def solveMIP(pb, surfaces, MIP = True, draw_scene = None, plot = True, initGuess = None, initGuessMip = None):  
    if not MIP_OK:
        print "Mixed integer formulation requires gurobi packaged in cvxpy"
        raise ImportError
        
    grb.setParam('LogFile', '')
    grb.setParam('OutputFlag', 0)
    #~ grb.setParam('OutputFlag', 1)
       
    A, b, E, e = pl1.convertProblemToLp(pb)   
    slackMatrix = pl1.slackSelectionMatrix(pb)    
    slackIndices = [i for i,el in enumerate (slackMatrix) if el > 0]
    numSlackVariables = len([el for el in slackMatrix if el > 0])
    
    
    model = grb.Model("mip")
    
    rdim = A.shape[1]
    
    #add continuous variables
    cVars = []
    for i in range(rdim):
        if slackMatrix[i] > 0:
            if MIP:
                cVars.append(model.addVar(name="slack%d" % i, obj = 0, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
            else:
                cVars.append(model.addVar(name="slack%d" % i, obj = 1, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
        else:
            cVars.append(model.addVar(name="c%d" % i, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
            
    
    
    # Update model to integrate new variables
    model.update()
    
    x = np.array(model.getVars(), copy=False)
    
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
    
    
    if MIP:    
                
        #create boolean variables  
        boolvars = []
        for i in range(numSlackVariables):
            boolvars.append(model.addVar(vtype=grb.GRB.BINARY,
                                 obj=0,
                                 name="boolVar%d" % i))
        model.update()
        #Big M value
        M = 1000.
        [model.addConstr(cVars[el] == M * boolvars[i], "boolAlpha%d" % i ) for i, el in enumerate(slackIndices)]  
    
        model.update()
    
        currentSum = []
        currentSum2 = []
        previousL = 0
        for i, el in enumerate(slackIndices):
            if i!= 0 and el - previousL > 2.:
                model.addConstr(grb.quicksum(currentSum) == len(currentSum) -1, "card%d" % i)
                assert len(currentSum) > 0      
                currentSum = [boolvars[i]]
                currentSum2 = [el]          
            elif el !=0:
                currentSum = currentSum + [boolvars[i]]
                currentSum2 = currentSum2 + [el]
            previousL  = el        
        if len(currentSum) > 1:
            model.addConstr(grb.quicksum(currentSum) == len(currentSum) -1, "card%d" % i)
    model.modelSense = grb.GRB.MINIMIZE
    
    if initGuess is not None:
        for (i,el) in initGuess:
            x[i].start = el
    
    if MIP and initGuessMip is not None:
        for (i,el) in initGuessMip:
            boolvars[i].start = el
        
    
    model.update()
    t1 = clock()
    model.optimize()
    t2 = clock()
    res = [el.x for el in cVars]
    print "time to solve MIP ", timMs(t1,t2)
    
    plot = plot and draw_scene is not None 
    if plot:
        ax = draw_scene(surfaces)
        pl1.plotQPRes(pb, res, ax=ax)
    
    #~ return timMs(t1,t2)
    
    
    if MIP:
        return res, [el.x for el in boolvars]
    else:
        return res
        
