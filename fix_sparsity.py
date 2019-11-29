import numpy as np


from sl1m.constants_and_tools import *
from sl1m import planner_l1 as pl1
from sl1m import planner    as pl

import qp


# try to import mixed integer solver
MIP_OK = False  
try:
    import gurobipy
    import cvxpy as cp
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
    # b += ones(b.shape)  
    C = identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    # res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    res = qp.solve_lp(c,A,b,E,e)
    t3 = clock()
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    if type(res) is int:
        print "LP fails"
        return False
    
    coms, footpos, allfeetpos = pl.retrieve_points_from_res(pb, res)
    
    plot = plot and draw_scene is not None 
    if plot:
        ax = draw_scene(surfaces)
        pl.plotQPRes(pb, res, ax=ax)
    
    return pb, coms, footpos, allfeetpos, res


### Calls the sl1m solver. Brute-forcedly tries to solve non fixed sparsity by handling the combinatorial.
### Ultimately calls solve which provides the approriate cost function
def solveL1(pb, surfaces, draw_scene = None, plot = True, convert = True):     
    A, b, E, e = pl1.convertProblemToLp(pb, convert)    
    C = identity(A.shape[1]) * 0.00001
    c = pl1.slackSelectionMatrix(pb)
        
    res = qp.quadprog_solve_qp(C, c,A,b,E,e)
        
    ok = pl1.isSparsityFixed(pb, res)
    solutionIndices = None
    solutionComb = None
    if not ok:
        print "not ok"
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
        print "ok ?", ok
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

def solveMIP(pb, surfaces, MIP = True, draw_scene = None, plot = True, convert = True):  
    if not MIP_OK:
        print "Mixed integer formulation requires gurobi packaged in cvxpy"
        raise ImportError
        
    gurobipy.setParam('LogFile', '')
    gurobipy.setParam('OutputFlag', 0)
       
    A, b, E, e = pl1.convertProblemToLp(pb, convert)   
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
                currentSum = []
            elif el !=0:
                currentSum = currentSum + [boolvars[i]]
            previousL  = el
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
    
    return pb, res#timMs(t1,t2)
        
