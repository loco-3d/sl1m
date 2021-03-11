import numpy as np


from sl1m.constants_and_tools import *
from sl1m import planner_l1 as pl1
from sl1m import planner    as pl

from sl1m import qp


np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})



### This solver is called when the sparsity is fixed. It assumes the first contact surface for each phase
### is the one used for contact creation.
def solve(pb,surfaces, draw_scene = None, plot = True ):  
        
    t1 = clock()
    A, b, E, e = pl.convertProblemToLp(pb)    
    C = identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    res = qp.quadprog_solve_qp(C,c,A,b,E,e)
    t3 = clock()
    
    print("time to set up problem" , timMs(t1,t2))
    print("time to solve problem"  , timMs(t2,t3))
    print("total time"             , timMs(t1,t3))
    
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
        
    res = qp.gurobi_solve_lp(c,A,b,E,e)
        
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
                res = qp.gurobi_solve_lp(c,A,b,E,e)
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
                print("unfeasible problem")
                pass
            
        t4 = clock()      
        
        print("time to solve combinatorial ", timMs(t3,t4))
    
    if ok:
        surfacesret, indices = pl1.bestSelectedSurfaces(pb, res)        
        for i, phase in enumerate(pb["phaseData"]): 
            phase["S"] = [surfaces[i][indices[i]]]
        if solutionIndices is not None:
            for i, idx in enumerate(solutionIndices):
                pb["phaseData"][idx]["S"] = [surfaces[idx][solutionComb[i]]]
        
        return solve(pb,surfaces, draw_scene = draw_scene, plot = True )  


############### MIXED-INTEGER SOLVER ###############

def solveMIP(pb, surfaces, draw_scene = None, plot = True):  
    A, b, E, e = pl1.convertProblemToLp(pb)   
    c = pl1.slackSelectionMatrix(pb)
    
    t1 = clock()
    res = qp.gurobi_solve_mip(c,A,b,E,e)
    t2 = clock()

    print("time to solve MIP ", timMs(t1,t2))
    
    plot = plot and draw_scene is not None 
    if plot:
        ax = draw_scene(surfaces)
        pl1.plotQPRes(pb, res, ax=ax)
    
    coms, footpos, allfeetpos = pl1.retrieve_points_from_res(pb, res)
    
    return pb, coms, footpos, allfeetpos, res
        
