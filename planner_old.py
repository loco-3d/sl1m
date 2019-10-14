import numpy as np


from mpcroc.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *

from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray
import qp



NUM_COM_PER_PHASES = 1

############### Problem definition #############
from mpcroc.problem_definition import *

def getTotalNumVariablesAndIneqConstraints(pb):
    cols    = sum([phase["S"][0].shape[1] + NUM_COM_PER_PHASES for  phase in pb["phaseData"] ])
    rows    = sum( [sum([k.shape[0]  for (_, k) in phase["K"][0] ]) + phase["relativeK"][0][1].shape[0] for  phase in pb["phaseData"] ])
    rowsOld = sum( [sum([k.shape[0]  for (_, k) in phase["K"][0] ])  for  phase in pb["phaseData"] ])
    rows += cols #positivity constraint
    return rows, cols
    
def getTotalNumEqualityConstraints(pb):
    return pb["nphases"]
  
def numVariablesForPhase(phase):
    return phase["S"][0].shape[1] + NUM_COM_PER_PHASES
      
# COM position is expressed as convex combination from contact surfaces for x and y,
# while z is directly given by a variable
def currentCOMExpressionMatrix(phaseDataT, startCol, endCol, comIdx = 0):
    lenVar = endCol - startCol
    ret = zeros((3, lenVar))
    #surfaces
    ret[:2,:-1] = phaseDataT["S"][0][:2,:]
    ret[-1, -NUM_COM_PER_PHASES + comIdx] = 1.  #retrieve com z 
    
    return ret
    
    
# Foot position matrix
# used to retrieve foot position in the previous phase essentially
def currentFootExpressionMatrix(phaseDataT, startCol, endCol):
    lenVar = endCol - startCol
    ret = zeros((3, lenVar))
    ret[:,:phaseDataT["S"][0].shape[1]] = phaseDataT["S"][0][:,:]    
    return ret
    
# z position relative to support foot
def currentCOMzmatrix(phaseDataT, startCol, endCol, comIdx = 0):
    ret = zeros((1, endCol - startCol))
    ret[-1,phaseDataT["S"][0].shape[1] + comIdx] = 1.
    return ret
    
    
def FixedFootKinConstraintInitPhase(pb, phaseDataT, comZMatrix, fixedFootMatrix, A, b, startCol, endCol):
    fixed = phaseDataT["fixed"]
    pos   = pb["p0"][fixed]
    K, k = phaseDataT["K"][0][fixed]
    idRow = K.shape[0]
    #TODO 
    A[:idRow, startCol:endCol] = K[:,-1:].dot(comZMatrix)
    b[:idRow                 ] = k + K[:,-1:].dot(pos[-1:])
    return idRow
    
def FixedFootKinConstraintVarPhase(pb, phaseDataT, comZMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["K"][0][fixed]
    idRow = startRow + K.shape[0]
    A[startRow:idRow, startCol:endCol] = K[:,-1:].dot(comZMatrix)
    #prev variable size
    prevSize  = fixedFootMatrix.shape[1]
    A[startRow:idRow, (startCol-prevSize):startCol] = -K[:,-1:].dot(fixedFootMatrix[-1:,:])
    b[startRow:idRow                 ] = k 
    return idRow
    
    
def FixedFootConstraint(pb, phaseDataT, comZMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow):
    if startRow == 0:
        if pb["p0"] is None:
            fixed = phaseDataT["fixed"] #0 constraints
            K, k = phaseDataT["K"][0][fixed]
            return K.shape[0]
        else:
            return FixedFootKinConstraintInitPhase(pb, phaseDataT, comZMatrix, fixedFootMatrix, A, b, startCol, endCol)
    else:
        return FixedFootKinConstraintVarPhase (pb, phaseDataT, comZMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow)
     
     
def FootComConstraintInitPhase(pb, phaseDataT, comMatrix, footMatrix, targetFoot, A, b, endCol, startRow):
    target = phaseDataT[targetFoot]
    pos   = pb["p0"][target]
    K, k = phaseDataT["K"][0][target]
    idRow = K.shape[0] + startRow
        
    #no initial constraint on com
    if comMatrix is None:
        return idRow
    
    #initial constraint on com
    if len(comMatrix.shape) == 1 and comMatrix.shape[0] == 3: 
        A[startRow:idRow, endCol-footMatrix.shape[1]:endCol] -= K.dot(footMatrix)
        b[startRow:idRow                 ] = k - K.dot(comMatrix)
        return idRow

    
    A[startRow:idRow, endCol-comMatrix.shape[1]:endCol] = K.dot(comMatrix)
    b[startRow:idRow                 ] = k + K.dot(pos)
    return idRow
     
def FootComConstraintVarPhase(pb, phaseDataT, comMatrix, footMatrix, targetFoot, A, b, endCol, startRow):
    target = phaseDataT[targetFoot]
    K, k = phaseDataT["K"][0][target]
    idRow = K.shape[0] + startRow
    A[startRow:idRow, endCol-comMatrix.shape[1]:endCol] = K.dot(comMatrix)
    prevSize  = footMatrix.shape[1]
    A[startRow:idRow, (endCol-prevSize):endCol] -= K.dot(footMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
     
def FootComConstraint(pb, phaseDataT, comMatrix, footMatrix, targetFoot, A, b, endCol, startRow, first):
    if first:
        if pb["p0"] is None:
            target = phaseDataT[targetFoot] #0 constraints
            K, k = phaseDataT["K"][0][target]
            return startRow + K.shape[0]
        else:
            return FootComConstraintInitPhase(pb, phaseDataT, comMatrix, footMatrix, targetFoot, A, b, endCol, startRow)
    else:
        return FootComConstraintVarPhase (pb, phaseDataT, comMatrix, footMatrix, targetFoot, A, b, endCol, startRow)
        
def FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["relativeK"][0]
    footMatrix = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
    pos   = pb["p0"][fixed]
    idRow = K.shape[0]
    A[:idRow, startCol:endCol] = K.dot(footMatrix)
    b[:idRow                 ] = k + K.dot(pos)
    return idRow + startRow
    
def FixedFootRelativeDistanceConstraintVarPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow):
    K, k = phaseDataT["relativeK"][0]
    footMatrix = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
    idRow = startRow + K.shape[0]
    A[startRow:idRow, startCol:endCol] = K.dot(footMatrix)
    #prev variable size
    prevSize  = fixedFootMatrix.shape[1]
    A[startRow:idRow, (startCol-prevSize):startCol] = -K.dot(fixedFootMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
    
    
def FixedFootConstraintRelativeDistance(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow, first):
    if first:
        if True:
            K, k = phaseDataT["relativeK"][0] #0 constraints
            return startRow + K.shape[0]
        else:
            return FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow)
    else:
        return FixedFootRelativeDistanceConstraintVarPhase (pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow)        
    
#ensure alpha vector defines a convex combination 
def ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRow):
    ret = ones((1, endCol - startCol))
    ret[-1,-1] = 0
    idRow = startRow + 1
    E[startRow:idRow, startCol:endCol] = ret
    e[startRow:idRow                 ] = 1 
    return idRow

#com is foot position at last phase + z from this phase
def computePreviousCOMMatrix(pb,comZMatrix, fixedFootMatrix):
    if fixedFootMatrix is None:
        return pb["c0"]            
    ret = zeros((3,comZMatrix.shape[1]+fixedFootMatrix.shape[1]))
    ret[:2,:fixedFootMatrix.shape[1]] = fixedFootMatrix[:-1]
    ret[ 2, fixedFootMatrix.shape[1]:fixedFootMatrix.shape[1]+comZMatrix.shape[1]] = comZMatrix
    return ret
    
def computeCurrentCOMMatrix(pb,comZMatrix, fixedFootMatrix):
    if fixedFootMatrix is None:
        return pb["c0"]            
    ret = zeros((3,comZMatrix.shape[1]))
    ret[:2,:fixedFootMatrix.shape[1]] = fixedFootMatrix[:-1]
    ret[ 2, :] = comZMatrix
    return ret

def convertProblemToLp(pb):
    
    #assert that either an initial state cimplete exists or nothing exists
    assert (pb["c0"] is None and pb["p0"] is None) or (pb["c0"] is not None and pb["p0"] is not None), "initial state must be either fully described or not at all"
    
    #define first problem
    #A u <= b
    nIneq, nvars  = getTotalNumVariablesAndIneqConstraints(pb)
    A = zeros((nIneq, nvars)); b = zeros(nIneq)
    #E u = b
    nEq = getTotalNumEqualityConstraints(pb)
    E = zeros((nEq, nvars)); e = zeros(nEq)
    
    startRow = 0;
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    fixedFootMatrix = None;  
    
    #handle the initial state
    
    previousCOMMatrix  = None
    currentCOMMatrix   = None
    previousFootMatrix = None;  
    currentFootMatrix  = None;  
    
    startRowEq = 0;
    for i, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        nVarPhase = numVariablesForPhase(phaseDataT) 
        first = i == 0
        endCol = startCol + nVarPhase
        
        currentComZMatrix  = currentCOMzmatrix(phaseDataT, startCol, endCol)
        previousCOMMatrix  = computePreviousCOMMatrix(pb,currentComZMatrix, previousFootMatrix); 
        currentFootMatrix  = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
        currentCOMMatrix   = computeCurrentCOMMatrix(pb,currentComZMatrix, currentFootMatrix);        
        startRow = FixedFootConstraintRelativeDistance (pb, phaseDataT, previousFootMatrix, A, b, startCol, endCol, startRow, first) 
        
        # adding zeros to show that previous position is actually considered
        if previousFootMatrix is not None:
            previousFootMatrix = hstack([previousFootMatrix,zeros((previousFootMatrix.shape[0], nVarPhase))])
        startRow = FootComConstraint(pb, phaseDataT, currentCOMMatrix , previousFootMatrix, "fixed" , A, b, endCol, startRow, first)
        startRow = FootComConstraint(pb, phaseDataT, previousCOMMatrix, currentFootMatrix , "moving", A, b, endCol, startRow, first)
        
        #equality
        startRowEq = ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)
        
        # for next phase
        previousFootMatrix = currentFootMatrix
        startCol   = endCol 
        
    #add positivity constraint
    A[startRow:startRow+nvars, : ] = -identity(nvars)
    b[-nvars:      ] = zeros(nvars)
    
    
    startRow += nvars
    
    assert startRowEq == E.shape[0]
    #~ assert startRow   == A.shape[0]
    
    return (A,b,E,e)
    
def least_square_cost_function(pb, nvars, comTarget, reg_term = 0.0001):
    C = identity(nvars) * reg_term
    lastPhase = pb["phaseData"][-1]
    nvarPhase = numVariablesForPhase(lastPhase)
    lastCom = currentCOMExpressionMatrix(lastPhase, 0, nvarPhase)
    C[-3:,nvars-nvarPhase:] = lastCom
    return C, concatenate([zeros(nvars-3),comTarget])
    
###########################" PLOTTING ################"
    
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def retrieve_points_from_res(pb, res):
    coms = []
    if pb["c0"] is not None:
        coms = [pb["c0"]]
    if pb["p0"] is not None:
        footPos = [[pb["p0"][LF]],[pb["p0"][RF]]]
        allFeetPos = [footPos[0][0], footPos[1][0]]
    else:
        footPos = [[],[]]
        allFeetPos = []
    
    col = 0
    for i, phaseDataT in enumerate(pb["phaseData"]):  
        moving = phaseDataT["moving"]
        fixed = phaseDataT["fixed"]
        nvarPhase = numVariablesForPhase(phaseDataT)
        posMatrix = currentFootExpressionMatrix(phaseDataT, 0, nvarPhase)
        pos = posMatrix.dot(res[col:col + nvarPhase])
        footPos[moving] = footPos[moving] + [pos]
        if len(footPos[moving]) <2:
            footPos[moving] = footPos[moving] + [pos]
        com = zeros(3);
        com[2] = res[col + nvarPhase - 1]
        if len(footPos[fixed]) > 0:
            footPos[fixed]  = footPos[fixed] + [footPos[fixed][-1]] #duplicate fixed foot
            com[:2] = footPos[fixed][-1][:2]
        coms += [com]
        allFeetPos += [pos]
        col += nvarPhase
    return coms, footPos, allFeetPos
    
def plotBezier(bez, ax, color = "b", label = None, linewidth = 2.0, D3 = True, mx = None):
        step = 1000.
        if mx is None:
                mx = bez.max()
        points1 =  np.array([(bez(i/step*mx)[0],bez(i/step*mx)[1],bez(i/step*mx)[2]) for i in range(int(step))])
        x = points1[:,0]
        y = points1[:,1]
        if(D3):                
                z = points1[:,2]
                z = [el[0] for el in z]
                ax.plot(x.tolist(),y.tolist(),z,color)
        else:
                ax.plot(x.tolist(),y.tolist() ,color,linewidth=linewidth, label=label)
        return points1
    
def plotPoints(ax, wps, color = "b", D3 = True, linewidth=2):
    x = array(wps)[:,0]
    y = array(wps)[:,1]
    if(D3):                
            z = array(wps)[:,2]
            ax.scatter(x, y, z, c=color, marker='o', linewidth = 5) 
    else:
            ax.scatter(x,y,color=color, linewidth = linewidth)  
   
from plot_plytopes import plot_polytope_H_rep
   
def plotConstraints(ax, pb, allfeetpos, coms):
    for i, phase in enumerate(pb["phaseData"][:]):
        if i <1 :
            continue
        fixed =   phase["fixed"]  
        moving = phase["moving"]   
        oldK, oldk = pb["phaseData"][i-1]["K"][0][fixed]
        oldK = oldK.copy()
        oldk = oldk.copy()
        oldk += oldK.dot(allfeetpos[i-1])
        K, k = phase["K"][0][moving]  
        K = K.copy()
        k = k.copy()
        pos =  allfeetpos[i]
        com = coms[i]
        relK, relk = pb["phaseData"][i-1]["relativeK"][0]
        relK = relK.copy()
        relk = relk.copy()
        relk += relK.dot(allfeetpos[i-1])
        
        k = k + K.dot(pos)
        resK = vstack([oldK,K])
        resk = concatenate([oldk, k]).reshape((-1,)) 
        if True:
            try :                
                plot_polytope_H_rep(relK,relk.reshape((-1,1)), ax = ax)
            except: 
                print "qhullfailed"
    
        
def plotQPRes(pb, res, linewidth=2, ax = None, plot_constraints = False, show = True):
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        
    
    plotPoints(ax, coms, color = "b")
    plotPoints(ax, footpos[RF], color = "r")
    plotPoints(ax, footpos[LF], color = "g")
    
    cx = [c[0] for c in coms]
    cy = [c[1] for c in coms]
    cz = [c[2] for c in coms]
    ax.plot(cx, cy, cz)
    px = [c[0] for c in allfeetpos]
    py = [c[1] for c in allfeetpos]
    pz = [c[2] for c in allfeetpos]
    ax.plot(px, py, pz)
    
    
    if plot_constraints:
        plotConstraints(ax, pb, allfeetpos, coms)
        
    if show:
        plt.show()
       

####################### BENCH ###################"

from time import  clock
def timMs(t1, t2):
    return (t2-t1) * 1000.
    
####################### MAIN ###################"

if __name__ == '__main__':
    from mpcroc.planner_scenarios.flat_ground import genFlatGroundProblem
    pb = genFlatGroundProblem([0.,0.05,0.],[0.,-0.05,0.], 10)
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    
    print 'A.shape', A.shape
    
    C, c = least_square_cost_function (pb, A.shape[1], comTarget = array([30.,10.,0.2]))
    t2 = clock()
    res = qp.solve_least_square(C,c,A,b,E,e)
    t3 = clock()
    
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    plotQPRes(pb, res, plot_constraints = False)
    #~ plotQPRes(pb, res, plot_constraints = True)
    
        
