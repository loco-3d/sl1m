
############### Problem definition #############
from sl1m.problem_definition import *

#~ NUM VARIABLES

def getTotalNumVariablesAndIneqConstraints(pb):
    cols    = sum([phase["S"][0].shape[1] + 1 for  phase in pb["phaseData"] ])
    rows    = sum( [sum([k.shape[0]  for (_, k) in phase["K"][0] ]) + phase["relativeK"][0][1].shape[0] for  phase in pb["phaseData"][1:] ])
    rowsOld = sum( [sum([k.shape[0]  for (_, k) in phase["K"][0] ])  for  phase in pb["phaseData"] ])
    rows += cols #positivity constraint
    return rows, cols
    
def getTotalNumEqualityConstraints(pb):
    return pb["nphases"]
  
def numVariablesForPhase(phase):
    return phase["S"][0].shape[1] + 1
      
#~ FOOT AND COM EXPRESSION MATRICES
    
# z position relative to support foot
def comZExpr(phaseDataT, startCol, endCol):
    ret = zeros((3, endCol - startCol))
    ret[-1,phaseDataT["S"][0].shape[1]] = 1.
    return ret
      
#com is foot position at last phase + z from this phase
def comExpr(phaseDataT, fixedFootMatrix, startCol, endCol):
    comZMatrix = comZExpr(phaseDataT, startCol, endCol)
    ret = zeros((3,comZMatrix.shape[1]+fixedFootMatrix.shape[1]))
    ret[:2,:fixedFootMatrix.shape[1]] = fixedFootMatrix[:-1]
    ret[ 2, fixedFootMatrix.shape[1]:fixedFootMatrix.shape[1]+comZMatrix.shape[1]] = comZMatrix[-1,:]
    return ret
      
# Foot position matrix
# used to retrieve foot position in the previous phase essentially
def footMatrixExpr(phaseDataT, startCol, endCol):
    lenVar = endCol - startCol
    ret = zeros((3, lenVar))
    ret[:,:phaseDataT["S"][0].shape[1]] = phaseDataT["S"][0][:,:]    
    return ret
    
    
#~ CONSTRAINTS    
def footFrameConstraint(pb, phaseDataT, constrainedMatrix, frameMatrix, Kk,  A, b, endCol,startRow):    
    K, k = Kk[0], Kk[1]
    idRow = startRow + K.shape[0]
    A[startRow:idRow, endCol - constrainedMatrix.shape[1]:endCol] = K.dot(constrainedMatrix)
    A[startRow:idRow, endCol - frameMatrix.shape[1]:endCol]      -= K.dot(frameMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
    
#ensure alpha vector defines a convex combination 
def ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRow):
    ret = ones((1, endCol - startCol))
    ret[-1,-1] = 0
    idRow = startRow + 1
    E[startRow:idRow, startCol:endCol] = ret
    e[startRow:idRow                 ] = 1 
    return idRow

def addInitialPosPhase(pb):
    if pb["p0"] is not None:
        phaseNext = pb["phaseData"][0]
        phaseNextNext = pb["phaseData"][1]
        phase0 = {"S" : [array([pb["p0"][-1]]).T], "moving" : phaseNext["fixed"], "fixed" : phaseNext["moving"], "relativeK" : phaseNextNext["relativeK"], "K" : phaseNext["K"]  }
        pb["phaseData"] = [phase0] + pb["phaseData"]
        pb["p0"] = None
        pb["nphases"] =  pb["nphases"] +1

def convertProblemToLp(pb):
    addInitialPosPhase(pb)
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
    startRowEq = 0
    
    # initialisation. Add equality constraint on first position
    phaseDataT = pb["phaseData"][0]
    endCol = startCol + numVariablesForPhase(phaseDataT) 
    startRowEq = ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)    
    fixedFootMatrix = footMatrixExpr(phaseDataT, startCol, endCol)
    startCol   = endCol 
    
    for i, phaseDataT in enumerate(pb["phaseData"][1:]):   
        #inequality
        fixed  = phaseDataT["fixed"]
        moving = phaseDataT["moving"]
        nVarPhase = numVariablesForPhase(phaseDataT) 
        endCol = startCol + nVarPhase
        comMatrix =  comExpr(phaseDataT, fixedFootMatrix, startCol, endCol);  
        footMatrix = footMatrixExpr(phaseDataT, startCol, endCol)               
        fixedFootMatrix = hstack([fixedFootMatrix,zeros((fixedFootMatrix.shape[0], nVarPhase))]) # offset fixedfootmatrix
        
        startRow = footFrameConstraint(pb, phaseDataT, comMatrix , fixedFootMatrix, phaseDataT["K"][0][fixed ], A, b, endCol, startRow)
        startRow = footFrameConstraint(pb, phaseDataT, comMatrix , footMatrix     , phaseDataT["K"][0][moving], A, b, endCol, startRow)
        startRow = footFrameConstraint(pb, phaseDataT, footMatrix, fixedFootMatrix, phaseDataT["relativeK"][0], A, b, endCol, startRow)
        
        #equality
        startRowEq = ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)
        
        # for next phase
        fixedFootMatrix = footMatrix
        startCol   = endCol 
        
    #add positivity constraint
    A[startRow:startRow+nvars, : ] = -identity(nvars)
    b[-nvars:      ] = zeros(nvars)    
    
    startRow += nvars    
    assert startRowEq == E.shape[0]
    assert startRow   == A.shape[0] 
    
    return (A,b,E,e)
    
###########################" PLOTTING ################"
    
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def retrieve_points_from_res(pb, res):
    coms = []
    if pb["c0"] is not None:
        coms = [pb["c0"]]
                
    footPos = [[],[]]
    allFeetPos = []
    
    col = 0
    for i, phaseDataT in enumerate(pb["phaseData"]):  
        moving = phaseDataT["moving"]
        fixed = phaseDataT["fixed"]
        nvarPhase = numVariablesForPhase(phaseDataT)
        posMatrix = footMatrixExpr(phaseDataT, 0, nvarPhase)
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
        
def plotPoints(ax, wps, color = "b", D3 = True, linewidth=2):
    x = array(wps)[:,0]
    y = array(wps)[:,1]
    if(D3):                
            z = array(wps)[:,2]
            ax.scatter(x, y, z, color=color, marker='o', linewidth = 5) 
    else:
            ax.scatter(x,y,color=color, linewidth = linewidth)       
        
def plotQPRes(pb, res, linewidth=2, ax = None, show = True):
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)
    ax.view_init(elev=8.776933438381377, azim=-99.32358055821186)
    
    plotPoints(ax, coms, color = "b")
    plotPoints(ax, footpos[RF], color = "r")
    plotPoints(ax, footpos[LF], color = "g")
    
    cx = [c[0] for c in coms]
    cy = [c[1] for c in coms]
    cz = [c[2] for c in coms]
    ax.plot(cx, cy, cz, color="black")
    px = [c[0] for c in allfeetpos]
    py = [c[1] for c in allfeetpos]
    pz = [c[2] for c in allfeetpos]
    ax.plot(px, py, pz, color="black")
            
    if show:
        plt.draw()
        plt.show()
       

####################### MAIN ###################"

if __name__ == '__main__':
    from sl1m.stand_alone_scenarios.escaliers import gen_stair_pb, draw_scene
    pb = gen_stair_pb()
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)    
    C= identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    res = qp.solve_least_square(C,c,A,b,E,e)
    t3 = clock()
    
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax = ax)
    
        
