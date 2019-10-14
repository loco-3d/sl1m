import numpy as np


from constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *

from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray
import qp



############### Problem definition #############


LF = 0
RF = 1

def normalize(Ab):
    A = Ab[0]
    b = Ab[1]
    Ares = zeros(A.shape)
    bres = zeros(b.shape)
    for i in range(A.shape[0]):
        n = norm(A[i,:])
        if n <= 0.000001:
            n = 1.
        Ares[i,:] = A[i,:] / n; bres[i] = b[i] / n
    #~ print Ares
    return Ares, bres

def genKinematicConstraints(normals = [z, z], min_height = None):
    res = [None, None]
    trLF = default_transform_from_pos_normal(zero3, normals[LF])
    trRF = default_transform_from_pos_normal(zero3, normals[RF])
    KLF = left_foot_hrp2_constraints  (trLF)
    KRF = right_foot_hrp2_constraints (trRF)
    #~ KLF = (zeros((3,3)), zeros(3))
    #~ KRF = (zeros((3,3)), zeros(3))
    if min_height is None:
        res [LF] = KLF
        res [RF] = KRF
    else:
        res [LF] = addHeightConstraint(KLF[0], KLF[1], min_height)
        res [RF] = addHeightConstraint(KRF[0], KRF[1], min_height)
    #~ res [LF] = normalize(res [LF])
    #~ res [RF] = normalize(res [RF])
    return res

def genFootRelativeConstraints(normals = [z, z]):
    res = [None, None]
    trLF = default_transform_from_pos_normal(zero3, normals[LF])
    trRF = default_transform_from_pos_normal(zero3, normals[RF])
    KLF = right_foot_in_lf_frame_hrp2_constraints  (trLF)
    KRF = left_foot_in_rf_frame_hrp2_constraints (trRF)
    
    #~ KLF = (zeros((3,3)), zeros(3))
    #~ KRF = (zeros((3,3)), zeros(3))
    res [LF] = KLF #constraints of right foot in lf frame. Same idea as COM in lf frame
    res [RF] = KRF
    return res
    
def copyKin(kC):
    return [(Kk[0].copy(), Kk[1].copy()) for Kk in kC]

# generate problem for flat ground scenarios, from init pos for feet p0LF and p0RF
# this method also defines the template definition for defining contact plan
def genFlatGroundProblem(p0LF, p0RF, nphases = 5):
    #~ for i in range(nphases)
    #~ kinematicConstraints = genKinematicConstraints(min_height = 0.4)
    kinematicConstraints = genKinematicConstraints()
    relativeConstraints = genFootRelativeConstraints()
    surfaces = [None, None];
    p0 = [[0., -0.06,0.], [0., 0.06,0.]];
    #~ surfaces[LF] = array([[-50.,0.04,0.], [50.,0.04,0.], [50.,40.,0.], [-50.,40.,0.] ]).T  # for the moment arbitrary surfaces to separate foot
    surfaces[LF] = array([[-50.,-40,0.], [50.,-40.,0.], [50.,-0.04,0.], [-50.,-0.04,0.] ]).T
    surfaces[RF] = array([[-50.,-40,0.], [50.,-40.,0.], [50.,-0.04,0.], [-50.,-0.04,0.] ]).T
    p0 [LF]  = p0LF;
    p0 [RF]  = p0RF;
    
    #~ res = { "p0" : p0, "nphases": nphases}
    res = { "p0" : None, "nphases": nphases}
    
    #TODO in non planar cases, K must be rotated
    #relative k is the constraint with respect to the previous frame
    phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : copyKin(kinematicConstraints), "S" : surfaces[i%2], "relativeK" : relativeConstraints[(i) % 2] } for i in range(nphases)]
    res ["phaseData"] = phaseData
    return res
    

def getTotalNumVariablesAndIneqConstraints(pb):
    cols    = sum([phase["S"].shape[1] + 1 for  phase in pb["phaseData"] ])
    rows    = sum( [sum([k.shape[0]  for (_, k) in phase["K"] ]) + phase["relativeK"][1].shape[0] for  phase in pb["phaseData"] ])
    rowsOld = sum( [sum([k.shape[0]  for (_, k) in phase["K"] ])  for  phase in pb["phaseData"] ])
    rows += cols #positivity constraint
    return rows, cols
    
def getTotalNumEqualityConstraints(pb):
    return pb["nphases"]
  
def numVariablesForPhase(phase):
    return phase["S"].shape[1] + 1
  
#global selection matrix for retrieving x and y values of COM
#~ select = array([[1., 0., 0.],[0., 1., 0.]])

    
# COM position is expressed as convex combination from contact surfaces for x and y,
# while z is directly given by a variable
def currentCOMExpressionMatrix(phaseDataT, startCol, endCol):
    lenVar = endCol - startCol
    ret = zeros((3, lenVar))
    #surfaces
    ret[:2,:-1] = phaseDataT["S"][:2,:]
    ret[-1, -1] = 1.  #retrieve com z 
    
    return ret
    
# get current alpha values
def currentAlphaExpressionMatrix(phaseDataT, startCol, endCol):
    lenVar = endCol - startCol
    ret = zeros((lenVar-1, lenVar))
    for i in range(lenVar-1):
        ret[i,i] = 1
    #~ ret[-1, -1] = 1.  #retrieve com z 
    
    return ret
    
    
# Foot position matrix
# used to retrieve foot position in the previous phase essentially
def currentFootExpressionMatrix(phaseDataT, startCol, endCol):
    lenVar = endCol - startCol
    ret = zeros((3, lenVar))
    ret[:,:-1] = phaseDataT["S"][:,:]    
    return ret
    
# z position relative to support foot
def currentCOMzmatrix(phaseDataT, startCol, endCol):
    ret = zeros((1, endCol - startCol))
    ret[-1,:-1] = -phaseDataT["S"][-1,:]    
    ret[-1,- 1] = 1.
    return ret
    
    
def FixedFootKinConstraintInitPhase(pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol):
    fixed = phaseDataT["fixed"]
    pos   = pb["p0"][fixed]
    K, k = phaseDataT["K"][fixed]
    idRow = K.shape[0]
    A[:idRow, startCol:endCol] = K.dot(comMatrix)
    b[:idRow                 ] = k + K.dot(pos)
    return idRow
    
def FixedFootKinConstraintVarPhase(pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["K"][fixed]
    idRow = startRow + K.shape[0]
    A[startRow:idRow, startCol:endCol] = K.dot(comMatrix)
    #prev variable size
    prevSize  = fixedFootMatrix.shape[1]
    A[startRow:idRow, (startCol-prevSize):startCol] = -K.dot(fixedFootMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
    
    
def FixedFootConstraint(pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow):
    if startRow == 0:
        if pb["p0"] is None:
            fixed = phaseDataT["fixed"] #0 constraints
            K, k = phaseDataT["K"][fixed]
            return K.shape[0]
        else:
            return FixedFootKinConstraintInitPhase(pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol)
    else:
        return FixedFootKinConstraintVarPhase (pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow)
        
def FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["relativeK"]
    footMatrix = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
    pos   = pb["p0"][fixed]
    idRow = K.shape[0]
    A[:idRow, startCol:endCol] = K.dot(footMatrix)
    b[:idRow                 ] = k + K.dot(pos)
    return idRow + startRow
    
def FixedFootRelativeDistanceConstraintVarPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow):
    K, k = phaseDataT["relativeK"]
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
        #~ if pb["p0"] is None:
        if True:
            K, k = phaseDataT["relativeK"] #0 constraints
            return startRow + K.shape[0]
        else:
            return FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow)
    else:
        return FixedFootRelativeDistanceConstraintVarPhase (pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow)
        
        
# com is always above foot, so we lony care about the z value
def MovingFootConstraint(phaseDataT, A, b, startCol, endCol, startRow):
    #get z selection matrix
    zCOMMatrix = currentCOMzmatrix(phaseDataT, startCol, endCol)
    
    moving = phaseDataT["moving"]
    K, k = phaseDataT["K"][moving]
    idRow = startRow + K.shape[0]
    A[startRow:idRow, startCol:endCol] = K[:,-1:].dot(zCOMMatrix)
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

def convertProblemToLp(pb):
    #define first problem
    #A u <= b
    nIneq, nvars  = getTotalNumVariablesAndIneqConstraints(pb)
    A = zeros((nIneq, nvars)); b = zeros(nIneq)
    #E u = b
    nEq = getTotalNumEqualityConstraints(pb)
    E = zeros((nEq, nvars)); e = zeros(nEq)
    
    startRow = 0;
    startCol = 0;
    endCol   = 0;
    fixedFootMatrix = None;
    
    startRowEq = 0;
    for i, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        comMatrix = currentCOMExpressionMatrix(phaseDataT, startCol, endCol)     
        startRow = FixedFootConstraint (pb, phaseDataT, comMatrix, fixedFootMatrix, A, b, startCol, endCol, startRow) 
        startRow = FixedFootConstraintRelativeDistance (pb, phaseDataT, fixedFootMatrix, A, b, startCol, endCol, startRow, first = i == 0) 
        startRow = MovingFootConstraint(phaseDataT, A, b, startCol, endCol, startRow)
        
        #equality
        startRowEq = ConvexConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)
        
        # for next phase
        fixedFootMatrix = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
        startCol   = endCol 
        
    #add positivity constraint
    A[startRow:, : ] = -identity(nvars)
    b[-nvars:      ] = zeros(nvars)
    
    #TODO quick hack for constraining COM height : remove me
    #~ sc = 0
    #~ for i, phaseDataT in enumerate(pb["phaseData"]): 
        #~ sc += numVariablesForPhase(phaseDataT)
        #~ b[startRow + sc-1] = -0.4
    
    startRow += nvars
    
    assert startRowEq == E.shape[0]
    assert startRow   == A.shape[0]
    
    return (A,b,E,e)
    
def L1NormConstraintVarPhase(pb, phaseDataT, prevAlphaMatrix, currentAlphaMatrix, A, b, startCol, endCol, startRow, l1varCol):
    assert prevAlphaMatrix.shape[0] == currentAlphaMatrix.shape[0]
    idRow = startRow + currentAlphaMatrix.shape[0]
    print "startRow", startRow
    print "startCol", startCol
    print "endCol", endCol
    print "idRow", idRow
    print "l1varCol", l1varCol
    A[startRow:idRow, startCol:endCol] = currentAlphaMatrix
    #prev variable size
    prevSize  = prevAlphaMatrix.shape[1]
    A[startRow:idRow, (startCol-prevSize):startCol] = -prevAlphaMatrix
    print "startCol prev", (startCol-prevSize)
    idRow2 = idRow + currentAlphaMatrix.shape[0]
    print "l1varCol2", l1varCol+idRow2-idRow
    A[idRow:idRow2, :] = -A[startRow:idRow,:]
    A[startRow:idRow, l1varCol:l1varCol+idRow2-idRow] = -identity(idRow2-idRow)
    A[idRow:idRow2  , l1varCol:l1varCol+idRow2-idRow] = -identity(idRow2-idRow)
    print "constraint ", A[startRow:idRow2, l1varCol:l1varCol+idRow2-idRow].shape
    print "constraint A ", A[startRow:idRow, (startCol-prevSize):startCol]
    print "constraint A ", A[startRow:idRow, startCol:endCol]
    print "constraint B ", A[idRow:idRow2, (startCol-prevSize):startCol]
    print "constraint B ", A[idRow:idRow2, startCol:endCol]
    return idRow2, l1varCol + idRow2 - idRow
    
def add_l1_norm_minimization(pb, Ao, bo, Eo, eo, Co, co, weight =1., reg_term = 0.0001):
    nphases     = pb["nphases"] -1 
    n_add_vars  = sum([phase["S"].shape[1] for  phase in pb["phaseData"][2:] ])
    
    A = hstack([Ao,zeros((Ao.shape[0], n_add_vars))])
    E = hstack([Eo,zeros((Eo.shape[0], n_add_vars))])
    c = concatenate([co, -ones(n_add_vars) * weight])
    C = identity((Co.shape[0]+n_add_vars)) * reg_term
    C [:Co.shape[0],:Co.shape[0]] = Co
    
    extra_constraints = n_add_vars * 2
    A = vstack([A,zeros((extra_constraints, A.shape[1]))])
    print "new shape ", A.shape
    b = concatenate([bo, zeros(extra_constraints)])
    
    startRow = Ao.shape[0];
    startCol = 0;
    endCol   = 0;
    startColNewVar  = Ao.shape[1];
    print "nstartColNewVar ", startColNewVar
    
    prevAlphaMatrix = None
    
    for i, phaseDataT in enumerate(pb["phaseData"]):
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        currentAlphaMatrix = currentAlphaExpressionMatrix(phaseDataT, startCol, endCol)
        if i > 1:
            startRow, startColNewVar = L1NormConstraintVarPhase(pb, phaseDataT, prevAlphaMatrix, currentAlphaMatrix, A, b, startCol, endCol, startRow, startColNewVar)        
        # for next phase
        prevAlphaMatrix = currentAlphaMatrix
        startCol   = endCol 
    
    return A, b, E, e, C, c
    
def least_square_cost_function(pb, nvars, comTarget, weight =1., reg_term = 0.0001):
    C = identity(nvars) * reg_term 
    lastPhase = pb["phaseData"][-1]
    nvarPhase = numVariablesForPhase(lastPhase)
    lastCom = currentCOMExpressionMatrix(lastPhase, 0, nvarPhase)
    C[-3:,nvars-nvarPhase:] = lastCom
    return C * weight, concatenate([zeros(nvars-3),comTarget])
    
###########################" PLOTTING ################"
    
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def retrieve_points_from_res(pb, res):
    coms = []
    if pb["p0"] is not None:
    #~ if False:
        footPos = [[pb["p0"][LF]],[pb["p0"][RF]]]
        allFeetPos = [footPos[0][0], footPos[1][0]]
        coms       = [footPos[0][0], footPos[1][0]]
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
            print "inf 2" , i
            footPos[moving] = footPos[moving] + [pos]
        if len(footPos[fixed]) > 0:
            print "sup 1" , i
            footPos[fixed]  = footPos[fixed] + [footPos[fixed][-1]] #duplicate fixed foot
        print "phase ", i , " ", len(footPos[0]), " ", len(footPos[1])
        print "phase ", i , " ", len(allFeetPos)
        #~ else:
            #~ footPos[fixed] = [None]
        com = pos.copy()
        com[2] = res[col + nvarPhase - 1]
        coms += [com]
        allFeetPos += [pos]
        col += nvarPhase
    return coms, footPos, allFeetPos
    
def plotBezier(bez, ax, color = "b", label = None, linewidth = 2.0, D3 = True, mx = None):
        step = 1000.
        if mx is None:
                mx = bez.max()
        #~ points1 =  np.array([(bez(i/step*mx)[0][0],bez(i/step*mx)[1][0],bez(i/step*mx)[2][0]) for i in range(int(step))])
        points1 =  np.array([(bez(i/step*mx)[0],bez(i/step*mx)[1],bez(i/step*mx)[2]) for i in range(int(step))])
        #~ print "? ", bez(bez.max())
        x = points1[:,0]
        y = points1[:,1]
        if(D3):                
                z = points1[:,2]
                z = [el[0] for el in z]
                #~ print "bezer"
                ax.plot(x.tolist(),y.tolist(),z,color)
                #~ plt.show()
        else:
                ax.plot(x.tolist(),y.tolist() ,color,linewidth=linewidth, label=label)
                #~ plt.show()
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
        #~ if i == 0:
        if i <1 :
            continue
        fixed =   phase["fixed"]  
        moving = phase["moving"]   
        oldK, oldk = pb["phaseData"][i-1]["K"][fixed]
        oldK = oldK.copy()
        oldk = oldk.copy()
        oldk += oldK.dot(allfeetpos[i-1])
        K, k = phase["K"][moving]  
        K = K.copy()
        k = k.copy()
        pos =  allfeetpos[i]
        com = coms[i]
        relK, relk = pb["phaseData"][i-1]["relativeK"]
        relK = relK.copy()
        relk = relk.copy()
        relk += relK.dot(allfeetpos[i-1])
        
        k = k + K.dot(pos)
        resK = vstack([oldK,K])
        resk = concatenate([oldk, k]).reshape((-1,)) 
        if True:
        #~ if i %2 ==0:
        #~ if i %2 ==1:
            try :                
                #~ plot_polytope_H_rep(resK,resk.reshape((-1,1)), ax = ax)
                plot_polytope_H_rep(relK,relk.reshape((-1,1)), ax = ax)
                #~ plot_polytope_H_rep(K,k.reshape((-1,1)), ax = ax)
            except: 
                print "qhullfailed"
        #~ oldpos = pos
    
        
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
    pb = genFlatGroundProblem([0.,0.05,0.],[0.,-0.05,0.], 10)
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    
    print 'A.shape', A.shape
    
    C, c = least_square_cost_function (pb, A.shape[1], comTarget = array([30.,10.,0.2]))
    t2 = clock()
    
    c = zeros(c.shape[0])
    #~ A, b, E, e, C, c = add_l1_norm_minimization(pb, A, b, E, e, C, c, weight =1., reg_term = 0.0001)
    
    #~ res = qp.solve_least_square(C,c,A,b,E,e)
    res = qp.solve_lp(-c,A,b,E,e)["x"]
    t3 = clock()
    
    
    #~ add_l1_norm_minimization
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    plotQPRes(pb, res, plot_constraints = False)
    #~ plotQPRes(pb, res, plot_constraints = True)
    
        
