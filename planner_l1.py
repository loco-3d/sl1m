

############### Problem definition #############

from sl1m.problem_definition import *
           
        
#### constraint specification ####  

DEFAULT_NUM_VARS = 4

#extra variables when multiple surface: a0, variable for inequality, a1 slack for equality constraint, then a2 = |a1|
#so vars are x, y, z, zcom, a0, a1
#extra constraints are a1 - a0 <= 0 ; -a1 - a0 <=  0;  a0 >= 0, a0 >=0 is implied by the first 2
NUM_SLACK_PER_SURFACE = 2
NUM_INEQ_SLACK_PER_SURFACE = 2 # 

def numVariablesForPhase(phase):
    ret = DEFAULT_NUM_VARS
    numSurfaces =len(phase["S"])
    if numSurfaces > 1:
        ret += numSurfaces * NUM_SLACK_PER_SURFACE
    return ret
    
def numIneqForPhase(phase):
    #surface -one because last is equality
    ret =  sum([k.shape[0]  for (_, k) in phase["K"][0] ]) + phase["relativeK"][0][1].shape[0] + sum([S[1].shape[0]-1 for S in  phase["S"]])
    numSurfaces =len(phase["S"])
    if numSurfaces > 1:
        ret += numSurfaces * NUM_INEQ_SLACK_PER_SURFACE
    return ret
    
def getTotalNumVariablesAndIneqConstraints(pb):
    nphases = pb["nphases"]
    cols = sum([numVariablesForPhase(phase) for phase in  pb["phaseData"]])
    rows = sum( [numIneqForPhase(phase) for  phase in pb["phaseData"] ])
    return rows, cols
    
def numEqForPhase(phase):
    return len(phase["S"])
    
def getTotalNumEqualityConstraints(pb):
    return sum( [numEqForPhase(phase) for  phase in pb["phaseData"] ])
  
  
    
#### Selection matrices ####  

# COM position given by the first two parameters (x, y of foot) and the fourth one (height)
def _COMExpressionMatrix():
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:2,:2] = identity(2)
    ret[2, -1] = 1
    return ret
    
    
# Foot position given by the first three parameters (x, y, z of foot)
def _footExpressionMatrix():
    return identity(3)
    
def _footExpressionMatrixXY():
    return identity(2)
    
def _footExpressionMatrixz():
    ret =  zeros((1,3))
    ret[0,-1] =  1
    return ret
    
# z position of COM
def _COMzmatrix():
    ret = zeros((1, DEFAULT_NUM_VARS))
    ret[-1,-1:] = [1]
    return ret
 
COMExpressionMatrix = _COMExpressionMatrix()
footExpressionMatrix = _footExpressionMatrix()
footExpressionMatrixz = _footExpressionMatrixz()
footExpressionMatrixXY = _footExpressionMatrixXY()
COMzmatrix = _COMzmatrix()
    
# com is always above foot, so we lony care about the z value
def FixedFootCOMKinConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol):
    fixed = phaseDataT["fixed"]
    pos   = pb["p0"][fixed]
    K, k = phaseDataT["K"][0][fixed]
    idRow =  K.shape[0]
    comOff = COMzmatrix.shape[1]
    
    A[:idRow, startCol:startCol+comOff] = K[:,-1:].dot(COMzmatrix)
    b[:idRow                 ] = k + K[:,-1:].dot(pos[-1:])
    return idRow
    
def FixedFootCOMKinConstraintVarPhase(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["K"][0][fixed]
    idRow = startRow + K.shape[0]
    offCOm  =  COMzmatrix.shape[1]
    offFoot = footExpressionMatrixz.shape[1]
    A[startRow:idRow, startCol:startCol+offCOm]         =  K[:,-1:].dot(COMzmatrix)
    A[startRow:idRow, previousStartCol:previousStartCol+offFoot] = -K[:,-1:].dot(footExpressionMatrixz)
    b[startRow:idRow                 ] = k 
    
    return idRow
    
    
def FixedFootCOMConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    if startRow == 0:
        if pb["p0"] is None:
            fixed = phaseDataT["fixed"] #0 constraints
            K, k = phaseDataT["K"][0][fixed]
            return K.shape[0]
        else:
            return FixedFootCOMKinConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol)
    else:
        return FixedFootCOMKinConstraintVarPhase (pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow)
        
        
#com position at this phase is contrained by moving foot position
def MovingFootCOMConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"] 
    pos = pb["p0"][fixed]
    moving = phaseDataT["moving"]
    K, k = phaseDataT["K"][0][moving]
    idRow = startRow + K.shape[0]
    offCOmZ  = COMzmatrix.shape[1]
    offFoot = footExpressionMatrix.shape[1]
    A[startRow:idRow, startCol:startCol+offCOmZ]  =  K[:,-1:].dot(COMzmatrix)
    A[startRow:idRow, startCol:startCol+offFoot] -=  K.dot(footExpressionMatrix)
    b[startRow:idRow                 ] = k - K.dot(pos)
    return idRow


def MovingFootCOMConstraintVarPhase(phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):    
    moving = phaseDataT["moving"]
    K, k = phaseDataT["K"][0][moving]
    idRow = startRow + K.shape[0]
    offCOmZ  = COMzmatrix.shape[1]
    offFootXY = footExpressionMatrixXY.shape[1]
    offFoot   = footExpressionMatrix.shape[1]
    A[startRow:idRow, previousStartCol:previousStartCol+offFootXY]    =  K[:,:2].dot(footExpressionMatrixXY)
    A[startRow:idRow, startCol:startCol+offCOmZ]         =  K[:,-1:].dot(COMzmatrix)
    A[startRow:idRow, startCol:startCol+offFoot] = -K.dot(footExpressionMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
    
def MovingFootCOMConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow, first):
    if first:
        if pb["p0"] is None:
            moving = phaseDataT["moving"] #0 constraints
            K, k = phaseDataT["K"][0][moving]
            return startRow + K.shape[0]
        else:
            return MovingFootCOMConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol, startRow)
    else:
        return MovingFootCOMConstraintVarPhase (phaseDataT, A, b, previousStartCol, startCol, endCol, startRow)
    
    
def FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    fixed = phaseDataT["fixed"]
    K, k = phaseDataT["relativeK"][0]
    footMatrix = currentFootExpressionMatrix(phaseDataT, startCol, endCol)
    pos   = pb["p0"][fixed]
    idRow = K.shape[0]
    offFoot = footMatrix.shape[1]
    A[:idRow, startCol:startCol+offFoot] = K.dot(footMatrix)
    b[:idRow                 ] = k + K.dot(pos)
    return idRow + startRow
    
def FixedFootRelativeDistanceConstraintVarPhase(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    K, k = phaseDataT["relativeK"][0]
    idRow = startRow + K.shape[0]
    offFoot = footExpressionMatrix.shape[1]
    A[startRow:idRow, startCol:startCol+offFoot] = K.dot(footExpressionMatrix)
    A[startRow:idRow, previousStartCol:previousStartCol+offFoot] = -K.dot(footExpressionMatrix)
    b[startRow:idRow                 ] = k 
    return idRow
    
    
def FixedFootConstraintRelativeDistance(pb, phaseDataT, A, b, previousCol, startCol, endCol, startRow, first):
    if first:
        if True:
            K, k = phaseDataT["relativeK"][0] #0 constraints
            return startRow + K.shape[0]
        else:
            return FixedFootRelativeDistanceConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol, startRow)
    else:
        return FixedFootRelativeDistanceConstraintVarPhase (pb, phaseDataT, A, b, previousCol, startCol, endCol, startRow)
        
    
def SurfaceConstraint(phaseDataT, A, b, startCol, endCol, startRow):
    sRow = startRow
    nSurfaces = len(phaseDataT["S"])
    idS = DEFAULT_NUM_VARS
    for (S,s) in phaseDataT["S"]:   
        idRow = sRow + S.shape[0]-1
        A[sRow:idRow, startCol:startCol+footExpressionMatrix.shape[1]] = S[:-1,:].dot(footExpressionMatrix)
        b[sRow:idRow                 ] = s[:-1]
        if nSurfaces >1:
            A[sRow:idRow, startCol+idS] = -ones(idRow-sRow) 
        sRow = idRow
        idS += NUM_SLACK_PER_SURFACE
    return idRow
    
def SlackPositivityConstraint(phaseDataT, A, b, startCol, endCol, startRow):    
    numSurfaces =len(phaseDataT["S"])
    idRow = startRow
    if numSurfaces > 1:
        varOffset = DEFAULT_NUM_VARS #foot and com positions
        for row in range(startRow, startRow + numSurfaces * NUM_INEQ_SLACK_PER_SURFACE, NUM_INEQ_SLACK_PER_SURFACE):
            i = row - startRow
            col = startCol+varOffset + i / NUM_INEQ_SLACK_PER_SURFACE * NUM_SLACK_PER_SURFACE
            A[row  , col:col+2] = [-1,  1];  #  - a0 + a1  <= 0
            A[row+1, col:col+2] = [-1, -1];  # -a0 - a1 <= 0
            idRow = row + NUM_INEQ_SLACK_PER_SURFACE       
    return idRow
    
def EqualityConstraint(phaseDataT, E, e, startCol, endCol, startRowEq):    
    sRow = startRowEq
    nSurfaces = len(phaseDataT["S"])
    if nSurfaces == 1:
        E[sRow, startCol:startCol+footExpressionMatrix.shape[1]] = phaseDataT["S"][0][0][-1]
        e[sRow] = phaseDataT["S"][0][1][-1]
        return sRow+1
    else:
        idS = DEFAULT_NUM_VARS
        for (S,s) in phaseDataT["S"]:   
            E[sRow, startCol:startCol+footExpressionMatrix.shape[1]] = S[-1,:].dot(footExpressionMatrix)
            #~ E[sRow, startCol+idS+1] = -1 # E x  + a1 = e
            E[sRow, startCol+idS+1] = -1 # E x  + a1 = e
            e[sRow] = s[-1]
            idS += NUM_SLACK_PER_SURFACE
            sRow += 1
        return sRow 
    
def convertProblemToLp(pb, convertSurfaces = True):    
    if convertSurfaces:
        replace_surfaces_with_ineq_in_problem(pb)
    #define first problem
    #A u <= b
    nIneq, nvars  = getTotalNumVariablesAndIneqConstraints(pb)
    A = zeros((nIneq, nvars)); b = zeros(nIneq)
    #E u = b
    nEq = getTotalNumEqualityConstraints(pb)
    E = zeros((nEq, nvars)); e = zeros(nEq)
    
    startRow = 0;
    startRowEq = 0;
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    #~ fixedFootMatrix = None;
    
    for i, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        startRow = FixedFootCOMConstraint (pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow) 
        startRow = FixedFootConstraintRelativeDistance (pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow, first = i == 0) 
        startRow = MovingFootCOMConstraint(pb,phaseDataT, A, b, previousStartCol, startCol, endCol, startRow, first = i == 0)
        startRow = SurfaceConstraint(phaseDataT, A, b, startCol, endCol, startRow)
        startRow = SlackPositivityConstraint(phaseDataT, A, b, startCol, endCol, startRow)
        startRowEq = EqualityConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)
        previousStartCol = startCol
        startCol   = endCol 
    
    assert startRowEq == E.shape[0]
    assert startRow   == A.shape[0]
    
    A,b = normalize([A,b])
    E,e = normalize([E,e])
    return (A,b,E,e)
        
def slackSelectionMatrix(pb):
    nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    c = zeros(nvars)
    cIdx = 0
    for i, phase in enumerate(pb["phaseData"]):
        phaseVars = numVariablesForPhase(phase)
        nslacks = phaseVars - DEFAULT_NUM_VARS
        startIdx = cIdx + DEFAULT_NUM_VARS
        for i in range (0,nslacks,NUM_SLACK_PER_SURFACE):
            c[startIdx + i] = 1
        cIdx += phaseVars
    assert cIdx == nvars
    return c
    
def num_non_zeros(pb, res):
    nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    indices = []
    cIdx = 0
    wrongsurfaces = []
    for i, phase in enumerate(pb["phaseData"]):  
        numSurfaces = len(phase["S"])
        phaseVars = numVariablesForPhase(phase)
        if numSurfaces > 1:
            startIdx = cIdx + DEFAULT_NUM_VARS
            betas = [res[startIdx+j] for j in range(0,numSurfaces*2,2) ]
            if array(betas).min() > 0.01:
                #~ print "wrong ", i, array(betas).min()
                indices += [i]
                sorted_surfaces = np.argsort(betas)
                #~ print "sorted_surfaces ",sorted_surfaces
                wrongsurfaces += [[[phase["S"][idx]] for idx in sorted_surfaces]  ]
                #~ print "lens ", len([[phase["S"][idx]] for idx in sorted_surfaces]  )
        cIdx += phaseVars
    return indices, wrongsurfaces
    
def isSparsityFixed(pb, res):
    indices, wrongsurfaces = num_non_zeros(pb, res)
    return len(indices) == 0
   
def genOneComb(pb,indices, surfaces, res):
    pb1 = pb.copy()
    for i, idx in enumerate(indices):
        pb1["phaseData"][idx]["S"] = surfaces[i][0]
    res += [pb1]
   
import itertools
import copy
   
def genCombinatorialRec(pb, indices, wrongsurfaces, res):
    lenss  = [len(surfs) for surfs in wrongsurfaces]
    all_indexes = [[el for el in range(lens)]  for lens in [len(surfs) for surfs in wrongsurfaces]]
    combs = [el for el in itertools.product(*all_indexes)]
    for comb in combs:
        pb1 = copy.deepcopy(pb)
        for i, idx in enumerate(indices):
            pb1["phaseData"][idx]["S"] = wrongsurfaces[i][comb[i]]
        res += [[pb1, comb, indices]]
        
    
    
def generateAllFixedScenariosWithFixedSparsity(pb, res):
    indices, wrongsurfaces = num_non_zeros(pb, res)
    all_len = [len(surfs) for surfs in wrongsurfaces]
    comb = 1
    for el in all_len:
        comb *= el  
    res = []
    if comb >1000:
        print "problem probably too big ", comb
    else:
        genCombinatorialRec(pb, indices, wrongsurfaces, res)
    return res
    
    
    
def bestSelectedSurfaces(pb, res):
    surfaces = []
    indices  = []
    cIdx = 0
    nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    for i, phase in enumerate(pb["phaseData"]):  
        numSurfaces = len(phase["S"])
        phaseVars = numVariablesForPhase(phase)
        if numSurfaces == 1:
            surfaces = surfaces + [phase["S"][0]]
            indices = indices + [0]
        else:
            startIdx = cIdx + DEFAULT_NUM_VARS
            betas = [res[startIdx+j] for j in range(0,numSurfaces*2,2) ]
            assert betas >= -0.00000001
            bestIdx = betas.index(array(betas).min())
            surfaces = surfaces + [phase["S"][bestIdx]]
            
            indices = indices + [bestIdx]
        cIdx += phaseVars
    assert cIdx == nvars
    return surfaces, indices
    
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
        pos = footExpressionMatrix.dot(res[col:col + footExpressionMatrix.shape[1]])
        footPos[moving] = footPos[moving] + [pos]
        com = zeros(3);
        com[2] = res[col + 3]
        if len(footPos[moving]) <2:
            footPos[moving] = footPos[moving] + [pos]
        if len(footPos[fixed]) > 0:
            footPos[fixed]  = footPos[fixed] + [footPos[fixed][-1]] #duplicate fixed foot
            com[:2] = footPos[fixed][-1][:2]
        coms += [com]
        allFeetPos += [pos]
        col += numVariablesForPhase(phaseDataT)
    return coms, footPos, allFeetPos
    
    
def plotPoints(ax, wps, color = "b", D3 = True, linewidth=2):
    x = array(wps)[:,0]
    y = array(wps)[:,1]
    if(D3):                
            z = array(wps)[:,2]
            ax.scatter(x, y, z, color=color, marker='o', linewidth = 5) 
    else:
            ax.scatter(x,y,color=color, linewidth = linewidth)  
   
from tools.plot_plytopes import plot_polytope_H_rep
   
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
        #~ if i %2 ==0:
        #~ if i %2 ==1:
            try :                
                #~ plot_polytope_H_rep(resK,resk.reshape((-1,1)), ax = ax)
                plot_polytope_H_rep(relK,relk.reshape((-1,1)), ax = ax)
                #~ plot_polytope_H_rep(K,k.reshape((-1,1)), ax = ax)
            except: 
                print "qhullfailed"
    
        
def plotQPRes(pb, res, linewidth=2, ax = None, plot_constraints = False, show = True):
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)
        
    ax.set_autoscale_on(False)
    ax.view_init(elev=8.776933438381377, azim=-99.32358055821186)
    
    #~ plotPoints(ax, coms, color = "b")
    plotPoints(ax, footpos[RF], color = "r")
    plotPoints(ax, footpos[LF], color = "g")
    
    cx = [c[0] for c in coms]
    cy = [c[1] for c in coms]
    cz = [c[2] for c in coms]
    #~ ax.plot(cx, cy, cz)
    px = [c[0] for c in allfeetpos]
    py = [c[1] for c in allfeetpos]
    pz = [c[2] for c in allfeetpos]
    ax.plot(px, py, pz, color="black")
        
    if show:
        plt.ion()
        plt.show()
       
    
####################### MAIN ###################"

if __name__ == '__main__':
    from sl1m.stand_alone_scenarios.complex import gen_stair_pb,  draw_scene
    pb = gen_stair_pb()    
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    
    A,b = normalize([A,b])
    C = identity(A.shape[1]) * 0.00001
    c =  slackSelectionMatrix(pb) * 100.
    t2 = clock()
    res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    t3 = clock()
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    surfaces, indices = bestSelectedSurfaces(pb, res)
    for i, phase in enumerate(pb["phaseData"]):  
        phase["S"] = [surfaces[i]]
        
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb, False)
    
    
    C = identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    t3 = clock()
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False)
        
    
    
    
    
        
