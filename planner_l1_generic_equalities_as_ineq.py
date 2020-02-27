

############### Problem definition #############

from sl1m.problem_definition import *
           
        
#### problem related global variables #### 

NUM_SLACK_PER_SURFACE = 1
NUM_INEQ_SLACK_PER_SURFACE = 1 # 
M = 100
M2= M*M
        
        
#### global variables, depending on number of effectors ####  
# You need to call initGlobals with the number of effectors before using sl1m generic

N_EFFECTORS                        = None
DEFAULT_NUM_VARS                   = None
COM_WEIGHTS_PER_EFFECTOR           = None
DEFAULT_NUM_EQUALITY_CONSTRAINTS   = None
DEFAULT_NUM_INEQUALITY_CONSTRAINTS = None

#### Selection matrices ####  
COM_XY_ExpressionMatrix  = None
COM_Z1_ExpressionMatrix  = None
COM_Z2_ExpressionMatrix  = None
COM1_ExpressionMatrix    = None
COM2_ExpressionMatrix    = None
foot_ExpressionMatrix    = None
foot_XY_ExpressionMatrix = None

#Slack variables startIndex
BETA_START  = None
GAMMA_START = None
W_START     = None
ALPHA_START = None



# breakdown of variables for each phase

# [ c_x't, c_y't, c_z_1't, c_z_2't, p_i't, b_i't, g_i't, w_i_t ,  a_i_l't]   
# [ 1    , 1    ,   1    ,    1   , 3*n_e, 2*n_e, 3*n_e,    n_e, [0, n_s]]   

# COM position
def _COM_XY_ExpressionMatrix():
    ret = zeros((2, DEFAULT_NUM_VARS))
    ret[:2,:2] = identity(2)
    return ret
    
def _COM_Z1_ExpressionMatrix():
    ret = zeros((1, DEFAULT_NUM_VARS))
    ret[0,2] = 1
    return ret
    
def _COM_Z2_ExpressionMatrix():
    ret = zeros((1, DEFAULT_NUM_VARS))
    ret[0,3] = 1
    return ret
    
def _COME1xpressionMatrix():
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:,:3] = identity(3)
    return ret
    
def _COME2xpressionMatrix():
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:2,:2] = identity(2)
    ret[2,  3] = 1
    return ret
    
def _footExpressionMatrix(footId):
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:, 4 + footId * 3, 4 + footId * 3 +3] = identity(3)
    return ret
    
def _foot_XY_ExpressionMatrix(footId):
    ret = zeros((2, DEFAULT_NUM_VARS))
    ret[:, 4 + footId * 3, 4 + footId * 3 +2] = identity(2)
    return ret

def initGlobals(nEffectors, comWeightsPerEffector=None):
    global DEFAULT_NUM_VARS
    global N_EFFECTORS
    global DEFAULT_NUM_EQUALITY_CONSTRAINTS
    global DEFAULT_NUM_INEQUALITY_CONSTRAINTS
    global W_START
    global COM_WEIGHTS_PER_EFFECTOR
    
    N_EFFECTORS = nEffectors
    if comWeightsPerEffector is None:
        COM_WEIGHTS_PER_EFFECTOR = [1./(N_EFFECTORS -1) for _ in range(N_EFFECTORS) ]
    else:
        COM_WEIGHTS_PER_EFFECTOR = comWeightsPerEffector
    # for one phase, for all i:
    # [ c_x't, c_y't, c_z_1't, c_z_2't, p_i't, b_i't, g_i't, w_i_t ,  a_i_l't]   
    # [ 1    , 1    ,   1    ,    1   , 3*n_e, 2*n_e, 3*n_e,    n_e, [0, n_s]]   
    DEFAULT_NUM_VARS = N_EFFECTORS * 9 + 4     
    BETA_START = 4 + 3* N_EFFECTORS
    GAMMA_START = BETA_START +  2 * N_EFFECTORS
    W_START = GAMMA_START +  3 * N_EFFECTORS
    ALPHA_START = W_START +  N_EFFECTORS
    assert ALPHA_START == DEFAULT_NUM_VARS
    
    # c_x't, c_y't are given by a convex combination of fixed weight of the non moving effectors. There are N_EFFECTORS possible combinations
    # by default each position is equal to the previous one, so that 3 * N_EFFECTORS equalities 
    DEFAULT_NUM_EQUALITY_CONSTRAINTS = (2 + 3) * N_EFFECTORS
    # wi <= 1 ;  - M_wi <= b_ix + y <= M w_i; - M_wi <= g_ix + y + z <= M w_i      wi> 0 implicit
    DEFAULT_NUM_INEQUALITY_CONSTRAINTS = (1 + 1 + 1) * N_EFFECTORS
    
    
    global COM_XY_ExpressionMatrix 
    global COM_Z1_ExpressionMatrix 
    global COM_Z2_ExpressionMatrix 
    global COM1_ExpressionMatrix    
    global COM2_ExpressionMatrix    
    global COM_XY_ExpressionMatrix 
    global foot_ExpressionMatrix   
    global foot_XY_ExpressionMatrix
    
    COM_XY_ExpressionMatrix  = _COM_XY_ExpressionMatrix()
    COM_Z1_ExpressionMatrix  = _COM_Z1_ExpressionMatrix()
    COM_Z2_ExpressionMatrix  = _COM_Z2_ExpressionMatrix()
    COM1_ExpressionMatrix     = _COME1xpressionMatrix()
    COM2_ExpressionMatrix     = _COME2xpressionMatrix()
    foot_ExpressionMatrix    = [_footExpressionMatrix    (footId) for footId in range(N_EFFECTORS)]
    foot_XY_ExpressionMatrix = [_foot_XY_ExpressionMatrix(footId) for footId in range(N_EFFECTORS)]
    

### helper functions to count number of variables, equalities and inequalities constraints in the problem ###

def numVariablesForPhase(phase):
    ret = DEFAULT_NUM_VARS
    numSurfaces =len(phase["S"])
    if numSurfaces > 1:
        ret += numSurfaces + numSurfaces * N_EFFECTORS
    return ret
    
# ~ def numIneqForPhase(phase, phase = -1 ):
def numIneqForPhase(phase, phaseId):
    #surface -one because last is equality
    #COM Kinematic constraints: summation over all effectors, times 2 because there are 2 height possible for the transition
    #except for first time
    ret = sum([k.shape[0]  for (_, k) in phase["K"][0] ])
    if phaseId != 0:
        ret *= 2
    # relative kinematic constraints between each effectors
    for footIdFrame, constraintsInFootIdFrame in enumerate(phase["allRelativeK"][0]):
        for (footId, Kk ) in  constraintsInFootIdFrame:
            ret += Kk[0].shape[0]    
    # all inequalities relative to each contact surface 
    # S + 1 because equality becomes 2 inequalities
    # the inequalities must always be written for each effector
    ret += sum([(S[1].shape[0]) * N_EFFECTORS for S in  phase["S"]])
    numSurfaces =len(phase["S"])
    if numSurfaces >1:
        # alpha > 0 
        ret += numSurfaces * NUM_INEQ_SLACK_PER_SURFACE 
    ret += DEFAULT_NUM_INEQUALITY_CONSTRAINTS
    return ret
    
def getTotalNumVariablesAndIneqConstraints(pb):
    nphases = pb["nphases"]
    cols = sum([numVariablesForPhase(phase) for phase in  pb["phaseData"]])
    rows = sum([numIneqForPhase(phase, i) for  i, phase in enumerate(pb["phaseData"]) ])
    return rows, cols
    
def numEqForPhase(phase, phaseId):
    #no equality constraints at first phase
    if phaseId == 0:
        return 0.
    return DEFAULT_NUM_EQUALITY_CONSTRAINTS
    
def getTotalNumEqualityConstraints(pb):
    raise ValueError # initphase ?
    return sum( [numEqForPhase(phase) for  phase in pb["phaseData"] ])
    

### Constraint functions ###
        
# for all effectors i , Ki (c2 - pi) <= ki
def FootCOM2KinConstraint(pb, phaseDataT, A, b, startCol, endCol, startrow):
    idRow = startRow
    for footId, (K, k) in enumerate(phaseDataT["K"][0]):
        consLen = K.shape[0]
        A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS] =  K.dot(COM_Z2_ExpressionMatrix - foot_ExpressionMatrix[footId])
        b[idRow:idRow+consLen] = k 
        idRow += consLen    
    return idRow
    
# for all effectors i , Ki (c1't - pi'(t-1)) <= ki
# this should not be called for first phase
def FootCOM1KinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startrow):
    idRow = startRow
    for footId, (K, k) in enumerate(phaseDataT["K"][0]):
        consLen = K.shape[0]
        A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS]                  =  K.dot(COM_Z1_ExpressionMatrix)
        A[idRow:idRow+consLen, previousStartCol:previousStartCol+DEFAULT_NUM_VARS]  = -K.dot(foot_ExpressionMatrix[footId])
        b[idRow:idRow+consLen] = k 
        idRow +=   consLen    
    return idRow
    
def FootCOMKinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    idRow = startRow
    if startRow != 0:
        idRow = FootCOM1KinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startrow)
    return FootCOM2KinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startrow)
        # ~ if pb["p0"] is None:
            # ~ fixed = phaseDataT["fixed"] #0 constraints
            # ~ K, k = phaseDataT["K"][0][fixed]
            # ~ return K.shape[0]
        # ~ else:
            # ~ return FixedFootCOMKinConstraintInitPhase(pb, phaseDataT, A, b, startCol, endCol)
        
# for all effectors i , for all j !=i Ki (pj - pi) <= ki   
# TODO REMOVE FOR FIRST PHASE
def RelativeDistanceConstrain(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):    
    idRow = startRow
    for footIdFrame, constraintsInFootIdFrame in enumerate(phaseDataT["allRelativeK"][0]):
        for (footId, Kk ) in  constraintsInFootIdFrame:
            K = Kk[0]; k = Kk[1]
            consLen = K.shape[0]
            A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS] = K.dot(foot_ExpressionMatrix[footId] - foot_ExpressionMatrix[footIdFrame])
            b[idRow:idRow+consLen] = k 
            idRow += consLen    
    return idRow       
    
def SurfaceConstraint(phaseDataT, A, b, startCol, endCol, startRow):    
    sRow = startRow
    nSurfaces = len(phaseDataT["S"])
    idS = ALPHA_START
    for (S,s) in phaseDataT["S"]:   
        for footId in range(N_EFFECTORS):
            idRow = sRow + S.shape[0]
            # Sl pi - M alphal <= sl + M2 (1 - w_t)
            # Sl pi - M alphal + M2 w_t <= sl + M2
            onesM2 = ones(idRow-sRow) * M2
            A[sRow:idRow, startCol:startCol+DEFAULT_NUM_VARS] = S.dot(foot_ExpressionMatrix[footId])
            A[sRow:idRow, startCol+W_START + footId] = onesM2
            b[sRow:idRow                 ] = s + onesM2
            if nSurfaces >1:
                A[sRow:idRow, startCol+idS] = -ones(idRow-sRow) * M
            sRow = idRow
        idS += 1
    return idRow
    
def SlackPositivityConstraint(phaseDataT, A, b, startCol, endCol, startRow):     
    idRow = startRow
    for footId in range(N_EFFECTORS):
        # -Mwi <= b_i x + b_i y <= M wi  
        # -Mwi  - b_i x - b_i y} <= 0 ; b_ix + b_iy- M wi <= 0
        A[idRow, startCol + W_START + footId                                               ] = [-M]; 
        A[idRow, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = [-1, -1];
        idRow += 1 
        A[idRow, startCol + W_START + footId                                               ] = [-M]; 
        A[idRow, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = [1, 1]; 
        idRow += 1 
        # -Mwi <= g_i x + g_i y + g_i_z <= M wi  
        # -Mwi  - g_i x - g_i y - g_i_z  <= 0 ; g_ix + g_iy + g_iz - M wi <= 0
        A[idRow, startCol + W_START + footId                                               ] = [-M]; 
        A[idRow, (startCol + GAMMA_START + footId*3):(startCol + GAMMA_START + footId*3) + 3 ] = [-1, -1, -1];         
        idRow += 1 
        A[idRow, startCol + W_START + footId                                               ] = [-M]; 
        A[idRow, (startCol + GAMMA_START + footId*3):(startCol + GAMMA_START + footId*3) + 3 ] = [1, 1, 1];        
        idRow += 1 
        # wi <= 1
        A[idRow, startCol + W_START + footId ] = 1; 
        b[idRow] = 1.; 
    # -al < 0
    nSurfaces = len(phaseDataT["S"])
    if nSurfaces > 1:
        for i in range(nSurfaces):
             A[idRow+i, startCol + ALPHA_START + i ] = -1;
         idRow += nSurfaces   
    return idRow    
    
def CoMWeightedEqualityConstraint(phaseDataT, E, e, startCol, endCol, startRow):
    idRow = startRow
    for flyingFootId in range(N_EFFECTORS):
        # 0 =  sum(j != i) o_j'i p_j't + [b_ix't, b_iy't]^T - c_x_y
        EqMat = -COM_XY_ExpressionMatrix[:]
        for otherFootId in range(N_EFFECTORS):
            if flyingFootId != otherFootId:
                EqMat += COM_WEIGHTS_PER_EFFECTOR[otherFootId] * footExpressionMatrixXY[otherFootId]
        EqMat[idRow:idRow+2, startCol + BETA_START + flyingFootId*2:(startCol + BETA_START + flyingFootId*2)+2] = identity(2)
        E[idRow:idRow+2, startCol:startCol+DEFAULT_NUM_VARS] = EqMat; #e = 0 
        idRow+=2
    return idRow    
    
#only applies after first step
def FootContinuityEqualityConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    idRow = startRow
    for footId in range(N_EFFECTORS):
        # 0 =  p_(i-1)'t - p_(i-1)'t + [g_ix't, b_iy't, g_iz't]^T
        E[idRow:idRow+3, startCol:startCol + DEFAULT_NUM_VARS ]         =  footExpressionMatrix[footId]
        E[idRow:idRow+3, previousStartCol:previousStartCol + DEFAULT_NUM_VARS] = -footExpressionMatrix[footId]
        E[idRow:idRow+3, startCol + GAMMA_START + footId*3:(startCol + BETA_START + flyingFootId*3)+3] = identity(3); #e = 0 
        idRow+=3
    return idRow    
        
    
def convertProblemToLp(pb, convertSurfaces = True):        
    assert DEFAULT_NUM_VARS is not None, "call initGlobals first"
    if convertSurfaces:
        replace_surfaces_with_ineq_in_problem(pb, eqAsIneq = True)
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
    # nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    indices = []
    cIdx = 0
    wrongsurfaces = []
    wrongsurfaces_indices = []
    for i, phase in enumerate(pb["phaseData"]):  
        numSurfaces = len(phase["S"])
        phaseVars = numVariablesForPhase(phase)
        if numSurfaces > 1:
            startIdx = cIdx + DEFAULT_NUM_VARS
            betas = [res[startIdx+j] for j in range(0,numSurfaces*2,2) ]
            if array(betas).min() > 0.01: ####
                # print "wrong ", i, array(betas).min()
                indices += [i]
                sorted_surfaces = np.argsort(betas)
                wrongsurfaces_indices += [sorted_surfaces]
                #~ print "sorted_surfaces ",sorted_surfaces
                wrongsurfaces += [[[phase["S"][idx]] for idx in sorted_surfaces]  ]
                #~ print "lens ", len([[phase["S"][idx]] for idx in sorted_surfaces]  )
        cIdx += phaseVars
    return indices, wrongsurfaces, wrongsurfaces_indices
    
def isSparsityFixed(pb, res):
    indices, wrongsurfaces, wrongsurfaces_indices = num_non_zeros(pb, res)
    return len(indices) == 0
   
def genOneComb(pb,indices, surfaces, res):
    pb1 = pb.copy()
    for i, idx in enumerate(indices):
        pb1["phaseData"][idx]["S"] = surfaces[i][0]
    res += [pb1]
   
import itertools
import copy
   
def genCombinatorialRec(pb, indices, wrongsurfaces, wrongsurfaces_indices, res):
    lenss  = [len(surfs) for surfs in wrongsurfaces]
    all_indexes = [[el for el in range(lens)]  for lens in [len(surfs) for surfs in wrongsurfaces]]
    wrong_combs = [el for el in itertools.product(*wrongsurfaces_indices)]
    combs = [el for el in itertools.product(*all_indexes)]
    for j, comb in enumerate(combs):
        pb1 = copy.deepcopy(pb)
        for i, idx in enumerate(indices):
            pb1["phaseData"][idx]["S"] = wrongsurfaces[i][comb[i]]
        res += [[pb1, wrong_combs[j], indices]]
        
    
    
def generateAllFixedScenariosWithFixedSparsity(pb, res):
    indices, wrongsurfaces, wrongsurfaces_indices = num_non_zeros(pb, res)
    all_len = [len(surfs) for surfs in wrongsurfaces]
    comb = 1
    for el in all_len:
        comb *= el  
    res = []
    if comb >1000:
        print ("problem probably too big ", comb)
        return 1
    else:
        genCombinatorialRec(pb, indices, wrongsurfaces, wrongsurfaces_indices, res)
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
            assert min(betas) >= -0.00000001
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
   
   
def plotConstraints(ax, pb, allfeetpos, coms):
    from tools.plot_plytopes import plot_polytope_H_rep
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
                print("qhullfailed")
    
        
def plotQPRes(pb, res, linewidth=2, ax = None, plot_constraints = False, show = True):
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)
        
    ax.set_autoscale_on(False)
    ax.view_init(elev=8.776933438381377, azim=-99.32358055821186)
    
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
    
    print("time to set up problem" , timMs(t1,t2))
    print("time to solve problem"  , timMs(t2,t3))
    print("total time"             , timMs(t1,t3))
    
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
    
    print("time to set up problem" , timMs(t1,t2))
    print("time to solve problem"  , timMs(t2,t3))
    print("total time"             , timMs(t1,t3))
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False)
        
    
    
    
    
        
