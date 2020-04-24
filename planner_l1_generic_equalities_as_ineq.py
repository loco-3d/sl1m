

############### Problem definition #############

from sl1m.problem_definition import *
           
        
#### problem related global variables #### 

NUM_SLACK_PER_SURFACE = 1
NUM_INEQ_SLACK_PER_SURFACE = 1 # 
M = 10.
M2= M*10
# ~ M2= M
        
        
#### global variables, depending on number of effectors ####  
# You need to call initGlobals with the number of effectors before using sl1m generic
N_EFFECTORS                             = None
DEFAULT_NUM_VARS                        = None
COM_WEIGHTS_PER_EFFECTOR                = None
DEFAULT_NUM_EQUALITY_CONSTRAINTS        = None
DEFAULT_NUM_EQUALITY_CONSTRAINTS_START  = None
DEFAULT_NUM_INEQUALITY_CONSTRAINTS      = None

#### Selection matrices ####  
COM_XY_ExpressionMatrix  = None
COM_Z1_ExpressionMatrix  = None
COM_Z2_ExpressionMatrix  = None
COM1_ExpressionMatrix    = None
COM2_ExpressionMatrix    = None
footExpressionMatrix    = None
footExpressionMatrixXY  = None
wExpressionMatrix       = None
gExpressionMatrix       = None
bExpressionMatrix       = None

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
          
def _COM1ExpressionMatrix():
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:,:3] = identity(3)
    return ret
    
def _COM2ExpressionMatrix():
    ret = zeros((3, DEFAULT_NUM_VARS))
    ret[:2,:2] = identity(2)
    ret[2,  3] = 1
    return ret
        
def _footExpressionMatrix(footId, dim):
    ret = zeros((dim, DEFAULT_NUM_VARS))
    ret[:, 4 + footId * 3:4 + footId * 3 +dim] = identity(dim)
    return ret
        
def _contactDecisionExpressionMatrix():
    ret = zeros(DEFAULT_NUM_VARS)
    ret[W_START:] = ones(N_EFFECTORS)
    return ret
    
def _gExpressionMatrix():
    ret = zeros(DEFAULT_NUM_VARS)
    ret[GAMMA_START:W_START] = ones(3 * N_EFFECTORS)
    return ret
    
def _bExpressionMatrix():
    ret = zeros(DEFAULT_NUM_VARS)
    ret[BETA_START:GAMMA_START] = ones(2 * N_EFFECTORS)
    return ret
    
def initGlobals(nEffectors, comWeightsPerEffector=None):
    global DEFAULT_NUM_VARS
    global N_EFFECTORS
    global DEFAULT_NUM_EQUALITY_CONSTRAINTS
    global DEFAULT_NUM_EQUALITY_CONSTRAINTS_START
    global DEFAULT_NUM_INEQUALITY_CONSTRAINTS
    global COM_WEIGHTS_PER_EFFECTOR    
    global BETA_START 
    global GAMMA_START
    global W_START    
    global ALPHA_START
    
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
    DEFAULT_NUM_EQUALITY_CONSTRAINTS_START = (2) * N_EFFECTORS
    # wi <= 1 ;  - M_wi <= b_ix _ y <= M w_i; - M_wi <= g_ix _ y _ z <= M w_i      wi> 0 implicit
    DEFAULT_NUM_INEQUALITY_CONSTRAINTS = (1 + 2 * 2 + 2 * 3) * N_EFFECTORS
    
    
    global COM_XY_ExpressionMatrix 
    global COM_Z1_ExpressionMatrix 
    global COM_Z2_ExpressionMatrix  
    global COM_XY_ExpressionMatrix 
    global footExpressionMatrix   
    global footExpressionMatrixXY   
    global COM1_ExpressionMatrix  
    global COM2_ExpressionMatrix  
    global wExpressionMatrix   
    global gExpressionMatrix   
    global bExpressionMatrix   
    
    COM_XY_ExpressionMatrix  = _COM_XY_ExpressionMatrix()
    COM_Z1_ExpressionMatrix  = _COM_Z1_ExpressionMatrix()
    COM_Z2_ExpressionMatrix  = _COM_Z2_ExpressionMatrix()
    COM1_ExpressionMatrix    = _COM1ExpressionMatrix()
    COM2_ExpressionMatrix    = _COM2ExpressionMatrix()
    wExpressionMatrix        = _contactDecisionExpressionMatrix()
    gExpressionMatrix        = _gExpressionMatrix()
    bExpressionMatrix        = _bExpressionMatrix()
    footExpressionMatrix     = [_footExpressionMatrix(footId, 3) for footId in range(N_EFFECTORS)]
    footExpressionMatrixXY   = [_footExpressionMatrix(footId, 2) for footId in range(N_EFFECTORS)]
    

### helper functions to count number of variables, equalities and inequalities constraints in the problem ###

def numVariablesForPhase(phase):
    ret = DEFAULT_NUM_VARS
    numSurfaces =len(phase["S"])
    if numSurfaces > 1:
        ret += numSurfaces
    return ret
    
# ~ def numIneqForPhase(phase, phase = -1 ):
def numIneqForPhase(phase, phaseId):
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
    # the inequalities must always be written for each effector
    ret += sum([(S[1].shape[0]) * N_EFFECTORS for S in  phase["S"]])
    numSurfaces =len(phase["S"])
    if numSurfaces >1:
        # alpha > 0 
        ret += numSurfaces * NUM_INEQ_SLACK_PER_SURFACE 
        # ~ ret += numSurfaces * 2 * N_EFFECTORS 
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
        return DEFAULT_NUM_EQUALITY_CONSTRAINTS_START
    return DEFAULT_NUM_EQUALITY_CONSTRAINTS
    
def getTotalNumEqualityConstraints(pb):
    return int(sum( [numEqForPhase(phase, phaseId) for  phaseId, phase in enumerate(pb["phaseData"]) ]))
    

### Constraint functions ###
        
# for all effectors i , Ki (c2 - pi) <= ki
def FootCOM2KinConstraint(pb, phaseDataT, A, b, startCol, endCol, startRow):
    idRow = startRow
    for footId, (K, k) in enumerate(phaseDataT["K"][0]):
        consLen = K.shape[0]
        A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS] =  K.dot(COM2_ExpressionMatrix - footExpressionMatrix[footId])
        b[idRow:idRow+consLen] = k 
        idRow += consLen    
    return idRow
    
# for all effectors i , Ki (c1't - pi'(t-1)) <= ki
# this should not be called for first phase
def FootCOM1KinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow):
    idRow = startRow
    for footId, (K, k) in enumerate(phaseDataT["K"][0]):
        consLen = K.shape[0]
        A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS]                  =  K.dot(COM1_ExpressionMatrix)
        A[idRow:idRow+consLen, previousStartCol:previousStartCol+DEFAULT_NUM_VARS]  = -K.dot(footExpressionMatrix[footId])
        b[idRow:idRow+consLen] = k 
        idRow +=   consLen    
    return idRow
    
def FootCOMKinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow, phaseId):
    idRow = startRow
    if phaseId != 0:
        idRow = FootCOM1KinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow)
    return FootCOM2KinConstraint(pb, phaseDataT, A, b, startCol, endCol, idRow)
        
# for all effectors i , for all j !=i Ki (pj - pi) <= ki   
# TODO REMOVE FOR FIRST PHASE
def RelativeDistanceConstraint(pb, phaseDataT, A, b, startCol, endCol, startRow, phaseId):    
    idRow = startRow
    # ~ from tools.plot_plytopes import plot_polytope_H_rep
    limbNames = ['LFleg', 'RFleg', 'LHleg', 'RHleg']
    for footIdFrame, constraintsInFootIdFrame in enumerate(phaseDataT["allRelativeK"][0]):
        # ~ fig = plt.figure()
        # ~ fig.suptitle(limbNames[footIdFrame], fontsize=16)
        # ~ ax = fig.add_subplot(111, projection="3d")
        for (footId, Kk ) in  constraintsInFootIdFrame:
            K = Kk[0]; k = Kk[1]
            consLen = K.shape[0]
            # ~ print ("KK", K.shape)
            # ~ succPlot = plot_polytope_H_rep(K,k.reshape((-1,1)), ax = ax)
            # ~ if not (succPlot):
                # ~ print ("indices ", footIdFrame, footId)
            A[idRow:idRow+consLen, startCol:startCol+DEFAULT_NUM_VARS] = K.dot(footExpressionMatrix[footId] - footExpressionMatrix[footIdFrame])
            b[idRow:idRow+consLen] = k[:]
            idRow += consLen   
    return idRow       
    
def SurfaceConstraint(phaseDataT, A, b, startCol, endCol, startRow):    
    idRow = startRow
    sRow = startRow
    nSurfaces = len(phaseDataT["S"])
    idS = ALPHA_START
    for (S,s) in phaseDataT["S"]:   
        for footId in range(N_EFFECTORS):
            idRow = sRow + S.shape[0]
            # Sl pi - M alphal <= sl + M2 (1 - w_t)
            # Sl pi - M alphal + M2 w_t <= sl + M2
            onesM2 = ones(idRow-sRow) * M2
            A[sRow:idRow, startCol:startCol+DEFAULT_NUM_VARS] = S.dot(footExpressionMatrix[footId])
            A[sRow:idRow, startCol+W_START + footId] = onesM2
            b[sRow:idRow                 ] = s + onesM2
            if nSurfaces >1:
                A[sRow:idRow, startCol+idS] = -ones(idRow-sRow) * M
            sRow = idRow
        idS += 1
    return idRow
    
def SlackPositivityConstraint(phaseDataT, A, b, startCol, endCol, startRow):     
    idRow = startRow
    nSurfaces = len(phaseDataT["S"])
    for footId in range(N_EFFECTORS):
        # -Mwi <= b_i x y <= M wi  
        # -Mwi  - b_i xy} <= 0 ; b_ixy- M wi <= 0
        # ~ A[idRow:idRow+2, startCol + W_START + footId                                               ] = ones(2)*-M; 
        # ~ A[idRow:idRow+2, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = -identity(2);
        # ~ idRow += 2 
        # ~ A[idRow:idRow+2, startCol + W_START + footId                                               ] = ones(2)*-M; 
        # ~ A[idRow:idRow+2, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = identity(2);
        # ~ idRow += 2
        
        # -M (1 - wi) <= b_i x y <= M (1 - wi)
        # -M  + Mwi <= b_i x y <= M - Mwi
        #  Mwi - b <= M;  <= Mwi + b <= M
        A[idRow:idRow+2, startCol + W_START + footId                                               ] = ones(2)*M; 
        A[idRow:idRow+2, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = -identity(2);
        b[idRow:idRow+2] = ones(2)*M; 
        idRow += 2 
        A[idRow:idRow+2, startCol + W_START + footId                                               ] = ones(2)*M; 
        A[idRow:idRow+2, (startCol + BETA_START + footId*2):(startCol + BETA_START + footId*2) + 2 ] = identity(2);
        b[idRow:idRow+2] = ones(2)*M; 
        idRow += 2
        # -Mwi <= g_i x + g_i y + g_i_z <= M wi  
        # -Mwi  - g_i x - g_i y - g_i_z  <= 0 ; g_ix + g_iy + g_iz - M wi <= 0
        A[idRow:idRow+3, startCol + W_START + footId                                                 ] = ones(3)*-M; 
        A[idRow:idRow+3, (startCol + GAMMA_START + footId*3):(startCol + GAMMA_START + footId*3) + 3 ] = -identity(3);       
        idRow += 3 
        A[idRow:idRow+3, startCol + W_START + footId                                                 ] = ones(3)*-M; 
        A[idRow:idRow+3, (startCol + GAMMA_START + footId*3):(startCol + GAMMA_START + footId*3) + 3 ] = identity(3);        
        idRow += 3 
        # wi <= 1
        A[idRow, startCol + W_START + footId ] = 1; 
        b[idRow] = 1.00;     
        idRow += 1 
        # ~ if nSurfaces > 1:            
            # -Mwi <= a_l <= M wi  
            # -Mwi  - a_l  <= 0 ; a_l - M wi <= 0
            # ~ for i in range(nSurfaces):
                 # ~ A[idRow  , startCol + ALPHA_START + i ]  = -1;
                 # ~ A[idRow  , startCol + W_START + footId ] = -M;
                 # ~ A[idRow+1, startCol + ALPHA_START + i ]  =  1;
                 # ~ A[idRow+1, startCol + W_START + footId ] = -M;
                 # ~ idRow += 2
    # -al < 0
    # -al < 0
    # ~ nSurfaces = len(phaseDataT["S"])
    if nSurfaces > 1:
        for i in range(nSurfaces):
             A[idRow+i, startCol + ALPHA_START + i ] = -1;
        idRow += nSurfaces   
    return idRow    
    
def CoMWeightedEqualityConstraint(phaseDataT, E, e, startCol, endCol, startRow):
    idRow = startRow
    for flyingFootId in range(N_EFFECTORS):
        # 0 =  sum(j != i) o_j'i p_j't + [b_ix't, b_iy't]^T - c_x_y
        # one beta per equality dum dum TODO
        EqMat = -COM_XY_ExpressionMatrix[:]
        for otherFootId in range(N_EFFECTORS):
            if flyingFootId != otherFootId:
                EqMat += COM_WEIGHTS_PER_EFFECTOR[otherFootId] * footExpressionMatrixXY[otherFootId]
        EqMat[:2, BETA_START + flyingFootId*2:(BETA_START + flyingFootId*2)+2] = identity(2)
        E[idRow:idRow+2, startCol:startCol+DEFAULT_NUM_VARS] = EqMat; #e = 0 
        idRow+=2
    return idRow    
    
#only applies after first step
def FootContinuityEqualityConstraint(pb, phaseDataT, E, e, previousStartCol, startCol, endCol, startRow, phaseId):
    idRow = startRow
    if phaseId !=0:
        for footId in range(N_EFFECTORS):
            # 0 =  p_(i-1)'t - p_(i-1)'t + [g_ix't, b_iy't, g_iz't]^T
            E[idRow:idRow+3, startCol:startCol + DEFAULT_NUM_VARS ]         =  footExpressionMatrix[footId]
            E[idRow:idRow+3, previousStartCol:previousStartCol + DEFAULT_NUM_VARS] = -footExpressionMatrix[footId]
            E[idRow:idRow+3, startCol + GAMMA_START + footId*3:(startCol + GAMMA_START + footId*3)+3] = identity(3); #e = 0 
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
    
    for phaseId, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        startRow = FootCOMKinConstraint(pb, phaseDataT, A, b, previousStartCol, startCol, endCol, startRow, phaseId)
        startRow = RelativeDistanceConstraint(pb, phaseDataT, A, b, startCol, endCol, startRow, phaseId)
        startRow = SurfaceConstraint(phaseDataT, A, b, startCol, endCol, startRow)
        startRow = SlackPositivityConstraint(phaseDataT, A, b, startCol, endCol, startRow)
        
        #equalities     
        #no weighted com on first phase   
        # ~ if phaseId != 0:
            # ~ startRowEq = CoMWeightedEqualityConstraint(phaseDataT, E, e, startCol, endCol, startRowEq)
        startRowEq = FootContinuityEqualityConstraint(pb, phaseDataT, E, e, previousStartCol, startCol, endCol, startRowEq, phaseId)
        previousStartCol = startCol
        startCol   = endCol 
        
    print ('startRow ', startRow)
    print ('A.shape[0] ', A.shape[0])
        
    # ~ assert endCol   == A.shape[1]
    # ~ assert startRowEq == E.shape[0]
    # ~ assert startRow == A.shape[0]
    
    A,b = normalize([A,b])
    E,e = normalize([E,e])
    return (A,b,E,e)
        
def pbSelectionMatrix(pb, selectionMatrix = wExpressionMatrix):
    nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    c = zeros(nvars)
    cIdx = 0
    for i, phase in enumerate(pb["phaseData"]):
        c[cIdx:cIdx+DEFAULT_NUM_VARS] = selectionMatrix[:]
        cIdx += numVariablesForPhase(phase)
    assert cIdx == nvars
    return c

#contact activations
def wSelectionMatrix(pb):
    return pbSelectionMatrix(pb, selectionMatrix = wExpressionMatrix)
    
#position continuity violation
def gSelectionMatrix(pb):
    return pbSelectionMatrix(pb, selectionMatrix = gExpressionMatrix)
    
#com centering
def bSelectionMatrix(pb):
    return pbSelectionMatrix(pb, selectionMatrix = bExpressionMatrix)

def alphaSelectionMatrix(pb):
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
    
def gammaSelectionMatrix(pb):
    nvars = getTotalNumVariablesAndIneqConstraints(pb)[1]
    c = zeros(nvars)
    cIdx = 0
    for i, phase in enumerate(pb["phaseData"]):
        c[cIdx:cIdx+DEFAULT_NUM_VARS] = g[:]
        print ("wExpressionMatrix", wExpressionMatrix)
        cIdx += numVariablesForPhase(phase)
    assert cIdx == nvars
    return c
    return c

###########################" PLOTTING ################"
    
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def retrieve_points_from_res(pb, res):
    coms = []
    # ~ if pb["c0"] is not None:
        # ~ coms = [pb["c0"]]
        
    # ~ if pb["p0"] is not None:
        # ~ footPos = [[pb["p0"][LF]],[pb["p0"][RF]]]
        # ~ allFeetPos = [footPos[0][0], footPos[1][0]]
    # ~ else:
    footPos = [[] for _ in range(N_EFFECTORS)]
    allFeetPos = []
        
    col = 0
    for i, phaseDataT in enumerate(pb["phaseData"]):          
        coms += [COM1_ExpressionMatrix.dot(res[col:col + COM1_ExpressionMatrix.shape[1]]), COM2_ExpressionMatrix.dot(res[col:col + COM2_ExpressionMatrix.shape[1]])]
        allFeetPhase = []
        
        for footId in range(N_EFFECTORS):
            posE = array(footExpressionMatrix[footId].dot(res[col:col + footExpressionMatrix[footId].shape[1]]) )
            pos = posE * res[col + W_START + footId]
            # ~ pos = pos +  array([11,11,0.])
            # ~ allFeetPos += [allFeetPhase]
            # ~ if norm(pos) > 0.01:
            # ~ allFeetPhase += [pos]
            # ~ footPos[footId] += [pos]
            footPos[footId] += [posE]
                # ~ allFeetPos += [pos]
            
            print("posE", posE)
            allFeetPhase += [posE]
        print ("allFeetPhase" , allFeetPhase)
        allFeetPos += [allFeetPhase]
        print ("allFeetPos" , allFeetPos)
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
    
from scipy.spatial import ConvexHull   
from sl1m.tools.plot_plytopes import plot_hull   
def plotQPRes(pb, res, linewidth=2, ax = None, plot_constraints = False, show = True, plotSupport = False):
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)    
    ax.set_xlim3d(0, 2)
    ax.set_ylim3d(0,2)
    ax.set_zlim3d(0,2)
        
    ax.set_autoscale_on(False)
    ax.view_init(elev=8.776933438381377, azim=-99.32358055821186)
    
    from sl1m.tools.plot_plytopes import plot_polytope_H_rep
    from collections import namedtuple
    Ineq = namedtuple("Inequality", "A b")
    
    colors = ["r","b","g","y"]
    Af = None
    bf = None
    
    weighted = []
    for footId, (K, k) in enumerate(pb["phaseData"][0]["K"][0]):
        Afoot = K
        bfoot = (k + K.dot(footpos[footId][0])).reshape([-1,1])
        if Af is None:
            Af = K; bf = bfoot
        else:
            Af = np.vstack([Af, Afoot])
            bf = np.vstack([bf, bfoot])
        # ~ print ("footposid", footpos[footId][0])
        # ~ bfoot = k.reshape([-1,1])
        # ~ print ("range", (K.dot([-1,0.,0.]) - k).max())
        pond = zeros(3); 
        for otherFootId, positions in enumerate(footpos):
            if otherFootId != footId:
                pond += 0.33 * positions[0]
        pond[-1] = 0.4
        weighted += [pond]
    plot_polytope_H_rep(Af,bf, color = "r", just_pts = False, ax = ax)
        # ~ plot(Ineq(Afoot,k), ax = ax, show = False, color = colors[footId])
    # ~ [print("wtf", (Af[:,:2].dot(el[:2])-bf).max() ) for el in weighted]
    
    
    Afreez = Af[:,2:3]
    for el in weighted:
        bFreez = bf.reshape((-1,)) - Af[:,:2].dot(el[:2])
        # ~ print ('AAAAAAAAAAAAAAAAAAA',bf.shape)
        from sl1m.qp import solve_lp
        z = solve_lp(zeros(1),Afreez, bFreez)
        el[-1] = z.x
        plotPoints(ax, [el], color = "b")
        
    
    # ~ plotPoints(ax, coms[1:], color = "b")
    
    footColors = ['r', 'g', 'y', 'b']
    
    # ~ plotPoints(ax, coms, )
    [plotPoints(ax, footpos[idx], color = footColors[idx]) for idx in range(N_EFFECTORS)]
    
    cx = [c[0] for c in coms[1:]]
    cy = [c[1] for c in coms[1:]]
    cz = [c[2] for c in coms[1:]]
    ax.plot(cx, cy, cz)
    idC = 0;
    if plotSupport:
        for el in allfeetpos[:]:
            # ~ print("points ", [el[:2] for el in allfeetpos])
            pts = [e[:2] for e in el]
            apts = array(pts)
            chul = ConvexHull(array([e[:2] for e in el]))
            plot_hull(chul, pts = pts, apts = apts, ax = ax, color = footColors[idC])
            idC = (idC + 1) % 4
        # ~ px = [c[0] for c in el]
        # ~ py = [c[1] for c in el]
        # ~ pz = [c[2] for c in el]
        # ~ ax.plot(px, py, pz)
    # ~ px = [c[0] for c in allfeetpos]
    # ~ py = [c[1] for c in allfeetpos]
    # ~ pz = [c[2] for c in allfeetpos]
    # ~ ax.plot(px, py, pz)
        
    if show:
        plt.ion()
        plt.show()
       
    
    
# ~ def addInitEndConstraint(pb, E, e, posInit= array([0.16, 0.585, 0.0]), posEnd = array([1.2, 0.37, 0.40])):
def addInitEndConstraint(pb, E, e, posInit, posEnd, initCom, endCom):
    sizeAdd = 0
    if posInit is not None:
        sizeAdd += len(posInit) * 3
    if posEnd is not None:
        sizeAdd += len(posEnd) * 3
    if initCom is not None:
        sizeAdd +=  3
    if endCom is not None:
        sizeAdd +=  3
        
    nE = zeros((E.shape[0] +sizeAdd, E.shape[1] ))
    ne = zeros(E.shape[0] +sizeAdd)
    idRow = E.shape[0]
    nE[:idRow,:E.shape[1]] = E  
    ne[:idRow] = e[:]    
    
    print ("heo", E.shape)
    
    idRow = E.shape[0]
    nVarEnd = numVariablesForPhase(pb["phaseData"][-1])
    
    if initCom is not None:
        nE[idRow:idRow+3,:DEFAULT_NUM_VARS] = COM2_ExpressionMatrix[:] 
        ne[idRow:idRow+3] = initCom
        idRow+=3
    if endCom is not None:
        print("nVarEnd", nVarEnd, E.shape[1]-nVarEnd, E.shape[1]-nVarEnd+DEFAULT_NUM_VARS)
        nE[idRow:idRow+3:,E.shape[1]-nVarEnd:E.shape[1]-nVarEnd+DEFAULT_NUM_VARS] = COM1_ExpressionMatrix[:] 
        ne[idRow:idRow+3] = endCom
        idRow+=3
    if posInit is not None:
        print("adding ", idRow)
        for idFoot, pos in enumerate(posInit):
            print("adding ", idRow )
            nE[idRow:idRow+3,:DEFAULT_NUM_VARS] = footExpressionMatrix[idFoot][:]  
            ne[idRow:idRow+3] = pos
            idRow+=3
    if posEnd is not None:        
        for idFoot, pos in enumerate(posEnd):
            print("adding end", idRow )
            nE[idRow:idRow+3:,E.shape[1]-nVarEnd:E.shape[1]-nVarEnd+DEFAULT_NUM_VARS] = footExpressionMatrix[idFoot][:]
            ne[idRow:idRow+3] = pos
            idRow+=3
    return nE, ne
    
# ~ def addInitEndConstraint(pb, E, e, posInit= array([0.0, 0.0, 0.4]), posEnd = array([-0.0, 0.0, 0.4])):
    # ~ nE = zeros((E.shape[0] +9, E.shape[1] ))
    # ~ ne = zeros(E.shape[0] +9)
    # ~ idRow = E.shape[0]
    # ~ nE[:idRow,:E.shape[1]] = E  
    # ~ ne[:idRow] = e[:]
    # ~ nE[idRow:idRow+3,:DEFAULT_NUM_VARS] = footExpressionMatrix[0][:]  
    # ~ ne[idRow:idRow+3] = posInit
    # ~ idRow += 3
    # ~ nE[idRow:idRow+3,:DEFAULT_NUM_VARS] = footExpressionMatrix[1][:]  
    # ~ nE[idRow:idRow+3,:DEFAULT_NUM_VARS] = footExpressionMatrix[0][:]  
    # ~ ne[idRow:idRow+3] = posInit 
    # ~ ne[idRow:idRow+3] = posInit + array([0., -.20, 0.0])
    # ~ nVarEnd = numVariablesForPhase(pb["phaseData"][-1])
    # ~ print ("nVarEnd", nVarEnd, E.shape[1]-nVarEnd, E.shape[1]-nVarEnd,E.shape[1]-nVarEnd+DEFAULT_NUM_VARS )
    # ~ nE[-3:,E.shape[1]-nVarEnd:E.shape[1]-nVarEnd+DEFAULT_NUM_VARS] = footExpressionMatrix[0][:]   
    # ~ ne[-3:] = posEnd
    # ~ return nE, ne
    
    
    
####################### MAIN ###################"

# try to import mixed integer solver
MIP_OK = False  
try:
    import gurobipy
    import gurobipy as grb
    import cvxpy as cp
    from scipy.sparse import csr_matrix
    MIP_OK = True

except ImportError:
    pass

def tovals(variables):
    return array([el.value for el in variables])
    
def solveMIP(pb, surfaces, MIP = True, draw_scene = None, plot = True):  
    if not MIP_OK:
        print("Mixed integer formulation requires gurobi packaged in cvxpy")
        raise ImportError
        
    gurobipy.setParam('LogFile', '')
    gurobipy.setParam('OutputFlag',0)
       
    A, b, E, e = convertProblemToLp(pb, True)   
    E,e = addInitEndConstraint(pb, E, e)
    slackMatrix = wSelectionMatrix(pb)
    surfaceSlackMatrix = alphaSelectionMatrix(pb)
    
    rdim = A.shape[1]
    varReal = cp.Variable(rdim)
    constraints = []
    constraintNormalIneq = A * varReal <= b
    constraintNormalEq   = E * varReal == e
    
    constraints = [constraintNormalIneq, constraintNormalEq]
    #creating boolean vars
    
    slackIndices = [i for i,el in enumerate (slackMatrix) if el > 0]
    slackIndicesSurf = [i for i,el in enumerate (surfaceSlackMatrix) if el > 0]
    numSlackVariables = len([el for el in slackMatrix if el > 0])
    numSlackVariablesSurf = len([el for el in surfaceSlackMatrix if el > 0])
    boolvars = cp.Variable(numSlackVariables, boolean=True)      
    obj = cp.Minimize(slackMatrix * varReal)
    
    if MIP:    
        # ~ constraints = constraints + [varReal[el] <= boolvars[i] for i, el in enumerate(slackIndices)]   
        constraints = constraints + [varReal[el] == boolvars[i] for i, el in enumerate(slackIndices)]   
        
        
        currentSum = []
        previousL = 0
        for i, el in enumerate(slackIndices):
            if i!= 0 and el - previousL > 1.:
                assert len(currentSum) > 0
                constraints = constraints + [sum(currentSum) == 1 ]
                currentSum = [boolvars[i]]
            elif el !=0:
                currentSum = currentSum + [boolvars[i]]
            previousL  = el
        if len(currentSum) > 1:
            constraints = constraints + [sum(currentSum) == 1 ]
        
        
        if numSlackVariablesSurf > 0:
            boolvarsSurf = cp.Variable(numSlackVariablesSurf, boolean=True)    
            constraints = constraints + [varReal[el] <= boolvarsSurf[i] for i, el in enumerate(slackIndicesSurf)] 
            currentSum = []
            previousL = 0
            for i, el in enumerate(slackIndicesSurf):
                if i!= 0 and el - previousL > 1.:
                    assert len(currentSum) > 0
                    constraints = constraints + [sum(currentSum) <= len(currentSum) -1 ]
                    currentSum = [boolvarsSurf[i]]
                elif el !=0:
                    currentSum = currentSum + [boolvarsSurf[i]]
                previousL  = el
            if len(currentSum) > 1:
                constraints = constraints + [sum(currentSum) <= len(currentSum) -1 ]
        
        obj = cp.Minimize(surfaceSlackMatrix * varReal)
        # ~ obj = cp.Minimize(ones(numSlackVariables) * boolvars)
    prob = cp.Problem(obj, constraints)
    t1 = clock()
    res = prob.solve(solver=cp.GUROBI, verbose=True )
    t2 = clock()
    res = tovals(varReal)
    print("time to solve MIP ", timMs(t1,t2))

    
    # return timMs(t1,t2)
    return pb, res, timMs(t1,t2)


# gurobi cost functions
def targetCom(pb, cVars, endCom):
    nVarEnd = numVariablesForPhase(pb["phaseData"][-1])
    cx_end_diff = cVars[-nVarEnd]   - endCom[0]
    cy_end_diff = cVars[-nVarEnd+1] - endCom[1]
    # ~ cz_end_diff = cVars[-nVarEnd+2] - endCom[2]
    #
    # ~ return cx_end_diff * cx_end_diff + cy_end_diff * cy_end_diff + cz_end_diff * cz_end_diff
    return cx_end_diff * cx_end_diff + cy_end_diff * cy_end_diff
    
# gurobi cost functions
def targetLegCenter(pb, cVars, endCom):
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    #~ fixedFootMatrix = None;
    footOffset = 4
    cost = 0
    
    fact = 1. / float(N_EFFECTORS)
    
    for phaseId, phaseDataT in enumerate(pb["phaseData"]):  
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        footCol = startCol + footOffset
        # expressing everything in terms of first foot
        x= 0; y = 0;
        for footId in range(1, N_EFFECTORS):
            x += cVars[footCol+footId * 3    ] * fact
            y += cVars[footCol+footId * 3 + 1] * fact
        startCol = endCol         
        if phaseId == len(pb["phaseData"]) - 1:
            cost = (x - endCom[0]) * (x - endCom[0]) + (y - endCom[1]) * (y - endCom[1])
    print ("cost ", cost)
    return cost

def posturalCost(pb, cVars, initPos, initCom):
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    #~ fixedFootMatrix = None;
    footOffset = 4
    
    refCost = [ initPos[i] - initPos[0] for i in range(1, N_EFFECTORS) ]
    cost = 0
    
    for phaseId, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        footCol = startCol + footOffset
        # expressing everything in terms of first foot
        for footId in range(1, N_EFFECTORS):
            cost += sum([(cVars[footCol+footId * 3 + el] - cVars[footCol + el] - refCost[footId-1][el])* (cVars[footCol+footId * 3 + el] - cVars[footCol + el] - refCost[footId-1][el])  for el in range(3)])
        startCol = endCol
    return cost
    
def stepSizeCost(pb, cVars, initPos, initCom):
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    #~ fixedFootMatrix = None;
    footOffset = 4
    
    cost = 0
    
    for phaseId, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        footCol = startCol + footOffset
        # expressing everything in terms of first foot
        for footId in range(1, N_EFFECTORS):
            cost += sum([(cVars[footCol+footId * 3 + el] - 0.15)* (cVars[footCol+footId * 3 + el] - 0.15)  for el in range(2)])
        startCol = endCol
    return cost
    
def maxStepSizeCost(pb, cVars, initPos, initCom):
    startCol = 0;
    previousStartCol = 0;
    endCol   = 0;
    #~ fixedFootMatrix = None;
    footOffset = 4
    
    cost = 0
    
    for phaseId, phaseDataT in enumerate(pb["phaseData"]):   
        #inequality
        endCol = startCol + numVariablesForPhase(phaseDataT)
        footCol = startCol + footOffset
        # expressing everything in terms of first foot
        for footId in range(1, N_EFFECTORS):
            cost += sum([(cVars[footCol+footId * 3 + el] - 0.15)* (cVars[footCol+footId * 3 + el] - 0.15)  for el in range(2)])
        startCol = endCol
    return cost


def solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, initGuess = None, initGuessMip = None, l1Contact = False, initPos = None,  endPos = None, initCom = None,  endCom = None):  
    if not MIP_OK:
        print ("Mixed integer formulation requires gurobi packaged in cvxpy")
        raise ImportError
        
    grb.setParam('LogFile', '')
    grb.setParam('OutputFlag', 0)
    
    grb.setParam('OutputFlag', 1)
       
    A, b, E, e = convertProblemToLp(pb)   
    print ("initpos ", initPos)
    # ~ E,e = addInitEndConstraint(pb, E, e, initPos, endPos, initCom, endCom)
    #todo is end constraint desirable ?
    E,e = addInitEndConstraint(pb, E, e, initPos, endPos, initCom, None)
    slackMatrix = wSelectionMatrix(pb)    
    slackIndices = [i for i,el in enumerate (slackMatrix) if el > 0]
    numSlackVariables = len([el for el in slackMatrix if el > 0])
    
    surfaceSlackMatrix = alphaSelectionMatrix(pb)
    slackIndicesSurf = [i for i,el in enumerate (surfaceSlackMatrix) if el > 0]
    numSlackVariablesSurf = len([el for el in surfaceSlackMatrix if el > 0])
    
    model = grb.Model("mip")
    
    rdim = A.shape[1]
    
    #add continuous variables
    cVars = []
    for i in range(rdim):
        if slackMatrix[i]>0:
            if MIP:
                cVars.append(model.addVar(name="slack%d" % i, obj = 0, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
            else:
                cVars.append(model.addVar(name="slack%d" % i, obj = 1, vtype=grb.GRB.CONTINUOUS, lb = -grb.GRB.INFINITY, ub = grb.GRB.INFINITY))
        elif surfaceSlackMatrix[i]>0:            
            if MIP and not l1Contact:
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
                                 obj=1,
                                 name="boolVar%d" % i))
        model.update()
        
            
        #Big M value
        M = 1000.
        [model.addConstr(cVars[el] == boolvars[i], "boolW%d" % i ) for i, el in enumerate(slackIndices)]  
    
    
               
        #adding constraint that two consecutive contacts can't be for the same effector
        for i in range(N_EFFECTORS):
            for slackId in range(N_EFFECTORS+i, numSlackVariables, N_EFFECTORS):
                print ("id ", slackId, slackId- N_EFFECTORS )
                model.addConstr(boolvars[slackId] + boolvars[slackId - N_EFFECTORS] <= 1)
        
        model.update()
    
        currentSum = []
        currentSum2 = []
        previousL = 0
        for i, el in enumerate(slackIndices):
            if i!= 0 and el - previousL > 1.:
                # ~ model.addConstr(grb.quicksum(currentSum) <= len(currentSum) -1, "card%d" % i)
                model.addConstr(grb.quicksum(currentSum) == 1, "card%d" % i)
                assert len(currentSum) > 0      
                currentSum = [boolvars[i]]
                currentSum2 = [el]          
            elif el !=0:
                currentSum = currentSum + [boolvars[i]]
                currentSum2 = currentSum2 + [el]
            previousL  = el        
        if len(currentSum) > 1:
            model.addConstr(grb.quicksum(currentSum) == 1, "card%d" % i)
        
        # ~ model.addConstr(grb.quicksum(boolvars) == 4, "card%d" % i)
        print ("bools ", len(boolvars))
         
        if not l1Contact:        
            boolvarsSurf = []
            for i in range(numSlackVariablesSurf):
                boolvarsSurf.append(model.addVar(vtype=grb.GRB.BINARY,
                                     obj=0,
                                     name="boolvarSurf%d" % i))
                
            [model.addConstr(cVars[el] <= boolvarsSurf[i], "boolAlpha%d" % i ) for i, el in enumerate(slackIndicesSurf)]  
            currentSum = []
            currentSum2 = []
            previousL = 0
            for i, el in enumerate(slackIndicesSurf):
                if i!= 0 and el - previousL > 1.:
                    model.addConstr(grb.quicksum(currentSum) == len(currentSum) -1, "card%d" % i)
                    assert len(currentSum) > 0      
                    currentSum = [boolvarsSurf[i]]
                    currentSum2 = [el]          
                elif el !=0:
                    currentSum = currentSum + [boolvarsSurf[i]]
                    currentSum2 = currentSum2 + [el]
                previousL  = el        
            if len(currentSum) > 1:
                model.addConstr(grb.quicksum(currentSum) == len(currentSum) -1, "card%d" % i)
            
            
    model.modelSense = grb.GRB.MINIMIZE
    # ~ model.Params.SOLUTION_LIMIT = 1
    # ~ model.Params.TIME_LIMIT = 10
    
    if initGuess is not None:
        for (i,el) in initGuess:
            x[i].start = el
    
    if MIP and initGuessMip is not None:
        for (i,el) in initGuessMip:
            boolvars[i].start = el
        
    
    if initPos is not None:
        obj = 0
        print (" initPos is not None !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        obj += 0.01 * posturalCost(pb, cVars,  initPos, initCom)
        # ~ obj += 0.01 * posturalCost(pb, cVars,  initPos, initCom)
        # ~ obj += stepSizeCost(pb, cVars,  initPos, initCom)
        # ~ obj += 10. * targetCom(pb, cVars,  endCom)
        obj += 10. * targetCom(pb, cVars,  endCom)
        # ~ obj = targetLegCenter(pb, cVars,  endCom)
        model.setObjective(obj)
    
    # ~ grb.setParam('SOLUTION_LIMIT', 1)
    # ~ grb.setParam('TIME_LIMIT', 10)
    model.update()
    t1 = clock()
    model.optimize()
    t2 = clock()
    res = [el.x for el in cVars]
    resbool = [el.x for el in boolvars]
    print ("resbool", resbool)
    print ("time to solve MIP ", timMs(t1,t2))
    
    
    #~ return timMs(t1,t2)
    
    return pb, res, timMs(t1,t2)
    # ~ if MIP:
        # ~ return res
    # ~ else:
        # ~ return res


if __name__ == '__main__':
    from sl1m.stand_alone_scenarios.escaliers import gen_stair_pb,  draw_scene, surfaces
    # ~ from sl1m.stand_alone_scenarios.complex import gen_stair_pb,  draw_scene, surfaces
    pb = gen_stair_pb()    
    
    t1 = clock()
    initGlobals(nEffectors = 2)
    # ~ A, b, E, e = convertProblemToLp(pb)
    # ~ E,e = addInitEndConstraint(pb, E, e)
    
    # ~ pb, res, time = solveMIP(pb, surfaces, MIP = True, draw_scene = None, plot = True)
    
   
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = None)
    # ~ pb, res, time = solveMIPGurobi(pb, surfaces, MIP = False, draw_scene = None, plot = True, l1Contact = True)
    # ~ pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = True)
    # ~ pb, res, time = solveMIPGurobi(pb, surfaces, MIP = False, draw_scene = None, plot = True)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False)
    # ~ pb = gen_stair_pb()  
    # ~ pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = True)
    # ~ ax = draw_scene(None)
    # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    # ~ C = identity(A.shape[1]) * 0.00001
    # ~ c =  wSelectionMatrix(pb) * 100.
    # ~ t2 = clock()
    # ~ res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    # ~ t3 = clock()
    
    # ~ print("time to set up problem" , timMs(t1,t2))
    # ~ print("time to solve problem"  , timMs(t2,t3))
    # ~ print("total time"             , timMs(t1,t3))
    
    # ~ if res.success:
        # ~ print ("success")
    # ~ else:
        # ~ print ('problem infeasible')
        # ~ assert False
    
    w = wSelectionMatrix(pb)
    b = bSelectionMatrix(pb)
    g = gSelectionMatrix(pb)
    al = alphaSelectionMatrix(pb)
    
    wR = [ el for el in w *res if abs(el) > 0.00001]
    bR = [ el for el in b *res if abs(el) > 0.0001]
    gR = [ el for el in g *res if abs(el) > 0.0001]
    alR = [ el for el in al *res if abs(el) > 0.00000001]
    
    print ("w", wR)
    print ("b", bR)
    print ("g", gR)
    print ("al", alR)
    
    # ~ ax = draw_scene(None)
    # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    # ~ surfaces, indices = bestSelectedSurfaces(pb, res)
    # ~ for i, phase in enumerate(pb["phaseData"]):  
        # ~ phase["S"] = [surfaces[i]]
        
    # ~ t1 = clock()
    # ~ A, b, E, e = convertProblemToLp(pb, False)
    
    
    # ~ C = identity(A.shape[1])
    # ~ c = zeros(A.shape[1])
    # ~ t2 = clock()
    # ~ res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    # ~ t3 = clock()
    
    # ~ print("time to set up problem" , timMs(t1,t2))
    # ~ print("time to solve problem"  , timMs(t2,t3))
    # ~ print("total time"             , timMs(t1,t3))
    
    # ~ coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    # ~ ax = draw_scene(None)
    # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False)
        
    
    
    
    
        
