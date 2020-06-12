############### Problem definition #############

from sl1m.problem_definition import *
           
        
#### problem related global variables #### 

NUM_INEQ_PER_SLACK = 1 #for all i,  sum_j a_i_j <= b_i
NUM_INEQ_PER_SURFACEPOINT = 1 #a_ij > 0
NUM_EQ_PER_SLACK = 0
NUM_EQ_PER_SURFACE = 0

#### surface
#a surface is a list of points

def numVariablesSurface(s):
    return len(s)
        
def numVariableSurfaces(S):
    return sum([numVariablesSurface(s) for s in S])

def numSlackVariables(S):
    if len(S) > 1:
        return len(S)
    return 0
    
def numVariablesForPhase(phase):
    S = phase["S"]
    return numVariableSurfaces(S) + numSlackVariables(S)

def numIneqForPhase(phase):
    S = phase["S"]
    return NUM_INEQ_PER_SURFACEPOINT * numVariableSurfaces(S) + NUM_INEQ_PER_SLACK * numSlackVariables(S)
    
def numEqForPhase(phase):
    return 1
    
def getTotalNumVariablesAndConstraintSize(pb):
    nphases = pb["nphases"]
    cols = sum([numVariablesForPhase(phase) for phase in  pb["phaseData"]])
    rowsIneq = sum( [numIneqForPhase(phase) for  phase in pb["phaseData"] ])
    rowsEq = sum( [numEqForPhase(phase) for  phase in pb["phaseData"] ])
    return cols, rowsEq,rowsIneq 
    


def slackSelectionMatrix(pb):
    cols, rowsEq,rowsIneq  =  getTotalNumVariablesAndConstraintSize(pb)
    c = zeros(cols)
    col = 0
    for phase in pb["phaseData"]: 
        hasSlack = len(phase["S"]) > 1
        for s in phase["S"]:
            col += numVariablesSurface(s)
            if hasSlack:
                c[col] = 1
                col +=1
    return c



pb = {"nphases" : 1,  "phaseData" : [ {"S" : [[ [1.,1.,0.], [1.,2.,0.], [2.,2.,0.]], [ [4.,2.,0.], [5.,2.,0.], [5.,3.,0.]]]} ] }

cols, rowsEq,rowsIneq  = getTotalNumVariablesAndConstraintSize(pb)


A = zeros((rowsIneq,cols))
b = zeros((rowsIneq))

E = zeros((rowsEq,cols))
e = zeros((rowsEq))

phase = pb["phaseData"][0]
S = pb["phaseData"][0]["S"]
rowIneq = 0
rowEq = 0
col = 0
endCol = 0
summationRow = 0
hasSlackVars = numSlackVariables(S) > 0
for s in S:
    print ("rowineq 0 ", rowIneq)
    #first equality vars then slacks
    nVarSurf = numVariablesSurface(s)
    print ("endCol ", endCol)
    endCol = col + nVarSurf
    #positivity constraints
    A[rowIneq:rowIneq+nVarSurf, col:endCol] = -identity(nVarSurf); 
    rowIneq += nVarSurf
    print ("rowineq ", rowIneq)
    #always updating the same row to say that the sum is equal to one
    E[rowEq,col:col+nVarSurf] = ones(nVarSurf)
    if hasSlackVars:
        A[rowIneq, col:endCol] = ones(nVarSurf)
        A[rowIneq, endCol] = -1
        endCol += 1
        rowIneq +=1
    col = endCol
E[0,:] = ones(E.shape[0]) - E[0]
e[rowEq]=1



c = slackSelectionMatrix(pb)

filt = ones(c.shape) - c

x_x = [s_in[0] for s in S for s_in in s]
x_y = [s_in[1] for s in S for s_in in s]
x_z = [s_in[2] for s in S for s_in in s]

x_val = zeros(c.shape)
y_val = zeros(c.shape)
z_val = zeros(c.shape)

idx = 0
for i, el in enumerate(filt):
    if el > 0.5:
        x_val[i] = x_x[idx]
        y_val[i] = x_y[idx]
        z_val[i] = x_z[idx]
        idx+=1

import sl1m.qp

unconst = qp.solve_lp_glpk(c,A,b,E,e)

#add random constraint
# x > 2.5
A = vstack([A,-x_val])
b = vstack([b.reshape([-1,1]),array([[-1.5]])]).reshape((-1,))
A = vstack([A,x_val])
b = vstack([b.reshape([-1,1]),array([[3.2]])]).reshape((-1,))

const = qp.solve_lp_glpk(c,A,b,E,e)
# ~ const = qp.solve_lp_glpk(c,A,b)
res = const.x

EPS_ = 0.01
#update cost function based on previous iteration
def reweight(x_i, c, eps = 0.01):
    return (1. / (x_i + ones(x_i.shape[0]) * eps)) * c
    
# ~ i = 0
# ~ MAX_ITER = 1
# ~ while i < MAX_ITER:
    # ~ i +=1
    # ~ c = reweight(array(res), c)
    # ~ res = qp.solve_lp_glpk(c,A,b,E,e).x
