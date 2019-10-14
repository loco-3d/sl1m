import numpy as np


#~ from constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *

from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray


import cvxpy as cp
from numpy.random import randn
from numpy import eye
eps =0.000001


from scipy.optimize import linprog

def solve_lp(q, G=None, h=None, C=None, d=None): 
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=None, method='interior-point', callback=None, options=None)
    return res    

def solve_qcqp(q, A=None, b=None, Eq=None, eq=None, Hsq = None, Asq = None, bsq = None):
    nvars = q.shape[0]
    x  = cp.Variable(nvars)
    constraints = []
    print "wtf", A * x
    
    if A is not None:
        constraints = [A * x <= b]    
    if Eq is not None:
        constraints = constraints + [Eq * x  == eq] 
    
    #~ if Hsq is not None:
        #~ if Asq is not None:
            #~ constraints = constraints + [ cp.quad_form(x,Hsq) + Asq * x <= bsq] 
        #~ else:
            #~ constraints = constraints + [ cp.quad_form(x,Hsq)           <= bsq] 
    
    cost = cp.quad_form(x,Hsq) + Asq*x
    
    #~ obj = cp.Minimize(q * x)
    obj = cp.Minimize(cost)
    prob = cp.Problem(obj, constraints)
    prob.solve()
    return prob, x
    
if __name__ == '__main__':
    #actually this is not used in this prototype, because the constraints are not there
    step1 = [[0.3, 0.6, 0.15],[0.3, -0.16, 0.15],[0.6, -0.16, 0.15],[0.6, 0.6, 0.15]]
    step2 = [[0.6, 0.6, 0.3 ],[0.6, -0.16, 0.3 ],[0.9, -0.16, 0.3 ],[0.9, 0.6, 0.3 ]]
    
    #LP of of dimension 10: find a point that is a convex combination of either the extremum points
    # of step 1 and step2 . This makes 8 variables, we add 2 variables that will be used to handle norm minimization
    
    nvars = 10
    # 8  + norm1(z) + z
    
    #~ dimEq = 2
    dimEq = 2
    Eq = zeros((dimEq,nvars)); eq = ones(dimEq) 
    #equality constraint, convex combination of all points is 1
    Eq[0,:] = ones(nvars); Eq[0,-2:] = zeros(2);
    #equality constraint, sum of normalization vars is one
    Eq[1,-2:] = 1; 
    Eq[,:4] = ones(4) * 2.; 
    # z = 1 - 2x
    
    #~ n_cons= nvars +4 + 1
    n_cons= nvars +4 
    A = zeros((n_cons, nvars))
    b = zeros(n_cons)
    
    currentRow = 0
    
    #inequality constraints: all variables are positive
    A[currentRow: nvars-1, :nvars-1    ] = -identity(nvars-1); 
    b[currentRow: nvars-1              ] = zeros(nvars-1)
    
    currentRow+= nvars - 1
    
    
    #norm1 constraints, norm of convex combination of step1 is variable 8
    # sum(x(i) * step1(i)) <= x [8]
    A[currentRow, :4 ] = ones(4);
    A[currentRow,  8 ] = -1.;
    currentRow+= 1
    # sum(x(i) * step1(i)) >= - x [8]
    A[currentRow, :4 ] = -ones(4);
    A[currentRow,  8 ] = -1.;
    currentRow+= 1
    
    
    #norm1 constraints, norm of convex combination of step2 is variable 9
    # sum(x(i+4) * step2(i)) <= x [9]
    A[currentRow, 4:8 ] = ones(4);
    A[currentRow, 9 ] = -1.;
    currentRow+= 1
    # sum(x(i+4) * step2(i)) >= - x [9]
    A[currentRow, 4:8 ] = -ones(4);
    A[currentRow, 9 ] = -1.;
    currentRow+= 1
    
    #~ A[currentRow, :8 ] = -ones(8);
    #~ b[currentRow ]     = -1;
    #~ currentRow+= 1
    
    assert currentRow == n_cons
    
    #cost is minimize l1 norm of variable sum, that is 8 and 9
    #~ c = zeros(nvars); c[-2:] = ones(2)
    #~ c = zeros(nvars); c[-2:] = ones(2); 
    c = zeros(nvars); c[-2:] = ones(2) * 0.001; 
    #~ c[-1] = -1
    #~ c = ones(nvars);c[-2:] = zeros(2)
    #~ c = zeros(nvars); c[-2:] = -ones(2); c[-1] = -c[-1]
    res = solve_lp(c,A,b,Eq,eq)
    #~ res = solve_lp(c,A,b)
    x = res["x"]
    
    print "first step variables", x[:4]
    print "second step variables", x[4:8]
    print "normalization variables", x[8:]
    
    
    #~ Eq = hstack([ones(4),zeros(4)]).reshape((1,-1)); eq = ones(1)
    #~ A = zeros((8,8)); b = zeros(8)
    #~ for i in range(0,8,2):
        #~ print ' 4+ i%2',  (i%2)
        #~ A[i  , i/2 ] = 1. ; A[i, 4+ i/2 ] = -1.;
        #~ A[i+1, i/2 ] = -1.; A[i, 4+ i/2 ] = -1.;
        
    #~ c = hstack([zeros(4),ones(4)]); 
    #~ res = solve_lp(c,A,b,Eq,eq)
    #~ x = res["xb
    
    #add l1 norm var
    #~ A = hstack([A,zeros((A.shape[0],2))])
    #~ lastPos =  zeros(A.shape[1]); 
    #~ lastPos[-2] = -1.
    #~ A = vstack([A,lastPos])
    #~ lastPos =  zeros(A.shape[1]); 
    #~ lastPos[-1] = -1.
    #~ A = vstack([A,lastPos])
    
    #~ b = concatenate([b,zeros(2)])
    
    
    #~ Eq = hstack([Eq,zeros((Eq.shape[0],2))])
    #~ lastEq =  zeros(Eq.shape[1]); 
    #~ lastPos[-2:] = ones(2)
    #~ eq = concatenate([eq,ones(1)])
    
    q = zeros(c.shape[0] + 1)
    q[-1] = -1
    
    Hsq = zeros((q.shape[0],q.shape[0]))
    Hsq[8:-1, 8:-1] = identity(2)
    #~ Hsq[-1,-1] = 0.
    Asq = zeros(q.shape[0])
    Asq[-1] = 1.
    bsq = 1.
    
    
    hsqp = zeros((c.shape[0],c.shape[0]))
    hsqp[8:, 8:] = identity(2)
    #~ hsqp[8:] = ones((2,1))
    
    bsqp = zeros(c.shape[0])
    
    #~ H, h = lsq(hsqp.T, bsqp)
    
    print "H ", hsqp
    print "h ", bsqp
    
    
    
    #~ pb, x2 = solve_qcqp(q=q, A=A, b=b, Eq=Eq, eq=eq, Hsq = Hsq, Asq = Asq, bsq = bsq)
    pb, x2 = solve_qcqp(q=c, A=A, b=b, Eq=Eq, eq=eq, Hsq = hsqp, Asq = bsqp, bsq = bsq)
