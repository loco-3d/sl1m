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

from scipy.optimize import linprog

def solve_lp(q, G=None, h=None, C=None, d=None): 
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=None, method='interior-point', callback=None, options=None)
    return res    


if __name__ == '__main__':
    #actually this is not used in this prototype, because the constraints are not there
    step1 = [[0.3, 0.6, 0.15],[0.3, -0.16, 0.15],[0.6, -0.16, 0.15],[0.6, 0.6, 0.15]]
    step2 = [[0.6, 0.6, 0.3 ],[0.6, -0.16, 0.3 ],[0.9, -0.16, 0.3 ],[0.9, 0.6, 0.3 ]]
    
    #LP of of dimension 10: find a point that is a convex combination of either the extremum points
    # of step 1 and step2 . This makes 8 variables, we add 2 variables that will be used to handle norm minimization
    
    nvars = 10
    
    #~ dimEq = 2
    dimEq = 1
    Eq = zeros((dimEq,nvars)); eq = ones(dimEq) * 1.2
    #equality constraint, convex combination of all points is 1
    Eq[0,:] = ones(nvars); Eq[0,-2:] = zeros(2);
    #equality constraint, sum of normalization vars is one
    #~ Eq[1,-2:] = ones(2); 
    
    #~ n_cons= nvars +4 + 1
    n_cons= nvars +4 
    A = zeros((n_cons, nvars))
    b = zeros(n_cons)
    
    currentRow = 0
    
    #inequality constraints: all variables are positive
    A[currentRow: nvars, :nvars    ] = -identity(nvars); 
    b[currentRow: nvars            ] = zeros(nvars)
    
    currentRow+= nvars
    
    
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
    #~ x = res["x"]
