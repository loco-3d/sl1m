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

def solve_lp(c, A=None, b=None, Eq=None, eq=None): 
    res = linprog(c, A_ub=A, b_ub=b, A_eq=Eq, b_eq=eq, bounds=None, method='interior-point', callback=None, options=None)
    return res    

    
if __name__ == '__main__':
    ## vars: x, |b|, |c|
    nvars = 3
    
    #defaut lp:
    #find, x, b, c
    # min |b| + |c|
    # st x <= -0.1 + b
    # st x >=  0.1 - c
    
    n_cons = 2 + 2 #c, b >= 0, x ineq
    A = zeros((n_cons,nvars))
    b = zeros(n_cons)
    A[:2,1:] = -identity(2) #b, c >= 0
    A[2,:] = [1,-1,0] ; b[2] = -0.1  # st x  - b <= -0.1 
    A[3,:] = [-1,0,1];  b[3] = -0.1  # st -x - c <= -0.1 
    
    c = array([0,1.,1.])
    res =  solve_lp(c, A, b)
    
    print "x , b, c, for nominal lp: [", res["x"][0] ,", " , res["x"][1], ", " , res["x"][2],"]"
    
    
    ## vars: x, y, |x|, |y|
    nvars = 4
    
    #defaut lp:
    #find, x, y
    # min |x| + |y|
    # st x + y = 1 
    
    # x + y = 1
    Eq = array([[1,1,0.,0.]])
    eq = array([1.])
    #~ Eq = None; eq = None
    
    n_cons = 2 + 2 #x, y >= 0, + norm constraint
    A = zeros((n_cons,nvars))
    b = zeros(n_cons)
    A[:2,:2] = -identity(2) #x, y > 0
    A[2,:] = [1,0,-1,0] 
    A[3,:] = [0,1,0,-1]
    
    c = array([0,0,1.,1.])
    res =  solve_lp(c, A, b, Eq, eq)
    
    print "x, y for nominal lp: [", res["x"][0] ,", " , res["x"][1], "]"
    
