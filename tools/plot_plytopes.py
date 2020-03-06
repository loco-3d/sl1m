import numpy as np
from hpp_centroidal_dynamics import *
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, vstack, hstack, asmatrix, identity, cross
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *
#~ from qp import solve_lp

import eigenpy
import cdd
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray

eigenpy.switchToNumpyArray()



from sl1m.constants_and_tools import *



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



Id = matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
z = array([0.,0.,1.])
zero3 = zeros(3) 


def generators(A,b, Aeq = None, beq = None):
    m = np.hstack([b,-A])
    matcdd = cdd.Matrix(m); matcdd.rep_type = cdd.RepType.INEQUALITY
    
    if Aeq is not None:
        meq = np.hstack([beq,-Aeq])
        matcdd.extend(meq.tolist(), True)
    
    H = cdd.Polyhedron(matcdd)
    g = H.get_generators()
    
    return [array(g[el][1:]) for el in range(g.row_size)], H
    
def filter(pts):
    hull = ConvexHull(pts, qhull_options='Q12')
    return [pts[i] for i in hull.vertices.tolist()]
    
def ineq(pts, canonicalize = False):
    apts = array(pts)
    m = np.hstack([ones((apts.shape[0],1)),apts])
    matcdd = cdd.Matrix(m); matcdd.rep_type = cdd.RepType.GENERATOR
    H = cdd.Polyhedron(matcdd)
    bmA = H.get_inequalities()
    if canonicalize:
        bmA.canonicalize()
    Ares = zeros((bmA.row_size,bmA.col_size-1))
    bres = zeros(bmA.row_size )
    for i in range(bmA.row_size):
        l = array(bmA[i])
        Ares[i,:] = -l[1:]
        bres[i]   =  l[0]
    return Ares, bres
    
def ineqQHull(hull):
    A = hull.equations[:,:-1]
    b = -hull.equations[:,-1]
    return A,b
    
    
def canon(A,b):
    m = np.hstack([b,-A])
    matcdd = cdd.Matrix(m); matcdd.rep_type = 1
    H = cdd.Polyhedron(matcdd)
    bmA = H.get_inequalities()
    #~ bmA.canonicalize()
    Ares = zeros((bmA.row_size,bmA.col_size-1))
    bres = zeros((bmA.row_size,1 ))
    for i in range(bmA.row_size):
        #~ print "line ", array(bmA[i])
        #~ print "A ", A[i][:]
        #~ print "b ", b[i]
        l = array(bmA[i])
        Ares[i,:] = -l[1:]
        bres[i]   =  l[0]
        #~ print "Ares ",Ares[i,:]
        #~ print "bres ",bres[i]
    return Ares, bres

def genPolytope(A,b):
    pts, H = generators(A,b)
    apts = array(pts)
    if len(apts) > 0:
        hull = ConvexHull(apts)
        return hull, pts, apts, H
    return None, None, None, None

off = 0

def plot_hull_in_subplot(hull, pts, apts, ax, color = "r", just_pts = False):
    global off
    off = (off + 1) % 5
    # Plot defining corner points
    #~ ax.plot(apts.T[0], apts.T[1], apts.T[2], "ko")
    if not just_pts:
        for s in hull.simplices:
            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
            if apts[0].shape[0] < 3:
                ax.plot(apts[s, 0], apts[s, 1], -off*0.05, color+"-")
            else:
                ax.plot(apts[s, 0], apts[s, 1], apts[s, 2], color+"-")


def plot_hull(hull, pts, apts, color = "r", just_pts = False, ax = None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    plot_hull_in_subplot(hull, pts, array(pts), ax, color, just_pts)
    #~ plt.show()

def plot_polytope_H_rep(A_in,b_in, color = "r", just_pts = False, ax = None):
    hull, pts, apts, cd = genPolytope(A_in,b_in)
    if hull is None:
        print ("empty polytope", )
        return False
    plot_hull(hull, pts, apts, color, just_pts, ax = ax)
    return True

def plot_polytope_V_rep(pts, color = "r", just_pts = False, ax = None):
    pts = [array(el) for el in pts]
    apts = array(pts)
    hull = ConvexHull(apts, qhull_options='Q12')
    plot_hull(hull, pts, apts, color, just_pts, ax = ax)
    
def plot_polytope_CDD_PolyHeron(H, color = "r", just_pts = False):
    g = H.get_generators()
    pts = [array(g[el][1:]) for el in range(g.row_size)]
    plot_polytope_V_rep(pts, color, just_pts)
    
