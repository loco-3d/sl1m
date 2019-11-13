import numpy as np
from hpp_centroidal_dynamics import *
from hpp_spline import *
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, vstack, hstack, asmatrix, identity, cross
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *
#~ from qp import solve_lp

import eigenpy
import cdd
#from curves import bezier3
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray

eigenpy.switchToNumpyArray()



from sl1m.constants_and_tools import *



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_hull_in_subplot(hull, pts, apts, ax, color = "r", just_pts = False):
    # Plot defining corner points
    #~ ax.plot(apts.T[0], apts.T[1], apts.T[2], "ko")
    if not just_pts:
        for s in hull.simplices:
            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
            ax.plot(apts[s, 0], apts[s, 1], apts[s, 2], color+"-")


def plot_hull(hull, pts, apts, color = "r", just_pts = False, ax = None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    plot_hull_in_subplot(hull, pts, array(pts), ax, color, just_pts)
    #~ plt.show()

def plot_polytope_H_rep(A_in,b_in, color = "r", just_pts = False, ax = None):
    hull, pts, apts, cd = genPolytope(A_in,b_in)
    plot_hull(hull, pts, apts, color, just_pts, ax = ax)

def plot_polytope_V_rep(pts, color = "r", just_pts = False, ax = None):
    pts = [array(el) for el in pts]
    apts = array(pts)
    hull = ConvexHull(apts, qhull_options='Q12')
    plot_hull(hull, pts, apts, color, just_pts, ax = ax)
    
def plot_polytope_CDD_PolyHeron(H, color = "r", just_pts = False):
    g = H.get_generators()
    pts = [array(g[el][1:]) for el in range(g.row_size)]
    plot_polytope_V_rep(pts, color, just_pts)
    
