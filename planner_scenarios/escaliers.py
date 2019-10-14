import numpy as np


from mpcroc.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *
#~ from qp import solve_lp

from mpcroc.planner import *


from plot_plytopes import *

floor = [ [0.16, 1., 0.], [-1.8, 1., 0.], [-1.8, -1., 0.], [0.16, -1., 0.]  ]

step1 = [[0.3, 0.6, 0.15],[0.3, -0.16, 0.15],[0.6, -0.16, 0.15],[0.6, 0.6, 0.15]]
step2 = [[0.6, 0.6, 0.3 ],[0.6, -0.16, 0.3 ],[0.9, -0.16, 0.3 ],[0.9, 0.6, 0.3 ]]
step3 = [[0.9, 0.6, 0.45],[0.9, -0.16, 0.45],[1.2, -0.16, 0.45],[1.2, 0.6, 0.45]]
step4 = [[1.2, 0.6, 0.6 ],[1.2, -0.16, 0.6 ],[1.5, -0.16, 0.6 ],[1.5, 0.6, 0.6 ]]


floor = [ [0.16, 1., 0.], [-1.8, 1., 0.], [-1.8, -1., 0.], [0.16, -1., 0.]  ]

step1 = [[0.3, 0.6, 0.1],[0.3, -0.16, 0.1],[0.6, -0.16, 0.1],[0.6, 0.6, 0.1]]
step2 = [[0.6, 0.6, 0.2 ],[0.6, -0.16, 0.2],[0.9, -0.16, 0.2 ],[0.9, 0.6, 0.2 ]]
step3 = [[0.9, 0.6, 0.3],[0.9, -0.16, 0.3],[1.2, -0.16, 0.3],[1.2, 0.6, 0.3]]
step4 = [[1.2, 0.6, 0.4 ],[1.2, -0.16, 0.4 ],[1.5, -0.16, 0.4 ],[1.5, 0.6, 0.4 ]]


all_surfaces = [floor, step1, step2, step3, step4]

afloor = array(floor).T 
astep1 = array(step1).T
astep2 = array(step2).T
astep3 = array(step3).T
astep4 = array(step4).T

surfaces = [[afloor], [afloor], [astep1,astep2,astep3],[astep2,astep3,astep1], [astep3,astep2,astep1,astep4], [astep3,astep4], [astep4],[astep4]]

def gen_stair_pb():
    #~ kinematicConstraints = genKinematicConstraints(min_height = 0.6)
    kinematicConstraints = genKinematicConstraints(min_height = None)
    relativeConstraints = genFootRelativeConstraints()
    nphases = len(surfaces)
    p0 = None
    p0 = [array([0.,0., 0.]), array([0.,0., 0.])];
    res = { "p0" : p0, "c0" : None, "nphases": nphases}
    
    #TODO in non planar cases, K must be rotated
    phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i]))], "relativeK" : [relativeConstraints[(i) % 2] for _ in range(len(surfaces[i]))], "S" : surfaces[i] } for i in range(nphases)]
    res ["phaseData"] = phaseData
    return res 
    
    
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import scipy as sp

def draw_rectangle(l, ax):
    #~ plotPoints(ax,l)
    lr = l + [l[0]]
    cx = [c[0] for c in lr]
    cy = [c[1] for c in lr]
    cz = [c[2] for c in lr]
    ax.plot(cx, cy, cz)
    
def draw_scene(ax = None, color = "p"):
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    [draw_rectangle(l,ax) for l in all_surfaces]
    return ax
    
    
    
############# main ###################    

if __name__ == '__main__':
    
    from mpcroc.fix_sparsity import solveL1
    
    pb = gen_stair_pb()
    
    pb, coms, footpos, allfeetpos, res = solveL1(pb, surfaces, draw_scene)
    #~ ax = draw_scene()
    #~ plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    
    
