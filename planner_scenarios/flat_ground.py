import numpy as np


from mpcroc.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *
#~ from qp import solve_lp

from mpcroc.planner import *


from mpcroc.plot_plytopes import *


all_surfaces = []

left = [[-50.,0.04,0.], [50.,0.04,0.], [50.,40.,0.], [-50.,40.,0.]]
right = [[-50.,-40,0.], [50.,-40.,0.], [50.,-0.04,0.], [-50.,-0.04,0.]]
allsurf = [[-50.,-40,0.], [50.,-40.,0.], [50.,40,0.], [-50.,40,0.]]

#~ all_surfaces = [left, right]
all_surfaces = [allsurf]

# generate problem for flat ground scenarios, from init pos for feet p0LF and p0RF
# this method also defines the template definition for defining contact plan
def genFlatGroundProblem(p0LF, p0RF, nphases = 5):
    #~ for i in range(nphases)
    #~ kinematicConstraints = genKinematicConstraints(min_height = 0.4)
    kinematicConstraints = genKinematicConstraints()
    relativeConstraints = genFootRelativeConstraints()
    surfaces = [None, None];
    p0 = [[0., -0.06,0.], [0., 0.06,0.]];
    #~ surfaces[LF] = array(right).T  # for the moment arbitrary surfaces to separate foot
    surfaces[LF] = [array(allsurf).T]  # for the moment arbitrary surfaces to separate foot
    #~ surfaces[LF] = array([[-50.,-40,0.], [50.,-40.,0.], [50.,-0.04,0.], [-50.,-0.04,0.] ]).T
    surfaces[RF] = [array(allsurf).T]
    p0 [LF]  = p0LF;
    p0 [RF]  = p0RF;
    
    res = { "p0" : p0, "c0" : array([0., -0.06,0.5]), "nphases": nphases}
    #~ res = { "p0" : None, "nphases": nphases}
    
    #TODO in non planar cases, K must be rotated
    #relative k is the constraint with respect to the previous frame    
    phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i%2]))], "relativeK" : [relativeConstraints[(i) % 2] for _ in range(len(surfaces[i%2]))], "S" : surfaces[i%2] } for i in range(nphases)]
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
    pb = genFlatGroundProblem([0.,0.05,0.],[0.,-0.05,0.], 20)
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    
    print 'A.shape', A.shape
    
    C, c = least_square_cost_function (pb, A.shape[1], comTarget = array([3.,10.,0.2]))
    t2 = clock()
    res = qp.solve_least_square(C,c,A,b,E,e)
    t3 = clock()
    
    
    print "time to set up problem" , timMs(t1,t2)
    print "time to solve problem"  , timMs(t2,t3)
    print "total time"             , timMs(t1,t3)
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    plotQPRes(pb, res, plot_constraints = False)
    #~ plotQPRes(pb, res, plot_constraints = True)
    
    
