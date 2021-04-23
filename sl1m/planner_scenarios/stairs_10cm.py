import numpy as np


from sl1m.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *
#~ from qp import solve_lp

from sl1m.planner import *


floor = [[-0.07, 0.3, 0.], [0.23, 0.3, 0. ], [0.23, 1.3, 0.  ], [-0.07, 1.3, 0. ], ]
step1 = [[0.23, 1.3, 0.1], [0.53, 1.3, 0.1], [0.53, 0.3, 0.1 ], [0.23, 0.3, 0.1 ], ]
step2 = [[0.53, 1.3, 0.2], [0.83, 1.3, 0.2], [0.83, 0.3, 0.2 ], [0.53, 0.3, 0.2 ], ]
step3 = [[0.83, 1.3, 0.3], [1.13, 1.3, 0.3], [1.13, 0.3, 0.3 ], [0.83, 0.3, 0.3 ], ]
step4 = [[1.13, 1.3, 0.4], [1.43, 1.3, 0.4], [1.43, 0.3, 0.4 ], [1.13, 0.3, 0.4 ], ]
step5 = [[1.43, 1.3, 0.5], [1.73, 1.3, 0.5], [1.73, 0.3, 0.5 ], [1.43, 0.3, 0.5 ], ]
step6 = [[1.73, 1.3, 0.6], [2.03, 1.3, 0.6], [2.03, 0.3, 0.6 ], [1.73, 0.3, 0.6 ], ]
step7 = [[1.73, 1.7, 0.6], [2.73, 1.7, 0.6], [2.73,-0.3, 0.6 ], [1.73,-0.3, 0.6 ], ]

all_surfaces = [floor, step1, step2, step3, step4,step5,step6]

afloor = array(floor).T 
astep1 = array(step1).T
astep2 = array(step2).T
astep3 = array(step3).T
astep4 = array(step4).T
astep5 = array(step5).T
astep6 = array(step6).T
#~ astep7 = array(step7).T

def gen_stair_pb():
    #~ for i in range(nphases)
    #~ kinematicConstraints = genKinematicConstraints(min_height = 0.6)
    kinematicConstraints = genKinematicConstraints(min_height = None)
    relativeConstraints = genFootRelativeConstraints()
    
    #~ plot_polytope_H_rep(kinematicConstraints[0][0], kinematicConstraints[0][1].reshape((-1,1)))
    
    #~ surfaces = [[afloor], [afloor], [astep1,astep2], [astep2,astep1,astep3,astep4,astep5], [astep3,astep1,astep2,astep4,astep5], [astep4,astep1,astep2,astep3,astep5],[astep5,astep1,astep2,astep3,astep4],[astep6,astep5,astep1,astep2,astep3,astep4],[astep6,astep5,astep1,astep2,astep3,astep4],[astep6,astep5,astep1,astep2,astep3,astep4],[astep6,astep5,astep1,astep2,astep3,astep4],[astep6],[astep6]]
    surfaces = [[afloor], [afloor], [astep1,astep2], [astep2,astep1,astep3,astep4,astep5], [astep3,astep1,astep2,astep4,astep5], [astep4,astep1,astep2,astep3,astep5],[astep5,astep1,astep2,astep3,astep4],[astep6],[astep6]]
    #~ for Ss in surfaces:
        #~ for el in Ss:
            #~ el[1,:] -=0.5
    #~ surfaces = [[el[0]] for el in surfaces]
    nphases = len(surfaces)
    p0 = None
    
    res = { "p0" : None, "nphases": nphases}
    
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

def solve():    
    
    pb = gen_stair_pb()
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    
    #~ C, c = least_square_cost_function (pb, A.shape[1], comTarget = array([30.,0.,0.5]), reg_term = 0.0000001)
    C = identity(A.shape[1])
    c = zeros(A.shape[1])
    t2 = clock()
    res = qp.solve_least_square(C,c,A,b,E,e)
    t3 = clock()
    
    
    print("time to set up problem" , timMs(t1,t2))
    print("time to solve problem"  , timMs(t2,t3))
    print("total time"             , timMs(t1,t3))
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    return pb, coms, footpos, allfeetpos, res

if __name__ == '__main__':
    
    from sl1m.planner_l1 import *
    #~ draw_rectangle(l,ax)
    #~ plt.show()
    
    pb = gen_stair_pb()
    
    t1 = clock()
    A, b, E, e = convertProblemToLp(pb)
    A,b = normalize([A,b])
    
    #~ C, c = least_square_cost_function (pb, A.shape[1], comTarget = array([30.,0.,0.5]), reg_term = 0.0000001)
    C = identity(A.shape[1]) * 0.00001
    #~ c = slackSelectionMatrix(pb)
    c = zeros(A.shape[1])
    t2 = clock()
    #~ res = qp.quadprog_solve_qp(C, c,A,b,E,e)
    #~ t2 = clock()
    #~ res = qp.solve_least_square(C,c,A,b,E,e)
    res = qp.solve_lp(c,A,b,E,e)["x"]
    t3 = clock()
    
    
    print("time to set up problem" , timMs(t1,t2))
    print("time to solve problem"  , timMs(t2,t3))
    print("total time"             , timMs(t1,t3))
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    
    
    ax = draw_scene()
    #~ plotQPRes(pb, res, ax=ax, plot_constraints=True)
    plotQPRes(pb, res, ax=ax, plot_constraints=False)
    
    
