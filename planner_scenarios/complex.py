import numpy as np


from sl1m.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *
#~ from qp import solve_lp

from sl1m.planner import *


from plot_plytopes import *

floor2  = [[-3., 0.4  , 0.  ], [-2.7 ,  0.4, 0. ], [-2.7 , 0.1, 0.  ], [-03., 0.1, 0.  ], ]
floor = [[-0.30, 0.54  , 0.  ], [0.01 ,  0.54, 0. ], [0.01 , -0.46, 0.  ], [-0.30, -0.46, 0.  ], ]
step1  = [[ 0.01, 0.54  , 0.1 ], [0.31 ,  0.54, 0.1], [0.31 , -0.46, 0.1 ], [ 0.01, -0.46, 0.1 ], ]
step2  = [[ 0.31, 0.54  , 0.2 ], [0.61 ,  0.54, 0.2], [0.61 , -0.46, 0.2 ], [ 0.31, -0.46, 0.2 ], ]
step3  = [[ 0.61, 0.54  , 0.3 ], [0.91 ,  0.54, 0.3], [0.91 , -0.46, 0.3 ], [ 0.61, -0.46, 0.3 ], ]
step4  = [[ 0.91, 0.54  , 0.4 ], [1.21 ,  0.54, 0.4], [1.21 , -0.46, 0.4 ], [ 0.91, -0.46, 0.4 ], ]
step5  = [[ 1.21, 0.54  , 0.5 ], [1.51 ,  0.54, 0.5], [1.51 , -0.46, 0.5 ], [ 1.21, -0.46, 0.5 ], ]
step6  = [[ 1.51, 0.54  , 0.6 ], [1.81 ,  0.54, 0.6], [1.81 , -0.46, 0.6 ], [ 1.51, -0.46, 0.6 ], ]
#~ step7  = [[ 1.51, 0.94  , 0.6 ], [2.51 ,  0.94, 0.6], [2.51 , -1.06, 0.6 ], [ 1.51, -1.06, 0.6 ], ]
step7  = [[ 1.51,-0.46  , 0.6 ], [1.81 , -0.46, 0.6], [1.81 , -0.76, 0.6 ], [ 1.51, -0.76, 0.6 ], ]
bridge = [[ 1.51, -0.46 , 0.6 ], [1.51 , -0.76, 0.6], [-1.49, -0.76, 0.6 ], [-1.49, -0.46, 0.6 ], ]
#~ platfo = [[-1.49, -0.06 , 0.6 ], [-1.49, -1.06, 0.6], [-2.49, -1.06, 0.6 ], [-2.49, -0.06, 0.6 ], ]
platfo = [[-1.49, -0.35, 0.6 ], [-1.49, -1.06, 0.6], [-2.49, -1.06, 0.6 ], [-2.49, -0.35, 0.6 ], ]
#~ step8 =  [[-1.49, -0.06 , 0.45], [-1.49,  0.24, 0.45],[-2.49,  0.24, 0.45], [-2.49, -0.06, 0.45], ]
#~ step9 =  [[-1.49,  0.24 , 0.30], [-1.49,  0.54, 0.30],[-2.49,  0.54, 0.30], [-2.49,  0.24, 0.30], ]
#~ step10=  [[-1.49,  0.54 , 0.15], [-1.49,  0.84, 0.15],[-2.49,  0.84, 0.15], [-2.49,  0.54, 0.15], ]
slope = [[-1.49, -0.06 , 0.6 ], [-1.49, 1.5, 0.], [-2.49, 1.5, 0. ], [-2.49, -0.06, 0.6 ], ]
rub2  = [[ -2.11, 0.30 , 0.05 ], [-2.45 , 0.30, 0.05 ],  [ -2.45, 0.53, 0.05 ], [-2.11, 0.53, 0.05  ], ]
rub3  = [[ -1.91, -0.15 , 0.1 ], [-2.25 , -0.15, 0.1 ],  [ -2.25, 0.18, 0.1 ], [-1.91, 0.18, 0.1  ], ]
rub4  = [[ -1.69, 0.19 , 0.15 ], [-2.03 , 0.19, 0.15 ],  [ -2.03, 0.53, 0.15 ], [-1.69, 0.53, 0.15  ], ]
rub5  = [[ -1.49, -0.15 , 0.2 ], [-1.83 , -0.15, 0.2 ],  [ -1.83, 0.18, 0.2 ], [-1.49, 0.18, 0.2  ], ]
rub6  = [[ -1.29, 0.19 , 0.2 ], [-1.63 , 0.19, 0.2 ],  [ -1.63, 0.53, 0.2 ], [-1.29, 0.53, 0.2  ], ]
rub7  = [[ -1.09, -0.15 , 0.15 ], [-1.43 , -0.15, 0.15],  [ -1.43, 0.18, 0.15], [-1.09, 0.18, 0.15 ], ]
rub75  = [[ -0.89, 0.19 , 0.1 ], [-1.23 , 0.19, 0.1],  [ -1.23, 0.53, 0.1], [-0.89, 0.53, 0.1 ], ]
rub8  = [[ -0.89, -0.15 , 0.025 ], [-1.02 , -0.15, 0.025],  [ -1.02, 0.18, 0.025], [-0.89, 0.18, 0.025 ], ]
rub9  = [[ -0.35, -0.15 , 0.025 ], [-0.86 , -0.15, 0.025], [-0.86, 0.52, 0.025 ], [ -0.35, 0.52, 0.025], ]
rub8  = [[ -0.89, -0.15 , 0.05 ], [-1.02 , -0.15, 0.05],  [ -1.02, 0.18, 0.05], [-0.89, 0.18, 0.05 ], ]
rub9  = [[ -0.35, -0.15 , 0.05 ], [-0.86 , -0.15, 0.05], [-0.86, 0.52, 0.05 ], [ -0.35, 0.52, 0.05], ]

all_surfaces = [floor2, floor, step1, step2, step3, step4,step5,step6, step7, bridge, platfo, rub8, rub9,rub7, rub75, rub6, rub5, rub4, rub3, rub2]

arub9  = array(rub9).T 
arub8  = array(rub8).T 
arub75  = array(rub75).T 
arub7  = array(rub7).T 
arub6  = array(rub6).T 
arub5  = array(rub5).T 
arub4  = array(rub4).T 
arub3  = array(rub3).T 
arub2  = array(rub2).T 
#~ arub1  = array(rub1).T 
afloor  = array(floor).T 
astep1  = array(step1).T
astep2  = array(step2).T
astep3  = array(step3).T
astep4  = array(step4).T
astep5  = array(step5).T
astep6  = array(step6).T
astep7  = array(step7).T
#~ astep8  = array(step8).T
#~ astep9  = array(step9).T
#~ astep10  = array(step10).T
abridge = array(bridge).T
aplatfo = array(platfo).T
aslope  = array(slope).T


allrub = [arub2,arub3,arub5,arub4,arub6,arub7,arub75,arub9]
allsteps = [astep2,astep1,astep3,astep4,astep5,astep6,astep7]
allstepsandfloor = allsteps + [arub9,afloor]
    

#~ surfaces0 = [[arub2,arub3],[arub3,arub2],[arub4,arub3,arub5],allrub,allrub,allrub,[arub75] ,[arub9,afloor],[arub9,afloor],[afloor,arub9],[astep1], [astep2,astep3,astep4,astep5],[astep3,astep2,astep4,astep5],[astep4,astep1,astep2,astep3,astep5], [astep5,astep4,astep1,astep2,astep3],[astep6,astep5,astep4],[astep6],[astep6],[astep6,],[astep7],[astep7],[abridge,aplatfo],[abridge,aplatfo],[abridge,aplatfo],[abridge,astep7,aplatfo]]

surfaces0 = [[arub2,arub3],[arub3,arub2],[arub4,arub3,arub5],allrub,allrub,allrub,[arub75] ,allstepsandfloor,allstepsandfloor,allstepsandfloor,allsteps, allsteps,allsteps,allsteps, allsteps,[astep6],[astep6,],[astep7],[astep7],[abridge,aplatfo],[abridge,aplatfo],[abridge,aplatfo],[abridge,aplatfo]]

surfaces1 = [  [abridge], [abridge],[abridge,aplatfo],[abridge,aplatfo],[abridge,aplatfo],[abridge,aplatfo],[aplatfo],[aplatfo],[aplatfo],[aplatfo]]
surfaces = surfaces0 + surfaces1

def gen_stair_pb():
    #~ for i in range(nphases)
    #~ kinematicConstraints = genKinematicConstraints(min_height = 0.6)
    kinematicConstraints = genKinematicConstraints(min_height = None)
    relativeConstraints = genFootRelativeConstraints()
        
    #~ arrayAllSurf = [array(el).T for el in all_surfaces]
    #~ surfaces = [[arub2,arub3],[arub3,arub2]] + [ arrayAllSurf for _ in range(30)] + [[aplatfo],[aplatfo]]
    
    nphases = len(surfaces)
    
    p0 = [array([-2.7805096486250154, 0.33499999999999996, 0.]), array([-2.7805096486250154, 0.145,0.])];
    
    #~ res = { "p0" : p0, "c0" : None, "nphases": nphases}
    res = { "p0" : None, "c0" : None, "nphases": nphases}
    
    #TODO in non planar cases, K must be rotated
    phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i]))], "relativeK" : [relativeConstraints[(i)%2] for _ in range(len(surfaces[i]))], "S" : surfaces[i] } for i in range(nphases)]
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
    
    from sl1m.fix_sparsity import solveL1
    
    pb = gen_stair_pb()
    
    pb, coms, footpos, allfeetpos, res = solveL1(pb, surfaces, draw_scene)
    #~ ax = draw_scene()
    #~ plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    
