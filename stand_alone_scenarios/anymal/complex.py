import numpy as np

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm


from sl1m.constants_and_tools import *
from sl1m.planner_l1_generic_equalities_as_ineq import *



from sl1m.stand_alone_scenarios.constraints_anymal import *


start  = [[-3., 0.4  , 0.  ], [-2.7 ,  0.4, 0. ], [-2.7 , 0.1, 0.  ], [-03., 0.1, 0.  ], ]
floor = [[-0.30, 0.54  , 0.  ], [0.01 ,  0.54, 0. ], [0.01 , -0.46, 0.  ], [-0.30, -0.46, 0.  ], ]
step1  = [[ 0.01, 0.54  , 0.1 ], [0.31 ,  0.54, 0.1], [0.31 , -0.46, 0.1 ], [ 0.01, -0.46, 0.1 ], ]
step2  = [[ 0.31, 0.54  , 0.2 ], [0.61 ,  0.54, 0.2], [0.61 , -0.46, 0.2 ], [ 0.31, -0.46, 0.2 ], ]
step3  = [[ 0.61, 0.54  , 0.3 ], [0.91 ,  0.54, 0.3], [0.91 , -0.46, 0.3 ], [ 0.61, -0.46, 0.3 ], ]
step4  = [[ 0.91, 0.54  , 0.4 ], [1.21 ,  0.54, 0.4], [1.21 , -0.46, 0.4 ], [ 0.91, -0.46, 0.4 ], ]
step5  = [[ 1.21, 0.54  , 0.5 ], [1.51 ,  0.54, 0.5], [1.51 , -0.46, 0.5 ], [ 1.21, -0.46, 0.5 ], ]
step6  = [[ 1.51, 0.54  , 0.6 ], [1.81 ,  0.54, 0.6], [1.81 , -0.46, 0.6 ], [ 1.51, -0.46, 0.6 ], ]
step7  = [[ 1.51,-0.46  , 0.6 ], [1.81 , -0.46, 0.6], [1.81 , -0.76, 0.6 ], [ 1.51, -0.76, 0.6 ], ]
# ~ bridge = [[ 1.51, -0.46 , 0.6 ], [1.51 , -0.76, 0.6], [-1.49, -0.76, 0.6 ], [-1.49, -0.46, 0.6 ], ]
bridge = [[ 1.51, -0.06 , 0.6 ], [1.51 , -0.96, 0.6], [-1.49, -0.96, 0.6 ], [-1.49, -0.06, 0.6 ], ]
platfo = [[-1.49, -0.35, 0.6 ], [-1.49, -1.06, 0.6], [-2.49, -1.06, 0.6 ], [-2.49, -0.35, 0.6 ], ]
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

all_surfaces = [start, floor, step1, step2, step3, step4,step5,step6, step7, bridge, platfo, rub8, rub9,rub7, rub75, rub6, rub5, rub4, rub3, rub2]

astart = array(start).T
arub9  = array(rub9).T 
arub8  = array(rub8).T 
arub75  = array(rub75).T 
arub7  = array(rub7).T 
arub6  = array(rub6).T 
arub5  = array(rub5).T 
arub4  = array(rub4).T 
arub3  = array(rub3).T 
arub2  = array(rub2).T 
astart  = array(start).T 
afloor  = array(floor).T 
astep1  = array(step1).T
astep2  = array(step2).T
astep3  = array(step3).T
astep4  = array(step4).T
astep5  = array(step5).T
astep6  = array(step6).T
astep7  = array(step7).T
abridge = array(bridge).T
aplatfo = array(platfo).T
aslope  = array(slope).T


allrub = [arub2,arub3,arub5,arub4,arub6,arub7,arub75,arub9]
allsteps = [astep2,astep1,astep3,astep4,astep5,astep6,astep7]
allstepsandfloor = allsteps + [arub9,afloor]
allrubfloorsteps = allrub + allsteps + [afloor]
end = [astep5, astep6, astep7, abridge, aplatfo ]

# ~ surfacesred = [[abridge],end,end,end, [aplatfo]]

surfaces = [[astep7],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],[abridge],end,end,end, [aplatfo]]
surfaces = [[astart],allrubfloorsteps,allrubfloorsteps,allrubfloorsteps,allrubfloorsteps, [arub7]]

def gen_flat_pb(surfs = surfaces):    
    kinematicConstraints = genCOMConstraints()
    relativeConstraints = genRelativeConstraints()
    nphases = len(surfs)
    p0 = None
    res = { "p0" : p0, "c0" : None, "nphases": nphases}    
    phaseData = [ {"K" : [copyKin(kinematicConstraints) for _ in range(len(surfs[i]))], "S" : surfs[i], "allRelativeK" : [relativeConstraints] } for i in range(nphases)]
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
    
def draw_scene(surfaces, ax = None, color = "p"):
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    [draw_rectangle(l,ax) for l in all_surfaces]
    return ax
    
    
def solve():
    initGlobals(nEffectors = 4)  
    pb = gen_flat_pb()  
    initPos = None
    endPos = None
    initCom = None
    endCom  = None
    # ~ initCom = [0.2, 0., 0.5]
    # ~ endCom  = [2, 0.7, 0.5]
    # ~ endCom  = [0.2, 1.2, 0.5]
    # ~ initPos = [array([0.5, 0.0, 0.0]),  array([-0.1, 0.0, -0.0]),  array([0.3, 0.5, 0.0]),  array([-0.0, 0.4, 0.0])]
    # ~ endPos = [array()]
    # ~ initPos = [[0.2, 0.6, 0.0] ];     initPos = [array(el) for el in initPos]
    # ~ endPos = [[0.7, 0.6, 0.0] ];  endPos = [array(el) for el in endPos]
    # ~ initCom = [0.2, 2., 0.35]
    # ~ endCom = [0.2, 2., 0.35]
    print ("initPos", initPos)
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    return pb, coms, footpos, allfeetpos, res
############# main ###################    

if __name__ == '__main__':
    
    
    from sl1m.planner_l1_generic_equalities_as_ineq import solveMIPGurobi, plotQPRes, initGlobals 
    # ~ from sl1m.fix_sparsity import solveL1, solveMIP
    
    
    initGlobals(nEffectors = 4)  
    # ~ pb = gen_flat_pb([[afloor],[afloor],[afloor],[afloor]])  
    pb = gen_flat_pb()  
    # ~ initPos = [[0.2, 0.6, 0.0],[-0.7, 0.1, 0.0], [-0.3, 1.2, 0.0], [-0.8, 1.0, 0.0] ]
    initPos = None
    endPos = None
    initCom = None
    endCom  = None
    # ~ initCom = [0.2, 0., 0.5]
    # ~ endCom  = [2, 0.7, 0.5]
    # ~ endCom  = [0.2, 1.2, 0.5]
    # ~ initPos = [array([0.5, 0.0, 0.0]),  array([-0.1, 0.0, -0.0]),  array([0.3, 0.5, 0.0]),  array([-0.0, 0.4, 0.0])]
    # ~ endPos = [array()]
    # ~ initPos = [[0.2, 0.6, 0.0] ];     initPos = [array(el) for el in initPos]
    # ~ endPos = [[0.7, 0.6, 0.0] ];  endPos = [array(el) for el in endPos]
    # ~ initCom = [0.2, 2., 0.35]
    # ~ endCom = [0.2, 2., 0.35]
    print ("initPos", initPos)
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False,l1Surface = False,initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
    
    
    off = 7
    
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    initPos = allfeetpos[-off]
    endCom = coms[-1]
    
    # ~ pb = gen_flat_pb(surfaces[-off:])  
    pb = gen_flat_pb()  
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False,l1Surface = False,initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
    
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False, show = True)
    
    # ~ pb2 = gen_flat_pb([a_all_surfaces,a_all_surfaces,a_all_surfaces,a_all_surfaces,a_all_surfaces,a_all_surfaces,a_all_surfaces, a_all_surfaces, a_all_surfaces, [astep3]])  
    # ~ pb2 = gen_flat_pb(surfaces[-off:])  
    pb2 = gen_flat_pb()  
    pb, res, time = solveMIPGurobi(pb2, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False,l1Surface = True,initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
    
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False, show = True)
    
    
    
    
