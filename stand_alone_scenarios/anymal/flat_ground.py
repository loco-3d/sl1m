import numpy as np

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm


from sl1m.constants_and_tools import *
from sl1m.planner_l1_generic_equalities_as_ineq import *



from sl1m.stand_alone_scenarios.constraints_anymal import *

floor = [ [-10, -10., 0.], [10, -10., 0.], [10, 10., 0.], [-10, 10., 0.]  ]


afloor = array(floor).T 
all_surfaces = [floor]

surfaces = [[afloor] for _ in range (10)]

def gen_flat_pb():    
    kinematicConstraints = genCOMConstraints()
    relativeConstraints = genRelativeConstraints()
    nphases = len(surfaces)
    p0 = None
    res = { "p0" : p0, "c0" : None, "nphases": nphases}    
    phaseData = [ {"K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i]))], "S" : surfaces[i], "allRelativeK" : [relativeConstraints] } for i in range(nphases)]
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
    
    
    
############# main ###################    

if __name__ == '__main__':
    
    
    from sl1m.planner_l1_generic_equalities_as_ineq import solveMIPGurobi, plotQPRes, initGlobals 
    # ~ from sl1m.fix_sparsity import solveL1, solveMIP
    
    
    initGlobals(nEffectors = 4)  
    pb = gen_flat_pb()  
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False)
    ax = draw_scene(None)
    plotQPRes(pb, res, ax=ax, plot_constraints = False)
    
    
    
    
