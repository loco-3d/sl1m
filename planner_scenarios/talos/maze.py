import numpy as np
print "Plan guide trajectory ..."
import scenarios.sandbox.talos_maze_path as tp
v = tp.v
ps = tp.ps
root_init = tp.q_init
print "Guide planned."
from tools.surfaces_from_path import getSurfacesFromGuideContinuous

from sl1m.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm
import random
#~ from scipy.spatial import ConvexHull
#~ from hpp_bezier_com_traj import *
#~ from qp import solve_lp

from sl1m.planner import *
from sl1m.planner_scenarios.talos.constraints import *

Z_AXIS = np.array([0,0,1]).T


def gen_pb(root_init,R, surfaces):

    nphases = len(surfaces)
    lf_0 = array(root_init[0:3]) + array([0, 0.085,-0.98]) # values for talos ! 
    rf_0 = array(root_init[0:3]) + array([0,-0.085,-0.98]) # values for talos ! 
    p0 = [lf_0,rf_0];
    print "p0 used : ",p0
    
    res = { "p0" : p0, "c0" : None, "nphases": nphases}
    #res = { "p0" : None, "c0" : None, "nphases": nphases}
    
    #print "surfaces = ",surfaces
    #TODO in non planar cases, K must be rotated
    phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : [genKinematicConstraints(left_foot_constraints,right_foot_constraints,index = i, rotation = R, min_height = 0.3) for _ in range(len(surfaces[i]))], "relativeK" : [genFootRelativeConstraints(right_foot_in_lf_frame_constraints,left_foot_in_rf_frame_constraints,index = i, rotation = R)[(i) % 2] for _ in range(len(surfaces[i]))], "rootOrientation" : R[i], "S" : surfaces[i] } for i in range(nphases)]
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

def plotSurface (points, ax, plt,color_id = -1):
    xs = np.append(points[0,:] ,points[0,0] ).tolist()
    ys = np.append(points[1,:] ,points[1,0] ).tolist()
    zs = (np.append(points[2,:] ,points[2,0] ) - np.ones(len(xs))*0.005*color_id).tolist()
    colors = ['r','g','b','m','y','c']
    if color_id == -1: ax.plot(xs,ys,zs)
    else: ax.plot(xs,ys,zs,colors[color_id])
    plt.draw()
        
def draw_scene(surfaces,ax = None):
    colors = ['r','g','b','m','y','c']
    color_id = 0
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    #[draw_rectangle(l,ax) for l in all_surfaces]
    for surfaces_phase in surfaces: 
      for surface in surfaces_phase:
        plotSurface(surface, ax, plt,color_id)
      color_id += 1
      if color_id >= len(colors):
        color_id = 0
    return ax
    
    
    
############# main ###################    
def solve():
    from sl1m.fix_sparsity import solveL1
    success = False
    maxIt = 50
    it = 0
    defaultStep = 0.5
    step = defaultStep
    variation = 0.4
    while not success and it < maxIt:
      if it > 0 :
        step = defaultStep + random.uniform(-variation,variation)
      R,surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder,tp.ps,tp.afftool,tp.pathId,tp.v,step,True,False)
      pb = gen_pb(tp.q_init,R,surfaces)
      try:
        pb, coms, footpos, allfeetpos, res = solveL1(pb, surfaces, None)
        success = True
      except :  
        print "## Planner failed at iter : "+str(it)+" with step length = "+str(step)
      it += 1
    if not success :
      raise RuntimeError("planner always fail.") 
    return pb, coms, footpos, allfeetpos, res

if __name__ == '__main__':
    pb, coms, footpos, allfeetpos, res = solve()

