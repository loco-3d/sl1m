import numpy as np

from sl1m.constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm


from sl1m.planner import *
from sl1m.planner_scenarios.talos.constraints import *

from sl1m.tools.plot_plytopes import *

Z_AXIS = np.array([0,0,1]).T

begin = array([[1.75, 1.65, 1.65, 1.75],
        [0.3, 0.3, 0.1, 0.1],
        [0.6, 0.6, 0.6, 0.6]])

platform = array([[2.5, 1.5, 1.5, 2.5],
        [0.9, 0.9, -1.1, -1.1],
        [0.6, 0.6, 0.6, 0.6]])

bridge = array([[-1.5, -1.5, 1.5, 1.5],
        [-0.5, -0.8, -0.8, -0.5],
        [0.6, 0.6, 0.6, 0.6]])

end = array([[-1.5, -2.0, -2.0, -1.5],
        [-0.4, -0.4, -1.1, -1.1],
        [0.6, 0.6, 0.6, 0.6]])

### END HARDCODED SURFACES ####

surfaces = []
surfaces += [[begin]]+[[begin]]
for i in range(20):
  surfaces += [[platform,bridge]]
surfaces += [[end]]+[[end]]


def gen_pb(surfaces):
    kinematicConstraints = genKinematicConstraints(left_foot_constraints,right_foot_constraints,min_height = 0.3)
    relativeConstraints = genFootRelativeConstraints(right_foot_in_lf_frame_constraints,left_foot_in_rf_frame_constraints)
    
    nphases = len(surfaces)
    p0 = [array([ 1.7  ,  0.285,  0.6  ]), array([ 1.7  ,  0.115,  0.6  ])]
    print "p0 used : ",p0
    
    #res = { "p0" : p0, "c0" : None, "nphases": nphases}
    res = { "p0" : None, "c0" : None, "nphases": nphases}
    
    print "surfaces = ",surfaces
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
    R,surfaces = getSurfacesFromGuide(tp.rbprmBuilder,tp.ps,tp.afftool,tp.pathId,tp.v,1.,False)

    pb = gen_pb(tp.q_init,R,surfaces)

    return solveL1(pb, surfaces, draw_scene)

if __name__ == '__main__':
    from sl1m.fix_sparsity import solveL1    

    pb = gen_pb(surfaces)
    """
    import pickle
    f = open("pb_platform_bridge_noGuide","w")
    pickle.dump(pb,f)
    f.close()
    """
    pb, coms, footpos, allfeetpos, res = solveL1(pb, surfaces, draw_scene)

    
"""    
import matplotlib.pyplot as plt    
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
 
def plotSurface (points, ax, plt, c = 'rand'):
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    xs = np.append(xs, xs[0]).tolist()
    ys = np.append(ys, ys[0]).tolist()
    zs = np.append(zs, zs[0]).tolist()
    if c == 'rand': ax.plot(xs,ys,zs)
    else: ax.plot(xs,ys,zs,c)
    plt.draw()
       
def draw_scene(ax = None, color = "p"):
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    #[draw_rectangle(l,ax) for l in all_surfaces]
    for surface in my_surfaces: plotSurface(surface[0], ax, plt)
    return ax

global my_surfaces
my_surfaces = tp.surfaces    
ax = draw_scene()
plt.show()
"""

