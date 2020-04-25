import numpy as np

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm


from sl1m.constants_and_tools import *
from sl1m.planner_l1_generic_equalities_as_ineq import *
from hpp.corbaserver.rbprm import  rbprmstate
from hpp.corbaserver.rbprm import  state_alg
from sl1m.stand_alone_scenarios.constraints_anymal import *

all_surfaces = None
all_surfaces_array = None
default_surfaces = None


surfaces_and_normals = [
# first platform
([[-0.58, -1.2, 0.39],
   [-0.6, 0., 0.07],
   [0.61, -1.21, 0.39]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
([[-0.6, 0., 0.07],
[0.61,-0.004,0.09],
   [0.61, -1.21, 0.39]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
# second platform
([[-.59,0.,0.06],
[-0.57,1.21,0.4],
   [0.61, -0.004, 0.09]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
([[-0.57,1.21,0.4],
[0.6, 1.18, 0.43],
   [0.61, -0.004, 0.09]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
# third platform
([[0.6,-1.21,0.07],
[0.62,-0.009,0.435],
   [1.84, -1.22, 0.1]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
([[0.62,-0.009,0.435],
[1.82, -0.009, 0.48],
   [1.84, -1.22, 0.1]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
# fourth platform
([[0.62,-0.009,0.435],
[0.6,1.22,0.1],
   [1.82, -0.009, 0.48]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
([[0.6,1.22,0.1],
[1.84, 1.22, 0.13],
   [1.82, -0.009, 0.48]],
  [-0.007255147224555951, -0.277401422155636, 0.9607267113101313]),
 ]



all_surfaces = [el[0] for el in surfaces_and_normals]

all_surfaces_array = [array(el).T for el in all_surfaces]

default_surfaces = [all_surfaces_array for _ in range (10)]


surfaces = None

def genSurf():
    global surfaces
    surfaces = [all_surfaces_array for _ in range (6)]

genSurf()


def draw_rectangle(l, ax):
    #~ plotPoints(ax,l)
    lr = l + [l[0]]
    cx = [c[0] for c in lr]
    cy = [c[1] for c in lr]
    cz = [c[2] for c in lr]
    ax.plot(cx, cy, cz)

def draw_scene(surfaces = all_surfaces, ax = None, color = "p"):
        if ax is None:        
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
        [draw_rectangle(l,ax) for l in surfaces]
        return ax

def overrideSurfaces(surfs):
    global all_surfaces
    global all_surfaces_array
    global default_surfaces
    global surfaces
    all_surfaces = [el[0] for el in surfs]
    all_surfaces_array = [array(el).T for el in all_surfaces]
    default_surfaces = [all_surfaces_array for _ in range (5)]
    genSurf()
    print (' HELLLOOO ')
    draw_scene(all_surfaces)
    draw_scene()
    plt.show()
    


def gen_pb():    
    kinematicConstraints = genCOMConstraints()
    relativeConstraints = genRelativeConstraints()
    nphases = len(surfaces)
    p0 = None
    res = { "p0" : p0, "c0" : None, "nphases": nphases}    
    phaseData = [ {"K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i]))], "S" : surfaces[i], "allRelativeK" : [relativeConstraints] } for i in range(nphases)]
    res ["phaseData"] = phaseData
    return res 
    
       
    
# ~ def solve(initCom = default_init_com, initPos = default_init_pos, endCom = None):
def solve(initCom = None, initPos = None, endCom = None, endPos = None, surfaces = surfaces, ax = None):
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!initPos ", initPos)
    # ~ print ("default_init_pos ", default_init_pos)
    # ~ if endCom is None:
        # ~ endCom  = default_init_com + array([0.2, 0.0, 0.0])
    initGlobals(nEffectors = 4)  
    
    
    pb = gen_pb()  
    # ~ endPos = None
    print ("endCom", endCom)
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom, costs = [(1, posturalCost),(1.5, targetCom)])
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is not None:
        plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True)
    
    # ~ off = 2
    # ~ Coms = coms[:off];  Footpos = footpos[:off]; Allfeetpos = allfeetpos[:off]
    Coms = coms[:];  Footpos = footpos[:]; Allfeetpos = allfeetpos[:]
    
    print ("coms ", Coms)
    print ("Allfeetpos ", Allfeetpos)
    
    for i in range(2):
    # ~ for i in range(1):
        # ~ print("round !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", i )
        pb = gen_pb()  
        initPos = None
        # ~ endPos = None
        initCom = Coms[-1].copy()
        endCom  = endCom
        initPos = Allfeetpos[-1].copy()
        pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom,)
        coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
        Coms += coms[1:];  Footpos += footpos[1:]; Allfeetpos += allfeetpos[1:]        
        # ~ if ax is not None:
    
    for i in range(3):
    # ~ for i in range(1):
        # ~ print("round !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", i )
        pb = gen_pb()  
        initPos = None
        # ~ endPos = None
        initCom = Coms[-1].copy()
        endCom  = Coms[-1].copy() + array([20, 0.0, 0.0])
        initPos = Allfeetpos[-1].copy()
        pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom,
        costs = [(1., targetCom)])
        coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
        Coms += coms[1:];  Footpos += footpos[1:]; Allfeetpos += allfeetpos[1:]        
        # ~ if ax is not None:
            # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True)
    
    return pb, Coms, Footpos, Allfeetpos, res
    
    
############# main ###################    

if __name__ == '__main__':    
    
    
    def plotPoints(ax, wps, color = "b", D3 = True, linewidth=2):
        x = array(wps)[:,0]
        y = array(wps)[:,1]
        if(D3):                
                z = array(wps)[:,2]
                ax.scatter(x, y, z, color=color, marker='o', linewidth = 5) 
        else:
                ax.scatter(x,y,color=color, linewidth = linewidth)  
    
    def draw_rectangle(l, ax):
        #~ plotPoints(ax,l)
        lr = l + [l[0]]
        cx = [c[0] for c in lr]
        cy = [c[1] for c in lr]
        cz = [c[2] for c in lr]
        ax.plot(cx, cy, cz)
    
    def draw_scene(surfaces = all_surfaces, ax = None, color = "p"):
        if ax is None:        
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
        [draw_rectangle(l,ax) for l in surfaces]
        return ax
        
    # ~ initPos = [array([-0.8, 0.2, 0.13]),
 # ~ array([-0.8, -0.2, 0.13]),
 # ~ array([-1.6, 0.2, 0.13]),
 # ~ array([-1.6, -0.2, 0.13])]
 
    initPos = [array([0.60700772, 0.35028565, 0.19173743]),
 array([0.38222138, 0.23795264, 0.15202783]),
 array([-0.36776191,  0.34361369,  0.23769671]),
 array([-0.25487645, -0.23401013,  0.25597898])]

    initPos = [array([0.38437629, 0.18974121, 0.13850503]),
 array([ 0.38437629, -0.18974121,  0.1326605 ]),
 array([-0.38437629,  0.18974121,  0.11856727]),
 array([-0.38437629, -0.18974121,  0.05008679])]

    
    # ~ initPos = [array([0.6, 0.3, 0.2]),
 # ~ array([0.6, -0.5, 0.2]),
 # ~ array([0.2, 0.3, 0.2]),
 # ~ array([0.5, -0.5, 0.2])]


    
    endCom = [0, 1, 0.4]
    
    ax = draw_scene()
    plotPoints(ax, initPos)
    # ~ plt.show()
    pb, coms, footpos, allfeetpos, res = solve(initPos = initPos, endCom=endCom, endPos = None, ax = ax)
    
    
    # ~ for i in range(15):
        # ~ print ("WTTTTTTTTTTFFFFFFFFFFFFFFf!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", i)
        # ~ coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
        # ~ pb = gen_pb()  
        # ~ initPos = None
        # ~ endPos = None
        # ~ initCom = coms[-1]
        # ~ endCom  = coms[-1] + array([6, 0.0, 0.0])
        # ~ initPos = allfeetpos[-1]
        # ~ pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
        # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True)
        
    
    
    
    # ~ plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True, ax = ax)
    plt.show(block = False)
    
    
    
    
