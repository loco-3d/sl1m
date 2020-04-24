import numpy as np

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm


from sl1m.constants_and_tools import *
from sl1m.planner_l1_generic_equalities_as_ineq import *
from hpp.corbaserver.rbprm import  rbprmstate
from hpp.corbaserver.rbprm import  state_alg
from sl1m.stand_alone_scenarios.constraints_anymal import *

from sl1m.stand_alone_scenarios.anymal.plinth_surfaces import all_surfaces, all_surfaces_array, default_surfaces, draw_scene, default_init_pos, default_init_com


surfaces = None

def genSurf():
    global surfaces
    surfaces = [all_surfaces_array for _ in range (5)]

genSurf()

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
    
    
    ax = draw_scene(all_surfaces)
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
def solve(initCom = None, initPos = default_init_pos, endCom = None, surfaces = surfaces, ax = None):
    print ("initPos ", initPos)
    print ("default_init_pos ", default_init_pos)
    if endCom is None:
        endCom  = default_init_com + array([0.2, 0.0, 0.0])
    initGlobals(nEffectors = 4)  
    
    
    pb = gen_pb()  
    endPos = None
    # ~ print ("initPos", initPos)
    pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
    coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
    if ax is not None:
        plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True)
    
    # ~ off = 2
    # ~ Coms = coms[:off];  Footpos = footpos[:off]; Allfeetpos = allfeetpos[:off]
    Coms = coms[:];  Footpos = footpos[:]; Allfeetpos = allfeetpos[:]
    
    print ("coms ", Coms)
    print ("Allfeetpos ", Allfeetpos)
    
    for i in range(10):
        print("round !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", i )
        pb = gen_pb()  
        initPos = None
        endPos = None
        initCom = Coms[-1].copy()
        endCom  = Coms[-1].copy() + array([0.2, 0.0, 0.0])
        initPos = Allfeetpos[-1].copy()
        pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
        coms, footpos, allfeetpos = retrieve_points_from_res(pb, res)
        Coms += coms[:];  Footpos += footpos[:]; Allfeetpos += allfeetpos[:]        
        if ax is not None:
            plotQPRes(pb, res, ax=ax, plot_constraints = False, show = False, plotSupport = True)
    
    return pb, Coms, Footpos, Allfeetpos, res
    
    
############# main ###################    

if __name__ == '__main__':    
    
    ax = draw_scene()
    pb, coms, footpos, allfeetpos, res = solve(ax = ax)
    
    
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
    
    
    
    
