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


surfaces_and_normals = [([[0.39550296660292134, -0.27720715031876186, 0.16499999165534973],
   [0.39550296660292134, -0.44479297175155064, 0.16499999165534973],
   [0.5644970119394066, -0.44479297175155064, 0.16499999165534973],
   [0.5644970119394066, -0.27720715031876186, 0.16499999165534973]],
  [0.0, 0.0, 1.0]),
 ([[0.3955029688071581, 0.20279284117516114, 0.16499999165534973],
   [0.3955029688071581, 0.035207045099176625, 0.16499999165534973],
   [0.5644970097351698, 0.035207045099176625, 0.16499999165534973],
   [0.5644970097351698, 0.20279284117516114, 0.16499999165534973]],
  [0.0, 0.0, 1.0]),
 ([[-0.5644969567065002, -0.27720715472723556, 0.16499999165534973],
   [-0.5644969567065002, -0.44479296734307694, 0.16499999165534973],
   [-0.39550296223118286, -0.44479296734307694, 0.16499999165534973],
   [-0.39550296223118286, -0.27720715472723556, 0.16499999165534973]],
  [0.0, 0.0, 1.0]),
 ([[-0.5644969545022632, 0.20279283676668725, 0.16499999165534973],
   [-0.5644969545022632, 0.03520704950765053, 0.16499999165534973],
   [-0.3955029644354199, 0.03520704950765053, 0.16499999165534973],
   [-0.3955029644354199, 0.20279283676668725, 0.16499999165534973]],
  [0.0, 0.0, 1.0]),
 ([[0.39550296660292134, -0.03720714078201872, 0.2449999898672104],
   [0.39550296660292134, -0.20479296221480744, 0.2449999898672104],
   [0.5644970119394066, -0.20479296221480744, 0.2449999898672104],
   [0.5644970119394066, -0.03720714078201872, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[0.3955029688071581, 0.44279285071190433, 0.2449999898672104],
   [0.3955029688071581, 0.27520705463591977, 0.2449999898672104],
   [0.5644970097351698, 0.27520705463591977, 0.2449999898672104],
   [0.5644970097351698, 0.44279285071190433, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.5644969567065002, -0.03720714519049237, 0.2449999898672104],
   [-0.5644969567065002, -0.2047929578063338, 0.2449999898672104],
   [-0.39550296223118286, -0.2047929578063338, 0.2449999898672104],
   [-0.39550296223118286, -0.03720714519049237, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.5644969545022632, 0.4427928463034304, 0.2449999898672104],
   [-0.5644969545022632, 0.2752070590443937, 0.2449999898672104],
   [-0.3955029644354199, 0.2752070590443937, 0.2449999898672104],
   [-0.3955029644354199, 0.4427928463034304, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.08449704809882819, -0.037207136373545634, 0.2449999898672104],
   [-0.08449704809882819, -0.20479296662328053, 0.2449999898672104],
   [0.08449704809882819, -0.20479296662328053, 0.2449999898672104],
   [0.08449704809882819, -0.037207136373545634, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.08449704589459173, 0.4427928551203777, 0.2449999898672104],
   [-0.08449704589459173, 0.2752070502274464, 0.2449999898672104],
   [0.08449704589459173, 0.2752070502274464, 0.2449999898672104],
   [0.08449704589459173, 0.4427928551203777, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.32449705763557135, -0.2772071161079664, 0.2449999898672104],
   [-0.32449705763557135, -0.44479294635770134, 0.2449999898672104],
   [-0.15550296143791498, -0.44479294635770134, 0.2449999898672104],
   [-0.15550296143791498, -0.2772071161079664, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.3244970554313349, 0.2027928753859569, 0.2449999898672104],
   [-0.3244970554313349, 0.03520707049302565, 0.2449999898672104],
   [-0.15550296364215144, 0.03520707049302565, 0.2449999898672104],
   [-0.15550296364215144, 0.2027928753859569, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[0.15550296143791498, -0.2772071161079664, 0.2449999898672104],
   [0.15550296143791498, -0.44479294635770134, 0.2449999898672104],
   [0.32449705763557135, -0.44479294635770134, 0.2449999898672104],
   [0.32449705763557135, -0.2772071161079664, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[0.15550296364215144, 0.2027928753859569, 0.2449999898672104],
   [0.15550296364215144, 0.03520707049302565, 0.2449999898672104],
   [0.3244970554313349, 0.03520707049302565, 0.2449999898672104],
   [0.3244970554313349, 0.2027928753859569, 0.2449999898672104]],
  [0.0, 0.0, 1.0]),
 ([[-0.6401920455917902, 0.28740094947165623, 0.12999999523162842],
   [-1.5621051651015692, 0.28740094947165623, 0.12999999523162842],
   [-1.5621051651015692, -0.45012955712622776, 0.12999999523162842],
   [-0.6401920455917902, -0.45012955712622776, 0.12999999523162842]],
  [0.0, -0.0, 1.0]),
 ([[1.5659601592927235, 0.46464466094067264, 0.15000000596046448],
   [0.6366708374113782, 0.46464466094067264, 0.15000000596046448],
   [0.6366708374113782, -0.46464466094067264, 0.15000000596046448],
   [1.5659601592927235, -0.46464466094067264, 0.15000000596046448]],
  [0.0, -0.0, 1.0])]


all_surfaces = [el[0] for el in surfaces_and_normals]

all_surfaces_array = [array(el).T for el in all_surfaces]

default_surfaces = [all_surfaces_array for _ in range (10)]


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
def solve(initCom = None, initPos = None, endCom = None, surfaces = surfaces, ax = None):
    print ("initPos ", initPos)
    # ~ print ("default_init_pos ", default_init_pos)
    # ~ if endCom is None:
        # ~ endCom  = default_init_com + array([0.2, 0.0, 0.0])
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
    
    for i in range(11):
    # ~ for i in range(1):
        # ~ print("round !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", i )
        pb = gen_pb()  
        initPos = None
        endPos = None
        initCom = Coms[-1].copy()
        endCom  = Coms[-1].copy() + array([.2, 0.0, 0.0])
        initPos = Allfeetpos[-1].copy()
        pb, res, time = solveMIPGurobi(pb, surfaces, MIP = True, draw_scene = None, plot = True, l1Contact = False, initPos = initPos, endPos = endPos, initCom = initCom, endCom=  endCom)
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
        
    initPos = [array([-0.8, 0.2, 0.13]),
 array([-0.8, -0.2, 0.13]),
 array([-1.6, 0.2, 0.13]),
 array([-1.6, -0.2, 0.13])]
    
    endCom = [10., 0., 0.4]
    
    ax = draw_scene()
    plotPoints(ax, initPos)
    # ~ plt.show()
    pb, coms, footpos, allfeetpos, res = solve(initPos = initPos, endCom=endCom, ax = ax)
    
    
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
    
    
    
    
