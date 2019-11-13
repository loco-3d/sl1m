import numpy as np
print "Plan guide trajectory ..."
import lp_slalom_debris_path as tp
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
Z_AXIS = np.array([0,0,1]).T
from sl1m.planner_scenarios.talos.constraints import *
"""
from sl1m.plot_plytopes import *
rubB1 =   [[-0.5079353885643417, 1.282387089736744, 0.0021999916093061456],  [-0.5079353885643416, 1.2444260597154655, 0.0021999916093061456],  [-0.4555532842024058, 1.244426054175804, 0.016235816474940087],  [-0.4555533048779398, 1.2823870841970821, 0.016235812378438574]]
rubB2 =   [[-0.49591913077286087, 1.1089893205877919, 0.004305334067727976],  [-0.4959192866539029, 1.0723211830469714, 0.014130502161835578],  [-0.4416887364162, 1.0723216355962044, 0.014130432868663537],  [-0.4416885369888308, 1.1089893095084922, 0.004305337036416214]]
rubB3 =   [[-0.48205335101433644, 0.8487340468277829, 0.016235827387680896],  [-0.42967124312103594, 0.8487340681209311, 0.00219999743865322],  [-0.4296712416379183, 0.8866951153800525, 0.002199995597725976],  [-0.4820533288556889, 0.886695094086905, 0.016235821450252305]]
rubB4 =   [[-0.4959189669199701, 0.6631073957571295, 0.01413046086543412],  [-0.49591898907856746, 0.6264397238524806, 0.004305370958339671],  [-0.4416886194638495, 0.6264397496095392, 0.004305374972883749],  [-0.4416885973052521, 0.6631074215141881, 0.014130464879978225]]
rubB5 =   [[-0.1286227022905665, 1.2865674070731667, 0.004305356276006811],  [-0.18285298595645427, 1.2865674159019864, 0.004305353910326992],  [-0.18285296418321512, 1.2498997897345068, 0.01413043425581112],  [-0.1286226893627057, 1.2498997809056869, 0.014130436621490938]]
rubB6 =   [[-0.19487018593892436, 1.0642733322579028, 0.002199924194792202],  [-0.19487018593892436, 1.026311899523103, 0.002199924194792202],  [-0.14248738910025185, 1.0263118953186507, 0.016235892193181743],  [-0.14248740407073668, 1.06427332805345 , 0.01623588639038652]]
rubB7 =   [[-0.1286224879028442, 0.8406861284205095, 0.014130488309530481],  [-0.18285313142627807, 0.8406861372493323, 0.014130490675199364],  [-0.18285310965307405, 0.8040180965091135, 0.004305348131801057], [-0.12862247497500423, 0.8040180876802904, 0.00430534576613209]]
rubB8 =   [[-0.16898820573177617, 0.6685813451790906, 0.016235902651238112],  [-0.16898819670073095, 0.6306198870269437, 0.016235902925035444],  [-0.11660536431136716, 0.6306198843501182, 0.002199925819797003],  [-0.11660538815137254, 0.6685813156694488, 0.002199931316793472]]

rubE1 =  [[0.47039856640115135, 1.2823872614660115, 0.016235877922373423],  [0.4180158210967855, 1.282387267005657, 0.002199926065377139],  [0.41801584177237405, 1.2444258935258445, 0.0021999301618671058],  [0.47039856640115135, 1.244425887986198, 0.016235877922373423]]
rubE2 =  [[0.4300330225380426, 1.1089893042736412, 0.004305344981069587],  [0.4300328231105139, 1.0723212368072526, 0.014130497185884688],  [0.4842633076431641, 1.0723216893563094, 0.014130427892686833],  [0.48426346352412575, 1.1089892931943353, 0.004305347949760364]]
rubE3 =  [[0.44389771410174234, 0.8866951647708724, 0.01623582554675447],  [0.4438977126186285, 0.8487341064324276, 0.01623582738768091],  [0.496279820511929, 0.8487341277255759, 0.002199997438653234],  [0.4962798426705747, 0.8866951860640212, 0.0021999915012257676]]
rubE4 =  [[0.43003313137984017, 0.6264398056157277, 0.00430537689572764],  [0.4842634139015866, 0.6264398313727851, 0.00430538091027495],  [0.48426343606019295, 0.6631074589602306, 0.014130458942590155],  [0.4300331535384464, 0.6631074332031729, 0.014130454928042817]]
rubE5 =  [[0.743098096698756, 1.28656744542173, 0.004305359084522087],  [0.7430980523815441, 1.249899816830634, 0.014130440079408849],  [0.7973283363876394, 1.2498998683447495, 0.01413043205031209],  [0.7973283807048511, 1.2865674969358454, 0.004305351055425384]]
rubE6 =  [[0.783463669083787, 1.0642734072201385, 0.016235889639458684],  [0.7310808794622428, 1.0642733764552545, 0.0021999259076780076],  [0.7310808986548348, 1.0263119698608123, 0.002199928163265169],  [0.7834636903246144, 1.0263119469600557, 0.01623589513753257]]
rubE7 =  [[0.7430979576180856, 0.8040181510934469, 0.004305346786616704],  [0.7973285946775992, 0.8040181510934469, 0.004305346786616704],  [0.7973285946775991, 0.8406861907949811, 0.014130489051701167],  [0.7430979576180857, 0.8406861907949811, 0.014130489051701167]]
rubE8 =  [[0.7569629079334563, 0.6685813931634638, 0.01623589835472744],  [0.7569629169645031, 0.6306199582518592, 0.016235898628522177],  [0.8093456996032325, 0.6306199589940333, 0.0021999313804602627],  [0.8093456885239412, 0.6685813670728219, 0.0021999343491562184]]




mid =  [[0.2525372748550794, 1.188345286893716, 0.019086446287424022],  [0.0617591727940063, 1.1883452947119464, 0.019086446287423415],  [0.06175921867215623, 0.7002973822902399, 0.019086484191625633],  [0.25253731563565707, 0.7002973744720097, 0.01908648419162624]]
beginFloor = [[-1.7 ,0.2 ,0.],[ -0.8,1.2 ,0.],[-0.8,0.5 ,0.],[-1.7 ,-0.4 ,0.]]
endFloor = [[ 1.,1.2 ,0.],[2 ,0.2 ,0.],[2 ,-0.4 ,0.],[1.,0.5 ,0.]]
final = [[2 ,0.1 ,0.],[2.1 ,0.1 ,0.],[2.1 ,-0.1 ,0.],[2 ,-0.1 ,0.]]

all_surfaces = [beginFloor,rubB1,rubB2,rubB3,rubB4,rubB5,rubB6,rubB7,rubB8,mid,rubE1,rubE2,rubE3,rubE4,rubE5,rubE6,rubE7,rubE8,endFloor]

arubB1  = array(rubB1).T
arubB2  = array(rubB2).T
arubB3  = array(rubB3).T
arubB4  = array(rubB4).T
arubB5  = array(rubB5).T
arubB6  = array(rubB6).T
arubB7  = array(rubB7).T
arubB8  = array(rubB8).T
arubE1  = array(rubE1).T
arubE2  = array(rubE2).T
arubE3  = array(rubE3).T
arubE4  = array(rubE4).T
arubE5  = array(rubE5).T
arubE6  = array(rubE6).T
arubE7  = array(rubE7).T
arubE8  = array(rubE8).T
amid  = array(mid).T
abeginFloor  = array(beginFloor).T
aendFloor  = array(endFloor).T


all_rub_begin = [arubB1,arubB2,arubB3,arubB4,arubB5,arubB6,arubB7,arubB8]
all_rub_and = [arubE1,arubE2,arubE3,arubE4,arubE5,arubE6,arubE7,arubE8]


surfaces = [[abeginFloor],[arubB3],[arubB6],[amid],[amid],[arubE3],[arubE6],[aendFloor],[aendFloor]] # work with p0x  -0.80


#~ surfaces = [[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[arubB2],[arubB7],[amid],[arubE3],[arubE6],[aendFloor],[aendFloor]] 
surfaces = [[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[abeginFloor],[arubB2],[arubB7],[amid],[arubE3],[arubE6],[aendFloor],[aendFloor]]
"""

final = [[2 ,0.15 ,0.],[2.1 ,0.15 ,0.],[2.1 ,-0.15 ,0.],[2 ,-0.15 ,0.]]
afinal = array(final).T

def gen_pb(root_init,R, surfaces):
    
    nphases = len(surfaces)
    lf_0 = array(root_init[0:3]) + array([0, 0.085,-0.98]) # values for talos ! 
    rf_0 = array(root_init[0:3]) + array([0,-0.085,-0.98]) # values for talos ! 
    #p0 = [array([-3.0805096486250154, 0.335, 0.]), array([-3.0805096486250154, 0.145,0.])];  ## FIXME : get it from planning too
    #p0 = [array([-0.1805096486250154, 0.335, 0.]), array([-0.1805096486250154, 0.145,0.])];  ## FIXME : get it from planning too
    p0 = [lf_0,rf_0];
    print "p0 used : ",p0
    
    res = { "p0" : p0, "c0" : None, "nphases": nphases}
    #res = { "p0" : None, "c0" : None, "nphases": nphases}
    
    #print "surfaces = ",surfaces
    #TODO in non planar cases, K must be rotated
    #phaseData = [ {"moving" : i%2, "fixed" : (i+1) % 2 , "K" : [copyKin(kinematicConstraints) for _ in range(len(surfaces[i]))], "relativeK" : [relativeConstraints[(i)%2] for _ in range(len(surfaces[i]))], "S" : surfaces[i] } for i in range(nphases)]
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
    defaultStep = 1.
    step = defaultStep
    variation = 0.2
    while not success and it < maxIt:
      if it > 0 :
        step = defaultStep + random.uniform(-variation,variation)
      R,surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder,tp.ps,tp.afftool,tp.pathId,None,step,True)
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
    from sl1m.fix_sparsity import solveL1
    step = 1.
    R,surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder,tp.ps,tp.afftool,tp.pathId,None,step,True)

    pb = gen_pb(tp.q_init,R,surfaces)

    pb, coms, footpos, allfeetpos, res = solveL1(pb, surfaces, draw_scene)

