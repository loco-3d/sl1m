import numpy as np


from constants_and_tools import *

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *
#~ from qp import solve_lp

import eigenpy
eigenpy.switchToNumpyArray()
from hpp_centroidal_dynamics import *
from hpp_spline import *
import cdd
from curves import bezier3
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray
from cwc_qhull import compute_CWC
import qp

eigenpy.switchToNumpyArray()



def ineqFromCdata(cData):
    (H, h) = cData.contactPhase_.getPolytopeInequalities()
    return (-H, h)


def quasi_flat(P):
    aPmin = [ array(p[:2] + [-1])  for p in P]
    aPmax = [ array(p[:2] + [3 ])  for p in P]
    hull = ConvexHull(array(aPmin + aPmax))
    return ineqQHull(hull)

def staticEquilirium3DPolytope(cData, P, addKin = True, qflat = False):
    m = cData.contactPhase_.getMass()
    #~ #get 6D polytope
    A = None
    b = None
    if qflat:
        (A, b) = quasi_flat(P)
    else:
        (H, h) = ineqFromCdata(cData)
    #project to the space where aceleration is 0
        D = zeros((6,3))
        D[3:,:] = m * gX
        
        d = zeros((6,))
        d[:3] = -m * g
        
        A = H.dot(D)
        b = h.reshape((-1,)) - H.dot(d)    
    #add kinematic polytope
    if addKin :
        (K,k) = (cData.Kin_[0], cData.Kin_[1].reshape(-1,))
        
        A = vstack([A, K])
        b = concatenate([b, k]).reshape((-1,))
    
    return (A, b)


def centroid(pts):
    return sum([array(pt) for pt in pts]) / float(len(pts))

def findIntermediateContactPhase(phase1, phase2):
    keys = phase1.keys() 
    assert phase1.keys() == phase2.keys()
    res = {}
    for k,v in phase1.iteritems():
        if norm(centroid(v[0]) - centroid(phase2[k][0])) < 0.01:
            res[k] = v
    return res
            

def genIntermediateContactPhases(contacts):
    res = []
    for i, contact in enumerate(contacts[:-1]):
        nContact = contacts[i+1]
        res += [contact, findIntermediateContactPhase(contact,nContact)]
    res += [contacts[-1]]
    return res
        
from hpp.corbaserver.rbprm.hrp2 import Robot
from cPickle import load
def loadContactPhases(filename='./data/contacts_flat.txt'):    
    f = open(filename,'r')
    contacts = load(f)
    f.close()
    #TODO FIX FIRST PHASE
    contacts[0][Robot.lLegId] = contacts[1][Robot.lLegId]
    allContacts = genIntermediateContactPhases(contacts)
    return allContacts


#~ def ContactPhaseFromData2PAC(contacts, data2PAC, addKin = True):

def ContactPhase(contacts, phaseId, pD = None, addToPlan = True, addKin = True):    
    contact = contacts[phaseId]

    P = []
    N = []
    Kin = zeros((0,3))
    kin = zeros(0)
        
    if contact.has_key(rLeg):
        PR = contact[rLeg][0]
        NR = contact[rLeg][1]
        P += PR
        N += NR
        tr = default_transform_from_pos_normal(centroid(PR),NR[0])
        KinR,kinR = right_foot_hrp2_constraints(tr)
        Kin =  vstack([Kin,KinR])
        kin =  concatenate([kin, kinR])
        
    if contact.has_key(lLeg):
        PL = contact[lLeg][0]
        NL = contact[lLeg][1]
        P += PL
        N += NL
        tl = default_transform_from_pos_normal(centroid(PL),NL[0])
        KinL,kinL = left_foot_hrp2_constraints(tl)
        Kin =  vstack([Kin,KinL])
        kin =  concatenate([kin, kinL])
    
    cData0 = ContactData(Equilibrium("test", 55.88363633, 4, SOLVER_LP_QPOASES, True, 10, False))
    cData0.contactPhase_.setNewContacts(array(P), array(N), mu, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
    if addKin:
        cData0.setKinematicConstraints(matrix(Kin), matrix(kin).T)
    if addToPlan:
        pD.addContact(cData0)
    return cData0, P, N, Kin, kin
    
    
def genConstraints(phaseFrom, phaseTo):
    A0,b0 = phaseFrom["staticCone"]
    A1,b1 = phaseTo["staticCone"]
    K0, k0 = phaseFrom["Kin"]
    K1, k1 = phaseTo["Kin"]
    return (vstack([A0, A1,K0,K1]), concatenate([b0, b1, k0, k1]).reshape((-1,)) )
   

def normalize(A,b):
    Ares = zeros(A.shape)
    bres = zeros(b.shape)
    for i in range(A.shape[0]):
        n = norm(A[i,:])
        if n <= 0.000001:
            n = 1.
        Ares[i,:] = A[i,:] / n; bres[i] = b[i] / n
    print Ares
    return Ares, bres
    
def genConstraintsRobust(phaseFrom, phaseTo):
    A0,b0 = phaseFrom["staticCone"]
    A1,b1 = phaseTo["staticCone"]
    K0, k0 = phaseFrom["Kin"]
    K1, k1 = phaseTo["Kin"]
    A = vstack([A0, A1, K0, K1])
    #~ print "K0.shape ", K0.shape
    #~ print "k1.shape ", k1.shape
    #~ print "A.shape ", A.shape
    A = hstack([A,ones((A.shape[0],1))])
    #~ print "A ", A
    posConstraint = zeros((1,4)); posConstraint[0,-1] = -1
    A = vstack([A, posConstraint])
    b = concatenate([b0, b1, k0, k1, zeros(1)]).reshape((-1,))
    return normalize(A,b )
    
from itertools import izip

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)
    
    
def gen2PacSolution(contacts, startCOM = [0.01247679248626212, 0.0015769611415203033, 0.8127583093263778]):
    resData = [{} for _ in contacts]
    #for each contact phase
    # dic with static cone, global cone, kinematic constraints, 2 pac solution
    # TODO: find minimum cost over all points
    for phaseId, contact in enumerate(contacts):
        #~ resData[i][] =
        cData, P, N, K, k = ContactPhase(contacts, phaseId, addToPlan = False, addKin = True)
        H, h = ineqFromCdata(cData)
        A,b = staticEquilirium3DPolytope(cData, P, addKin = True, qflat = False)
        resData[phaseId] = {"cData" : cData, "P" : P, "N" : N, "Kin" : (K,k), "cone6D" : (H, h), "staticCone": (A,b) }
        
    resData[0]["waypoints"] = array(startCOM)
    #init first point with initial solution
    for i, data0, in enumerate(resData[:-1]):
        data1 = resData[i+1]
        #cost in minimum distance
        #~ A,b = genConstraints(data0, data1)
        #~ costTarget = data0["waypoints"]
        #~ res = qp.solve_least_square(A = identity(3),b = data0["waypoints"] ,G=A, h=b, C=None, d=None)
        # cost is maximize safety
        A,b = genConstraintsRobust(data0, data1)
        lcost = zeros(4); lcost[-1] = -1
        res = qp.solve_lp(lcost ,G=A, h=b, C=None, d=None)
        #~ print "b ", data0["waypoints"]
        print "res ", res
        resData[i+1]["waypoints"] = res[:3]
        
    return resData
    

#all colors
#~ from random import randint
#~ numColors = 20
#~ colors = [str(float(i) / float(numColors)) for i in range(numColors)]


def crocpb_from_quasi_static(data2PAC, phase1, phase2, flag = None):
    pD = ProblemData()
    for i in range(phase1,phase2 +1):
        #TODO avoid recomputation of polytopes..
        #~ ContactPhase(contacts, i, pD, addKin = False)
        #~ ContactPhase(contacts, i, pD, addKin = True)
        pD.addContact(data2PAC[i]["cData"])
    
    #~ c = Constraints()
    #~ c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL | ConstraintFlag.END_VEL | ConstraintFlag.END_POS
    #~ c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL 
    #~ pD.constraints_.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL |   ConstraintFlag.END_VEL | ConstraintFlag.END_POS 
    if flag is None:
        pD.constraints_.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL |   ConstraintFlag.END_VEL | ConstraintFlag.END_POS 
    else:
        pD.constraints_.flag_ = flag
    #~ pD.constraints_.flag_ = ConstraintFlag.INIT_POS| ConstraintFlag.INIT_VEL |ConstraintFlag.END_POS
    pD.c0_ = array(data2PAC[phase1]["waypoints"]).reshape((-1,1))
    pD.c1_ = array(data2PAC[phase2]["waypoints"]).reshape((-1,1))
    return pD
    
def crocpb_from_bezier_and_quasi_static(data2PAC, phase1, phase2, b):
    pD = ProblemData()
    for i in range(phase1,phase2 +1):
        #TODO avoid recomputation of polytopes..
        #~ ContactPhase(contacts, i, pD, addKin = False)
        #~ ContactPhase(contacts, i, pD, addKin = True)        
        pD.addContact(data2PAC[i]["cData"])
    
    #~ c = Constraints()
    c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL | ConstraintFlag.END_VEL | ConstraintFlag.END_POS
    #~ c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL 
    #~ pD.constraints_.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL |  ConstraintFlag.INIT_ACC | ConstraintFlag.END_VEL | ConstraintFlag.END_POS | ConstraintFlag.END_ACC
    
    #~ pD.constraints_.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL |   ConstraintFlag.END_VEL | ConstraintFlag.END_POS 
    #~ pD.constraints_.flag_ = ConstraintFlag.INIT_POS| ConstraintFlag.INIT_VEL |ConstraintFlag.END_POS
    pD.c0_  = b.derivate(0.,0)
    pD.dc0_ = b.derivate(0.,1)
    pD.ddc0_ = b.derivate(0.,1)
    pD.c1_  = array(data2PAC[phase2]["waypoints"]).reshape((-1,1))
    return pD
  
    
def plotBezier(bez, ax, color = "b", label = None, linewidth = 2.0, D3 = True, mx = None):
        step = 1000.
        if mx is None:
                mx = bez.max()
        #~ points1 =  np.array([(bez(i/step*mx)[0][0],bez(i/step*mx)[1][0],bez(i/step*mx)[2][0]) for i in range(int(step))])
        points1 =  np.array([(bez(i/step*mx)[0],bez(i/step*mx)[1],bez(i/step*mx)[2]) for i in range(int(step))])
        #~ print "? ", bez(bez.max())
        x = points1[:,0]
        y = points1[:,1]
        if(D3):                
                z = points1[:,2]
                z = [el[0] for el in z]
                #~ print "bezer"
                ax.plot(x.tolist(),y.tolist(),z,color)
                #~ plt.show()
        else:
                ax.plot(x.tolist(),y.tolist() ,color,linewidth=linewidth, label=label)
                #~ plt.show()
        return points1
    
def plotPoints(ax, wps, color = "b", D3 = True, linewidth=2):
    x = array(wps)[:,0]
    y = array(wps)[:,1]
    if(D3):                
            z = array(wps)[:,2]
            ax.scatter(x, y, z, c=color, marker='o', linewidth = 5) 
    else:
            ax.scatter(x,y,color=color, linewidth = linewidth)  
    
def plot(data2PAC, color ="r", D3 = True, linewidth=2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    #~ ax.set_xlim3d(-0, 1.7)
    #~ ax.set_ylim3d(-0.2, 0.2)
    ax.set_zlim3d(0, 1)
    coms = [data["waypoints"] for data in data2PAC]
    plotPoints(ax, coms, color = "r")
    
    for data in data2PAC:
        plotPoints(ax, data["P"], color = np.random.rand(3,))
    #~ cPoints = [] 
    #~ [ cPoints.append(p) for data in data2PAC for p in  data["P"] ]
    #~ for data0, data1 in pairwise(data2PAC):
    #~ plotPoints (ax, cPoints)
    #~ return False
    cx = [c[0] for c in coms]
    cy = [c[1] for c in coms]
    cz = [c[2] for c in coms]
    ax.plot(cx, cy, cz)
    return ax
    #~ plt.show()

if __name__ == '__main__':
    
    
    import time
    contacts = loadContactPhases()
    from hpp.corbaserver.rbprm.hrp2 import Robot
    
    from plot_plytopes import *
    rLeg = Robot.rLegId
    lLeg = Robot.lLegId
    #contact data from platform scenario


    data2PAC = gen2PacSolution(contacts)

    #TODO make sure constraints associated to right leg and orientation

    pD = ProblemData()

    phaseId = 0
    cData0, P0, N0, K0, k0 = ContactPhase(contacts, phaseId, pD)
        
    phaseId = 1
    cData1, P1, N1, K1, k1 = ContactPhase(contacts, phaseId, pD)
    
    phaseId = 2
    cData2, P2, N2, K2, k2 = ContactPhase(contacts, phaseId, pD)
    
        
    c = Constraints()
    #~ c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL | ConstraintFlag.END_VEL
    c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL 
    pD.constraints_.flag_
    
    startCOM = [0.1, 0.25000237362450967, 0.7522584665645766]
    pD.c0_ = matrix(startCOM).T

start = 0

pD = crocpb_from_quasi_static(data2PAC,start,start+2)


allbezier = []
startConstraints = None
rg = 20.
rg2 = rg / 2.
times2 = [matrix([float(i) / rg2 ,float(j) / rg2]).T for i in range(2,int(rg)) for j in range(2,int(rg))]
times3 = [matrix([float(i) / rg2 ,float(j) / rg2, float(k) / rg2]).T for i in range(2,int(rg)) for j in range(2,int(rg)) for k in range(2,int(rg))]
times4 = [matrix([float(i) / rg2 ,float(j) / rg2, float(k) / rg2, float(m) / rg2]).T for i in range(2,int(rg)) for j in range(2,int(rg)) for k in range(2,int(rg)) for m in range(3,int(rg))]
times4.reverse()
#~ times3.reverse()
res = None
#~ for i, el in enumerate(times4[30:]):
for i, el in enumerate(times3[30:]):
    #~ print i
    res = computeCOMTraj(pD, el, -1); 
    if res.success:
        b = res.c_of_t
        allbezier+= [b.extract(b.min(),el[-2,0])]
        startConstraints = b.extract(el[-2,0], b.max())
        print 'success times', el
        print 'time index ', i
        break
        
ax  = plot(data2PAC[:], color ="r", D3 = True, linewidth=2)
        

plotBezier(res.c_of_t, ax, color = "g", label = None, linewidth = 2.0, D3 = True, mx = None)

from time import  clock

t1 = clock()
k = 0
for i in range(start+2, len(data2PAC) - 3):
    pD = crocpb_from_bezier_and_quasi_static(data2PAC,i,i+2, startConstraints)
    res = None
    for j, el in enumerate(times3[:]):
        k += 1 
        res = computeCOMTraj(pD, el, -1); 
        if res.success:
            b = res.c_of_t
            allbezier+= [b.extract(b.min(),el[-2,0])]
            startConstraints = b.extract(el[-2,0], b.max())
            #~ print 'success times', el
            #~ print 'time index ', j
            break
    if not res.success:
        print "fail"
        break
    else:    
        plotBezier(startConstraints, ax, color = "b", label = None, linewidth = 2.0, D3 = True, mx = None)
     
t2 = clock()
print "time per call", (t2 - t1) / float(k) * 1000.
print "total time (ms)",  (t2 - t1) * 1000.

[ plotBezier(b, ax, color = "r", label = None, linewidth = 2.0, D3 = True, mx = None) for b in allbezier]
            
plt.show()


#~ from time import  clock
#~ pD = crocpb_from_bezier_and_quasi_static(data2PAC,2,4, startConstraints)
#~ times4.reverse()
#~ res = None
#~ t1 = clock()
#~ for i, el in enumerate(times3):
    #~ res = computeCOMTraj(pD, el, -1); 
    #~ if res.success:
        #~ t2 = clock()
        #~ print "time per call", (t2 - t1) / float(i) * 1000.
        #~ print 'success times', el
        #~ print 'time index ', i
        #~ break
        
