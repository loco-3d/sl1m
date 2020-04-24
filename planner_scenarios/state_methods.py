from hpp.gepetto import Viewer
from hpp.corbaserver import Client
from hpp.corbaserver.problem_solver import ProblemSolver
from hpp.corbaserver.rbprm import  rbprmstate
from hpp.corbaserver.rbprm import  state_alg
from hpp.corbaserver.rbprm.rbprmstate import StateHelper
from tools.display_tools import *

from scipy.optimize import linprog

from numpy import array, zeros, ones
import time

z = array([0.,0.,1.])

# ~ def getContactsFromConfig(q, limbs = Robot.limbs_names):
def getContactsFromConfig(fullBody, q, limbs):
    s = rbprmstate.State(fullBody, q = q, limbsIncontact = limbs)
    res = {}
    for limb in limbs:
        rLegPn = s.getContactPosAndNormalsForLimb(limb)
        res[limb] = (rLegPn[0][0], rLegPn[1][0])
    return res

def initConstraintsFrom_q_init(fullBody, q_init, limbNames):
    initContacts = getContactsFromConfig(fullBody, q_init, limbNames)
    initPos = [array(initContacts[el][0][0]) for el in limbNames]
    initCom = array(fullBody.getCenterOfMass())
    return initContacts, initPos, initCom

#static eq is com is convex combination of pos (projected)
def staticEq(positions, com):
    sizeX = len(positions)
    E = zeros((3,sizeX))
    for i, pos in enumerate(positions):
        E[:2,i] = pos[:2]
    e = array([com[0], com[1], 1.])
    E[2,:] = ones(sizeX)
    try:
        res = linprog(ones(sizeX), A_ub=None, b_ub=None, A_eq=E, b_eq=e, bounds=[(0.,1.) for _ in range(sizeX)], method='interior-point', callback=None, options={'presolve': True})
        return res['success']
    except:
        return False

def gen_state(pb, coms, allfeetpos,fullBody, limbNames, s, pId,  num_max_sample = 1, first = False, normal = z, newContact = True , useCom = False ):
    #~ pId = 6
    feetPos = allfeetpos[pId]
    sres, succ = StateHelper.cloneState(s)
    for limbId, pos in zip(limbNames, feetPos):
        sres, succ = state_alg.addNewContact(sres, limbId, pos.tolist(), normal.tolist(), num_max_sample= 200)         
        # ~ v(sres.q())   
        if not succ:
            print("succes ?", succ)
    succCom = False
    ite = 0
    # ~ #~ if first:
        # ~ #~ print "FIRST "
        # ~ #~ com[2] -= 0.25
    # ~ #~ while(not succCom and ite < 11):
    com = coms[pId*2].tolist()
    # ~ print("com ", com)
    #check whether com is feasible
    if staticEq(feetPos, com):
        print ("Config equilibrium")
        succCom = True
    if useCom:
        while(not succCom and ite < 30):
            succCom =  fullBody.projectStateToCOM(sres.sId ,com, num_max_sample)
            com[2] += 0.05
            ite += 1    
            if succCom:
                sres2, succ = StateHelper.cloneState(sres)
                q = sres2.q()
                q[3:7] = [0.,0.,0.,1.]
                q[3:] = fullBody.referenceConfig[3:]
                sres2.setQ(q)
                succCom2 =  fullBody.projectStateToCOM(sres2.sId ,com, num_max_sample)
                if succCom2:
                    print ("SUCCESS COM2 ", succCom)
                    sres = sres2
        print ("SUCCESS COM ", succCom)
    # ~ v(sres.q())
    fullBody.setCurrentConfig(sres.q())
    return sres

def applyOffset(allfeetpos, offset):
    for i in range(len(allfeetpos)):
        for j in range(len(allfeetpos[i])):
            allfeetpos[i][j][2]+=offset

#### show results #### 

def run(fullBody, states, allfeetpos, limbNames, coms, pb, offset, useCom = False):
    applyOffset(allfeetpos, offset)
    for i, _ in enumerate(allfeetpos[1:]):
        states.append(gen_state(pb,coms,allfeetpos,fullBody, limbNames, states[-1],i+1, useCom = useCom))
        
def play(states, v, dt = 0.5):
    for s in states:
        v(s.q())
        time.sleep(dt)
     
#### export contact sequence #### 
        
#now saving sequence to mlp
from pinocchio.utils import *
import importlib
import multicontact_api
from multicontact_api import ContactSequence
from mlp.utils.cs_tools import addPhaseFromConfig, setFinalState
from mlp.viewer.display_tools import initScene, displaySteppingStones
from pinocchio.utils import matrixToRpy
from pinocchio import Quaternion, SE3
# ~ from hpp.corbaserver.rbprm.tools.surfaces_from_path import getSurfacesFromGuideContinuous
import random
from mlp.utils.requirements import Requirements
multicontact_api.switchToNumpyArray()
from numpy.linalg import  norm


extraDof = [0 for _ in range(6)]


def exportCS(fullBody, q_init, states, pb, allfeetpos, limbNames, effectorNames):
    configs = [ st.q() + extraDof for st in states[:]]; i = 0
    cs = ContactSequence(0)
    addPhaseFromConfig(fullBody, cs, q_init, limbNames[:])
    rot = Quaternion.Identity() #todo update
    for pId in range(1, len(states)):
        print ('phase ', pId)
        # ~ rot = Quaternion.Identity()
        prevFeetPos = allfeetpos[pId-1]
        feetPos = allfeetpos[pId]
        switch = False
        for limbId, effId, pos, prevPos in zip(limbNames, effectorNames,  feetPos, prevFeetPos):
            if norm(pos - prevPos) > 0.001:
                if switch:
                    print ("already one change, error")
                switch = True
                placement = SE3()
                placement.translation = np.array(pos).T
                placement.rotation = rot.matrix()
                cs.moveEffectorToPlacement(effId, placement)  
        if not switch:
            print ("no switch at phase")
        q_end = configs[-1][:]
        fullBody.setCurrentConfig(q_end[:-6])
        com = fullBody.getCenterOfMass()
        setFinalState(cs, com, q=q_end)
    return cs
