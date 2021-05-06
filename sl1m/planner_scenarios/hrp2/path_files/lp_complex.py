from hpp.corbaserver.rbprm.hrp2 import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time
print("Plan guide trajectory ...")
from . import lp_complex_path as tp
print("Done.")
import time

pId = tp.ps.numberPaths() -1
fullBody = Robot ()

# Set the bounds for the root
fullBody.setJointBounds ("root_joint", [-10.135,2, -20, 20, 0, 2.8])
#fullBody.setConstrainedJointsBounds()
#fullBody.setJointBounds('leg_left_1_joint',[-0.1,0.2])
#fullBody.setJointBounds('leg_right_1_joint',[-0.2,0.1])
# add the 6 extraDof for velocity and acceleration (see *_path.py script)
fullBody.client.robot.setDimensionExtraConfigSpace(0)
#~ fullBody.client.robot.setExtraConfigSpaceBounds([-tp.vMax,tp.vMax,-tp.vMax,tp.vMax,0,0,-tp.aMax,tp.aMax,-tp.aMax,tp.aMax,0,0])
ps = tp.ProblemSolver( fullBody )
#~ ps.setParameter("Kinodynamic/velocityBound",tp.vMax)
#~ ps.setParameter("Kinodynamic/accelerationBound",tp.aMax)
#load the viewer
v = tp.Viewer (ps,viewerClient=tp.r.client, displayCoM = True)

# load a reference configuration
#~ q_ref = fullBody.referenceConfig[::]+[0]*6
q_ref = fullBody.referenceConfig[::]
q_init = q_ref[::] 
fullBody.setReferenceConfig(q_ref)
fullBody.setCurrentConfig (q_init)
fullBody.setPostureWeights(fullBody.postureWeights[::])
#~ fullBody.setPostureWeights(fullBody.postureWeights[::])
#fullBody.usePosturalTaskContactCreation(True)

print("Generate limb DB ...")
tStart = time.time()
# generate databases : 

nbSamples = 10000
fullBody.addLimb(fullBody.rLegId,fullBody.rleg,fullBody.rfoot,fullBody.rLegOffset,fullBody.rLegNormal, fullBody.rLegx, fullBody.rLegy, nbSamples, "fixedStep1", 0.01,kinematicConstraintsPath=fullBody.rLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.rLegId, "ReferenceConfiguration", True)
fullBody.addLimb(fullBody.lLegId,fullBody.lleg,fullBody.lfoot,fullBody.lLegOffset,fullBody.rLegNormal, fullBody.lLegx, fullBody.lLegy, nbSamples, "fixedStep1", 0.01,kinematicConstraintsPath=fullBody.lLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.lLegId, "ReferenceConfiguration", True)


tGenerate =  time.time() - tStart
print("Done.")
print("Databases generated in : "+str(tGenerate)+" s")

#define initial and final configurations : 
configSize = fullBody.getConfigSize() -fullBody.client.robot.getDimensionExtraConfigSpace()



#~ q_init[0:7] = tp.ps.configAtParam(pId,0.)[0:7] # use this to get the correct orientation
q_init[:3] = [0.1, -0.82, 0.648702]
q_init[7:] = [ 0.0, 0.0, 0.0, 0.0,                                                  # CHEST HEAD 7-10
        0.261799388,  0.174532925, 0.0, -0.523598776, 0.0, 0.0, 0.17, 		 # LARM       11-17
        0.261799388, -0.174532925, 0.0, -0.523598776, 0.0, 0.0, 0.17, 		 # RARM       18-24
        0.0, 0.0, -0.453785606, 0.872664626, -0.41887902, 0.0,               # LLEG       25-30
        0.0, 0.0, -0.453785606, 0.872664626, -0.41887902, 0.0,               # RLEG       31-36
        ];
v(q_init)

#fullBody.resetJointsBounds()
from hpp.corbaserver.rbprm import  rbprmstate
from hpp.corbaserver.rbprm import  state_alg

def getContactsFromConfig(q, limbs = [Robot.rLegId, Robot.lLegId]):
    s = rbprmstate.State(fullBody, q = q, limbsIncontact = limbs)
    rLegPn = s.getContactPosAndNormalsForLimb(Robot.rLegId)
    lLegPn = s.getContactPosAndNormalsForLimb(Robot.lLegId)
    return { Robot.rLegId : (rLegPn[0][0], rLegPn[1][0]),Robot.lLegId : (lLegPn[0][0], lLegPn[1][0]) }
    
#~ s = getContactsFromConfig ( q = fullBody.getCurrentConfig())


from sl1m.planner_scenarios.complex import *

from hpp.corbaserver.rbprm.rbprmstate import StateHelper

pb, coms, footpos, allfeetpos, res = solve() 

def gen_state(s, pId, com , num_max_sample = 1, first = False, normal = z, newContact = True  ):
    #~ pId = 6
    phase =  pb["phaseData"][pId]
    moving = phase["moving"]
    movingID = fullBody.lLegId
    if moving == RF:
        movingID = fullBody.rLegId
    pos = allfeetpos[pId]; pos[2]+=0.01
    if not first:
        com[2] += 0.5
    #~ com[2] += 0.4
    #~ print "com" , com
    #~ print "pos" , pos.tolist()
    #~ q = fullBody.getCurrentConfig()
    #~ s, succ = state_alg.addNewContact(s, fullBody.rLegId, rfPos.tolist(), z.tolist())
    if newContact:
        sres, succ = state_alg.addNewContact(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= 200)
    else:
        sres, succ = StateHelper.cloneState(s)
    if not succ:
        print("succes ?", succ)
    succCom = False
    ite = 0
    #~ if first:
        #~ print "FIRST "
        #~ com[2] -= 0.25
    #~ while(not succCom and ite < 11):
    while(not succCom and ite < 30):
        succCom =  fullBody.projectStateToCOM(sres.sId ,com, num_max_sample)
        com[2] -= 0.05
        ite += 1    
        if succCom:
            q = sres.q()
            q[3:7] = [0.,0.,0.,1.]
            q[3:] = fullBody.referenceConfig[3:]
            sres.setQ(q)
            succCom =  fullBody.projectStateToCOM(sres.sId ,com, num_max_sample)
            if not succCom:
                print("refail")
    v(sres.q())
    return sres
    
    
    
q = fullBody.getCurrentConfig()
q[:3] = [-2.69, 0.24, 0.649702]
v(q)
s = rbprmstate.State(fullBody, q = q, limbsIncontact = [fullBody.lLegId, fullBody.rLegId])

idfirst = 2
coms[0] = array(s.getCenterOfMass())
#~ coms[0] = coms[1]
#~ coms[1] = array(s.getCenterOfMass())

def normal(phase):
    s = phase["S"][0]
    n = cross(s[:,1] - s[:,0], s[:,2] - s[:,0])
    n /= norm(n)
    if n[2] < 0.:
        for i in range(3):
            n[i] = -n[i]
    #~ print "normal ", n
    return n
    
all_states = [s]
sprev = s
#~ for i in range(0, len(pb["phaseData"])-1):    
for i in range(0, 5):    
#~ for i in range(0, 3):    
#~ for i in range(2, 5):    
#~ for i in range(2, 3):    
#~ for i in range(1, 2):    
    #~ com = (coms[i-i]  ).tolist()
    com = (coms[i]).tolist()
    #~ if i > 2:
        #~ com = (coms[i-1] + (coms[i] - coms[i-1]) *0.8).tolist()
    #get normal
    n = normal(pb["phaseData"][i])
    snew = gen_state(sprev, i+2 , com, num_max_sample = 200, first = False, normal = n )
    all_states += [snew]   
    sprev = snew
    #~ com2 = coms[i+1].tolist()
    #~ snew2 = gen_state(snew, i , com2, num_max_sample = 20, normal = n, newContact = False )
    #~ all_states += [snew2]
    #~ sprev = snew2
    
#~ raise ValueError
    
#~ all_states = all_states[:-1]
extraDof = [0] * 6
configs = [ st.q() + extraDof for st in all_states[:]]; i = 0

print("SID ", [s.sId for s in all_states])

beginId = 0

#~ configs = configs[:-1]

#~ configs = [q1, q2, q3, q4, q5]
#~ from cPickle import dump

#~ f = open('contacts_plateformes.txt','w')
#~ dump(contacts,f)
#~ f.close()
ax = draw_scene()
    #~ plotQPRes(pb, res, ax=ax, plot_constraints=True)
#~ plotQPRes(pb, res, ax=ax, plot_constraints=False)

paths = []

def play_int(f_r = 100):
    for (pid, pl) in paths:
        for i in range(f_r):
            frame = float(i) / float(f_r) * pl
            v((ps.configAtParam(pid,frame)))
            

#~ def ppath(pid):
    #~ length = 

#~ raise ValueError


comTrajs = []
stateTrajs = []

def comTraj(c0, c1, t = 1.):
    zero = [0.,0.,0]
    return fullBody.generateComTraj([c0,c1],[zero],[zero],1.)


from hpp.gepetto import PathPlayer
pp = PathPlayer (v)

from time import sleep
from hpp.corbaserver.rbprm.rbprmstate import  StateHelper

def dq(t = 1, qs = configs):
    for q in qs: 
        v(q); sleep(t)
    
def projCom(state, c, qref):
    s, succ = StateHelper.cloneState(state)
    succCom = succ and fullBody.projectStateToCOM(s.sId ,c, 0)
    if succCom:
        return s
    #~ else:
        #~ c[2] += 0.2
    #~ v(qref)
    for i in range(10):
        s.setQ(qref)
        succCom = succ and fullBody.projectStateToCOM(s.sId ,c, 0)
        c[2] -= .05
        if succCom:
            return s
    #~ else:
    print("fail to project com")
    #~ return state.q()
    return state
        
def projectComPath(state,c0,c1,directPathId, fps = 24.):
    l = ps.pathLength(directPathId)
    f = l/float(fps)
    frame = 0.
    s = state
    qs = []
    ac0 = array(c0)
    ac1 = array(c1)
    #~ print "c0", c0
    #~ print "c1", c1
    while frame < l:
        #~ print "frame", frame
        c = (ac0 + frame * (ac1 - ac0)).tolist()
        #~ print "c", c
        frame += f
        s = projCom(s, c, ps.configAtParam(directPathId,frame))
        qs += [s.q()]
    return qs

def nil():
    global paths
    global stateTrajs
    global comTrajs
    #~ for j in range(2,len(all_states)-2,2):    
    for j in range(0,len(all_states),2):    
    #~ for j in range(4,10,2):    
    #~ for j in range(4,50,2):    
        print("state ", j)
        pathId = fullBody.limbRRT(all_states[j].sId, all_states[j+1].sId)
        #~ print "state ", j
        #~ pathId2 = fullBody.limbRRT(all_states[j+1].sId, all_states[j+2].sId)
        #~ print "state ", j
        paths += [(pathId, ps.pathLength(pathId))]
        #~ print "j", j 
        if j+2 < len(all_states):
            nPid = ps.directPath(all_states[j+1].q(),all_states[j+2].q(),False)[1]
            paths += [(nPid, ps.pathLength(nPid))]
            
            #~ print "coms j", all_states[j].getCenterOfMass()
            #~ print "coms j +1", all_states[j+1].getCenterOfMass()
            #~ print "coms j +2", all_states[j+2].getCenterOfMass()
            
            p0 = comTraj(all_states[j].getCenterOfMass(),all_states[j].getCenterOfMass())
            p1 = comTraj(all_states[j].getCenterOfMass(),all_states[j+1].getCenterOfMass())
            p2 = comTraj(all_states[j+1].getCenterOfMass(),all_states[j+2].getCenterOfMass())
            #~ qs = projectComPath(all_states[j+1],all_states[j+1].getCenterOfMass(), all_states[j+2].getCenterOfMass(),nPid)
            #~ pp(pathId)
            #~ pp(pathId2)
            #~ dq(1./24., qs)
            comTrajs += [(p0,p1,p2)]
            stateTrajs += [(all_states[j].sId, all_states[j+2].sId)]
            #~ comTrajs += [comTraj(all_states[j+1].getCenterOfMass(),all_states[j+2].getCenterOfMass())]
            #~ stateTrajs += [(all_states[j+1].sId, all_states[j+2].sId)]
        #~ if j > 30:
            #~ play_int()    
    #~ for stateTraj, comTtraj in zip(stateTrajs,comTrajs):
        #~ paths+=[fullBody.effectorRRTFromPosBetweenState(stateTraj[0],stateTraj[1],comTtraj[0],comTtraj[1],comTtraj[2])[-1]]
        #~ pp(int(paths[-1]))
        
def comrrt():
    global paths
    global stateTrajs
    global comTrajs
    #~ for j in range(2,len(all_states)-2,2):    
    for j in range(0,len(all_states),2):    
        #~ print "state ", j
        #~ pathId = fullBody.limbRRT(all_states[j].sId, all_states[j+1].sId)
        #~ paths += [(pathId, ps.pathLength(pathId))]
        #~ print "j", j 
        if j+2 < len(all_states):
            #~ nPid = ps.directPath(all_states[j+1].q(),all_states[j+2].q(),False)[1]
            #~ paths += [(nPid, ps.pathLength(nPid))]
            
            p0 = comTraj(all_states[j].getCenterOfMass(),all_states[j].getCenterOfMass())
            p1 = comTraj(all_states[j].getCenterOfMass(),all_states[j+1].getCenterOfMass())
            p2 = comTraj(all_states[j+1].getCenterOfMass(),all_states[j+2].getCenterOfMass())
            comTrajs += [(p0,p1,p2)]
            stateTrajs += [(all_states[j].sId, all_states[j+2].sId)]
            #~ comTrajs += [comTraj(all_states[j+1].getCenterOfMass(),all_states[j+2].getCenterOfMass())]
            #~ stateTrajs += [(all_states[j+1].sId, all_states[j+2].sId)]
        #~ play_int()    
    for stateTraj, comTtraj in zip(stateTrajs,comTrajs):
        paths+=[fullBody.comRRTFromPosBetweenState(stateTraj[0],stateTraj[1],comTtraj[0],comTtraj[1],comTtraj[2])[-1]]
        pp(int(paths[-1]))
        
#~ nil()

#~ comrrt()
    
#~ play_int()
#~ nil()
#~ beginId = 0
