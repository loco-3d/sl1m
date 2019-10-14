from hpp.corbaserver.rbprm.hrp2 import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time
print "Plan guide trajectory ..."
import lp_complex_path as tp
print "Done."
import time
DEFAULT_COM_HEIGHT = 0.72
pId = tp.ps.numberPaths() -1
fullBody = Robot ()

# Set the bounds for the root
fullBody.setJointBounds ("root_joint", [-10.135,2, -20, 20, 0, 2.8])
fullBody.client.robot.setDimensionExtraConfigSpace(6)
fullBody.client.robot.setExtraConfigSpaceBounds([0,0,0,0,0,0,0,0,0,0,0,0])

ps = tp.ProblemSolver( fullBody )

#load the viewer
v = tp.Viewer (ps,viewerClient=tp.r.client, displayCoM = True)

# load a reference configuration

q_ref = fullBody.referenceConfig[::] + [0]*6
q_init = q_ref[::] 
fullBody.setReferenceConfig(q_ref)
fullBody.setCurrentConfig (q_init)
fullBody.setPostureWeights(fullBody.postureWeights[::] + [0]*6)
#fullBody.usePosturalTaskContactCreation(True)

print "Generate limb DB ..."
tStart = time.time()
# generate databases : 

nbSamples = 10000
fullBody.addLimb(fullBody.rLegId,fullBody.rleg,fullBody.rfoot,fullBody.rLegOffset,fullBody.rLegNormal, fullBody.rLegx, fullBody.rLegy, nbSamples, "fixedStep1", 0.01,kinematicConstraintsPath=fullBody.rLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.rLegId, "ReferenceConfiguration", True)
fullBody.addLimb(fullBody.lLegId,fullBody.lleg,fullBody.lfoot,fullBody.lLegOffset,fullBody.rLegNormal, fullBody.lLegx, fullBody.lLegy, nbSamples, "fixedStep1", 0.01,kinematicConstraintsPath=fullBody.lLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.lLegId, "ReferenceConfiguration", True)


tGenerate =  time.time() - tStart
print "Done."
print "Databases generated in : "+str(tGenerate)+" s"

#define initial and final configurations : 
configSize = fullBody.getConfigSize() -fullBody.client.robot.getDimensionExtraConfigSpace()



#~ q_init[0:7] = tp.ps.configAtParam(pId,0.)[0:7] # use this to get the correct orientation
q_init[:3] = [0.1, -0.82, 0.648702]
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


from mpcroc.planner_scenarios.complex import *

from hpp.corbaserver.rbprm.rbprmstate import StateHelper

pb, coms, footpos, allfeetpos, res = solve() 


def computeCenterOfSupportPolygon(s):
    com = np.zeros(3)
    numContacts = float(len(s.getLimbsInContact()))
    for limbId in s.getLimbsInContact():
        com += np.array(s.getCenterOfContactForLimb(limbId)[0])
    com /= numContacts
    com[2] += DEFAULT_COM_HEIGHT 
    return com.tolist()

def projectCoMInSupportPolygon(s):
    desiredCOM = computeCenterOfSupportPolygon(s)
    #print "try to project state to com position : ",desiredCOM
    success = False
    maxIt = 20
    while not success and maxIt > 0:
      success = s.fullBody.projectStateToCOM(s.sId ,desiredCOM, maxNumSample = 0)
      maxIt -= 1 
      desiredCOM[2] -= 0.005
    return success
      

def gen_state(s, pId, com , num_max_sample = 0, first = False, normal = z, newContact = True  ,projectCOM = True):
    #~ pId = 6
    phase =  pb["phaseData"][pId]
    moving = phase["moving"]
    movingID = fullBody.lLegId
    if moving == RF:
        movingID = fullBody.rLegId
    print "# gen state for phase Id = ",pId
    #print "current config q=",s.q()
    #print "move limb ",movingID
    pos = allfeetpos[pId+2]; # +2 because it contains also the 2 feet pos at the init config
    pos[2] += 0.01
    if newContact:
        sres, succ = StateHelper.addNewContact(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= num_max_sample)
    else:
        sres, succ = StateHelper.cloneState(s)
    if not succ:
      print "Cannot project config q = ",sres.q()
      print "To new contact position for "+movingID+" = "+str(pos.tolist())+" : n = "+str(normal.tolist())
      raise RuntimeError("Cannot project feet to new contact position") # try something else ?? 
    if projectCOM :
        #print "config before projecting to com q1=",sres.q()
        successCOM = projectCoMInSupportPolygon(sres)
        if not successCOM:
            # is it really an issue ? 
            print "Unable to project CoM in the center of the support polygone"
    v(sres.q())
    return sres
    
    
    
q = fullBody.getCurrentConfig()
q[:3] = [-2.69, 0.21, 0.649702]
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
for i in range(0, len(pb["phaseData"])):    
#for i in range(0, 5):    
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
    snew = gen_state(sprev, i , com, num_max_sample = 0, first = False, normal = n )
    all_states += [snew]   
    sprev = snew
    #~ com2 = coms[i+1].tolist()
    #~ snew2 = gen_state(snew, i , com2, num_max_sample = 20, normal = n, newContact = False )
    #~ all_states += [snew2]
    #~ sprev = snew2
    
#~ raise ValueError
    
#~ all_states = all_states[:-1]
configs = [ st.q() for st in all_states[:]]; i = 0

print "SID ", [s.sId for s in all_states]

beginId = 0
"""
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
    print "fail to project com"
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
        print "state ", j
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
"""
