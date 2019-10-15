from hpp.corbaserver.rbprm.talos import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time

import numpy as np
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

import sl1m.planner_scenarios.talos.complex1 as lp
tp = lp.tp
pb, coms, footpos, allfeetpos, res = lp.solve(tp.R, tp.seqs, tp.surfaces, plot = True) 


>>>>>>> Stashed changes
print "Done."
import time
#Robot.urdfSuffix+="_safeFeet"

DEFAULT_COM_HEIGHT = 0.86

pId = tp.ps.numberPaths() -1
fullBody = Robot ()

# Set the bounds for the root
fullBody.setJointBounds ("root_joint", [-10.135,2, -20, 20, 0, 2.8])
fullBody.client.robot.setDimensionExtraConfigSpace(6)
fullBody.client.robot.setExtraConfigSpaceBounds([0,0,0,0,0,0,0,0,0,0,0,0])
# add the 6 extraDof for velocity and acceleration (see *_path.py script)
ps = tp.ProblemSolver( fullBody )

#load the viewer
v = tp.Viewer (ps,viewerClient=tp.r.client, displayCoM = True)

# load a reference configuration
q_ref = fullBody.referenceConfig_elbowsUp[::] +[0]*6
q_init = q_ref[::] 
fullBody.setReferenceConfig(q_ref)
fullBody.setCurrentConfig (q_init)
fullBody.setPostureWeights(fullBody.postureWeightsRootRotationConstrained[::] + [0]*6)
#fullBody.usePosturalTaskContactCreation(True)

print "Generate limb DB ..."
tStart = time.time()
# generate databases : 

nbSamples = 10000
fullBody.addLimb(fullBody.rLegId,fullBody.rleg,fullBody.rfoot,fullBody.rLegOffset,fullBody.rLegNormal, fullBody.rLegx, fullBody.rLegy, nbSamples, "static", 0.01,kinematicConstraintsPath=fullBody.rLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.rLegId, "ReferenceConfiguration", True)
fullBody.addLimb(fullBody.lLegId,fullBody.lleg,fullBody.lfoot,fullBody.lLegOffset,fullBody.rLegNormal, fullBody.lLegx, fullBody.lLegy, nbSamples, "static", 0.01,kinematicConstraintsPath=fullBody.lLegKinematicConstraints,kinematicConstraintsMin = 0.3)
fullBody.runLimbSampleAnalysis(fullBody.lLegId, "ReferenceConfiguration", True)


tGenerate =  time.time() - tStart
print "Done."
print "Databases generated in : "+str(tGenerate)+" s"

configSize = fullBody.getConfigSize() -fullBody.client.robot.getDimensionExtraConfigSpace()
q_init[:2] = [0.1, -0.82]
v(q_init)

fullBody.resetJointsBounds()
from hpp.corbaserver.rbprm.rbprmstate import State,StateHelper


from sl1m.planner_scenarios.talos.complex1 import *
pb, coms, footpos, allfeetpos, res = solve(plot = False) 

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
q[:2] = [-2.62, 0.24] # init root position, should match init feet positions defined in complex1:l93 ("p0")
v(q)
s = State(fullBody, q = q, limbsIncontact = [fullBody.lLegId, fullBody.rLegId])

idfirst = 2
#coms[0] = array(s.getCenterOfMass())
coms[0][2] -= 0.02 


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
    com = (coms[i]).tolist()
    n = normal(pb["phaseData"][i])
    snew = gen_state(sprev, i , com, num_max_sample = 0, first = False, normal = n )
    all_states += [snew]   
    sprev = snew
    
configs = [ st.q() for st in all_states[:]]; i = 0
#displayContactSequence(v,configs)


print "SID ", [s.sId for s in all_states]

beginId = 0


"""
import tools.display_tools as disp
disp.createSphere("c",v,color = v.color.green)

disp.moveSphere('c',v,pos)
"""

