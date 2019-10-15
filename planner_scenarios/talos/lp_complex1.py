from hpp.corbaserver.rbprm.talos import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time
from pinocchio import Quaternion
import numpy as np
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm

import sl1m.planner_scenarios.talos.complex1 as lp
tp = lp.tp
pb, coms, footpos, allfeetpos, res = lp.solve() 


print "Done."
import time
#Robot.urdfSuffix+="_safeFeet"

DEFAULT_COM_HEIGHT = 0.86
USE_ORIENTATION = True
pId = tp.ps.numberPaths() -1
fullBody = Robot ()

# Set the bounds for the root
fullBody.setJointBounds ("root_joint", [-20,20 ,-20, 20, 0, 2.8])
fullBody.client.robot.setDimensionExtraConfigSpace(6)
fullBody.client.robot.setExtraConfigSpaceBounds(tp.extraDofBounds)
# add the 6 extraDof for velocity and acceleration (see *_path.py script)
ps = tp.ProblemSolver( fullBody )

#load the viewer
v = tp.Viewer (ps,viewerClient=tp.v.client, displayCoM = True)

# load a reference configuration
q_ref = fullBody.referenceConfig_elbowsUp[::] +[0]*6
fullBody.setCurrentConfig (q_ref)
fullBody.setReferenceConfig(q_ref)
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
q_init = q_ref[::] 
q_init[0:7] = tp.q_init[0:7]
q_init[2] += q_ref[2] - 0.98 # 0.98 is the q_ref[2] for the rom device
v(q_init)

fullBody.resetJointsBounds()
from hpp.corbaserver.rbprm.rbprmstate import State,StateHelper

def quatConfigFromMatrix(m):
  quat = Quaternion(m)
  return [quat.x,quat.y,quat.z,quat.w]

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
      

def gen_state(s, pId, num_max_sample = 0, first = False, normal = lp.Z_AXIS, newContact = True  ,projectCOM = True):
    #~ pId = 6
    phase =  pb["phaseData"][pId]
    moving = phase["moving"]
    movingID = fullBody.lLegId
    if moving == lp.RF:
        movingID = fullBody.rLegId
    print "# gen state for phase Id = ",pId
    if USE_ORIENTATION:
      if pId < len(pb["phaseData"])-1:
        if phase["moving"] == lp.RF: 
          rot = quatConfigFromMatrix(pb["phaseData"][pId+1]["rootOrientation"]) # rotation of the root, from the guide
        else :
          rot = quatConfigFromMatrix(pb["phaseData"][pId]["rootOrientation"]) # rotation of the root, from the guide
        #quat0 = Quaternion(pb["phaseData"][pId]["rootOrientation"])
        #quat1 = Quaternion(pb["phaseData"][pId+1]["rootOrientation"])
        #rot = quatConfigFromMatrix((quat0.slerp(0.5,quat1)).matrix())
      else:
        rot = quatConfigFromMatrix(phase["rootOrientation"]) # rotation of the root, from the guide
    else :    
      rot = [0,0,0,1]
    #rot = quatConfigFromMatrix(phase["rootOrientation"]) # rotation of the root, from the guide
    #print "current config q=",s.q()
    #print "move limb ",movingID
    pos = allfeetpos[pId];
    print "Try to add contact for "+movingID+" pos = "+str(pos.tolist()+rot)
    disp.moveSphere('c',v,pos.tolist()+rot)
    if newContact:
        sres, succ = StateHelper.addNewContact(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= num_max_sample,rotation = rot)
        #sres, succ = StateHelper.removeContact(s,movingID)
        #assert succ
        #succ = sres.projectToRoot(s.q()[0:3]+rot)
        #sres, succ = StateHelper.addNewContact(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= num_max_sample)
    else:
        sres, succ = StateHelper.cloneState(s)
    if not succ:
      print "Cannot project config q = ",sres.q()
      print "To new contact position for "+movingID+" = "+str(pos.tolist()+rot)+" ; n = "+str(normal.tolist())
      raise RuntimeError("Cannot project feet to new contact position") # try something else ?? 
    if projectCOM :
        #print "config before projecting to com q1=",sres.q()
        successCOM = projectCoMInSupportPolygon(sres)
        if not successCOM:
            # is it really an issue ? 
            print "Unable to project CoM in the center of the support polygone"
        
    v(sres.q())
    return sres
    
import tools.display_tools as disp
disp.createSphere("c",v,color = v.color.green)
v.addLandmark("c",0.1)

    
    
q = fullBody.getCurrentConfig()
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
for i in range(2, len(pb["phaseData"])):    
    n = normal(pb["phaseData"][i])
    snew = gen_state(sprev, i , num_max_sample = 10000, first = False, normal = n )
    all_states += [snew]   
    sprev = snew
    
configs = [ st.q() for st in all_states[:]]; i = 0
#displayContactSequence(v,configs)


print "SID ", [s.sId for s in all_states]

beginId = 0


"""
import tools.display_tools as disp
disp.createSphere("c",v,color = v.color.green)
v.addLandmark("c",0.1)

disp.moveSphere('c',v,pos)
"""

