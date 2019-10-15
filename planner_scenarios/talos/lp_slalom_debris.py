from hpp.corbaserver.rbprm.talos import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time
from pinocchio import Quaternion
import numpy as np
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate,isclose
from numpy.linalg import norm
import random
from math import isnan
import sl1m.planner_scenarios.talos.slalom_debris as lp
tp = lp.tp
pb, coms, footpos, allfeetpos, res = lp.solve() 


print "Done."
import time
Robot.urdfSuffix+="_safeFeet"

DEFAULT_COM_HEIGHT = 0.86
Z_AXIS = lp.Z_AXIS
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
fullBody.setPostureWeights(fullBody.postureWeights[::] + [0]*6)
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
  return quatToConfig(quat)

def quatToConfig(quat):
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
    print "project state to com : ",desiredCOM
    q_save = s.q()[::]
    while not success and maxIt > 0:
      success = s.fullBody.projectStateToCOM(s.sId ,desiredCOM, maxNumSample = 0)
      maxIt -= 1 
      desiredCOM[2] -= 0.005
    print "success = ",success
    print "result = ",s.q()
    if success and isnan(s.q()[0]): # FIXME why does it happen ?
      success = False
      s.setQ(q_save)
    return success
 
def tryCreateContactAround(s, eff_id, pos, normal, num_max_sample= 0,rotation = Quaternion.Identity(),num_try = 1000):
    sres, succ = StateHelper.addNewContact(s, eff_id, pos, normal, num_max_sample= num_max_sample,rotation = quatToConfig(rotation))     
    i_try = 0
    x_bounds = [-0.1,0.1]
    y_bounds = [-0.05,0.05]
    print "try create contact around, first try success : ",succ
    print "result = ",sres.q()
    while not succ and i_try < num_try:
      print "try create contact around, try : ",i_try
      x = random.uniform(x_bounds[0],x_bounds[1])
      y = random.uniform(y_bounds[0],y_bounds[1])
      offset = np.matrix([x,y,0]).T
      offset = rotation.matrix() * offset
      new_pos = pos[::]
      for i in range(3):
        new_pos[i] += offset[i,0]
      sres, succ = StateHelper.addNewContact(s, eff_id, new_pos, normal, num_max_sample= num_max_sample,rotation = quatToConfig(rotation))     
      i_try += 1

    return sres,succ

def gen_state(s, pId, num_max_sample = 0, first = False, normal = lp.Z_AXIS, newContact = True  ,projectCOM = True):
    #~ pId = 6
    phase =  pb["phaseData"][pId]
    moving = phase["moving"]
    movingID = fullBody.lLegId
    if moving == lp.RF:
        movingID = fullBody.rLegId
    print "# gen state for phase Id = ",pId
    if False and pId < len(pb["phaseData"])-1:
      quat0 = Quaternion(pb["phaseData"][pId]["rootOrientation"])
      quat1 = Quaternion(pb["phaseData"][pId+1]["rootOrientation"])
      qrot = (quat0.slerp(0.5,quat1)).matrix()
    else:
      qrot = Quaternion(phase["rootOrientation"]) # rotation of the root, from the guide
    
    #q_n = Quaternion().FromTwoVectors(np.matrix(Z_AXIS).T,np.matrix(normal).T)
    #rot = quatToConfig(qrot * q_n)
    if not isclose(normal,Z_AXIS).all():
      qrot = Quaternion().FromTwoVectors(np.matrix(Z_AXIS).T,np.matrix(normal).T)
      # ignore guide orientation when normal is not z ...
    #rot = quatToConfig(qrot)
    pos = allfeetpos[pId];
    pos[2] += 0.002
    pose = pos.tolist()+quatToConfig(qrot)
    print "Try to add contact for "+movingID+" pos = "+str(pose)
    disp.moveSphere('c',v,pose)
    if newContact:
        sres, succ = tryCreateContactAround(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= num_max_sample,rotation = qrot)
        #sres, succ = StateHelper.removeContact(s,movingID)
        #assert succ
        #succ = sres.projectToRoot(s.q()[0:3]+rot)
        #sres, succ = StateHelper.addNewContact(s, movingID, pos.tolist(), normal.tolist(), num_max_sample= num_max_sample)
    else:
        sres, succ = StateHelper.cloneState(s)
    if not succ:
      print "Cannot project config q = ",sres.q()
      print "To new contact position for "+movingID+" = "+str(pose)+" ; n = "+str(normal.tolist())
      raise RuntimeError("Cannot project feet to new contact position") # try something else ?? 
    if projectCOM :
        #sfeet, _ = StateHelper.cloneState(sres)
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
    i = 0 
    n = np.zeros(3)
    while i < s.shape[1]-3 and norm(n) < 1e-4 :
      n = cross(s[:,i+1] - s[:,i+0], s[:,i+2] - s[:,i+0])
      i+=1
    if norm(n) < 1e-4:
      return Z_AXIS
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
    snew = gen_state(sprev, i , num_max_sample = 100, first = False, normal = n )
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

disp.moveSphere('c',v,p)
"""

