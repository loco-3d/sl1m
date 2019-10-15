from hpp.corbaserver.rbprm.hrp2 import Robot
from hpp.gepetto import Viewer
from tools.display_tools import *
import time
print "Plan guide trajectory ..."
import lp_stair_bauzil_hrp2_path as tp
print "Done."
import time

pId = tp.ps.numberPaths() -1
fullBody = Robot ()

# Set the bounds for the root
fullBody.setJointBounds ("root_joint", [-10.135,2, -2, 2, 0, 2.2])
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


from sl1m.planner_scenarios.escaliers import *

pb, coms, footpos, allfeetpos, res = solve() 

def gen_state(s, pId, com , num_max_sample = 1, first = False ):
    #~ pId = 6
    phase =  pb["phaseData"][pId]
    moving = phase["moving"]
    movingID = fullBody.lLegId
    if moving == RF:
        movingID = fullBody.rLegId
    pos = allfeetpos[pId]; pos[2]+=0.01
    com[2] += 0.5
    print "com" , com
    print "pos" , pos.tolist()
    #~ q = fullBody.getCurrentConfig()
    #~ s, succ = state_alg.addNewContact(s, fullBody.rLegId, rfPos.tolist(), z.tolist())
    sres, succ = state_alg.addNewContact(s, movingID, pos.tolist(), z.tolist(), num_max_sample= num_max_sample)
    print "succes ?", succ
    succCom = False
    ite = 0
    #~ if first:
        #~ print "FIRST "
        #~ com[2] -= 0.25
    #~ while(not succCom and ite < 11):
    while(not succCom and ite < 17):
        succCom =  fullBody.projectStateToCOM(sres.sId ,com, num_max_sample)
        com[2] -= 0.05
        ite += 1
    print "COM?", succCom
    v(sres.q())
    return sres
    
    
q = fullBody.getCurrentConfig()
q[:2] = allfeetpos[1][:2].tolist()
v(q)
s = rbprmstate.State(fullBody, q = q, limbsIncontact = [])
all_states = [s]
sprev = s
idfirst = 2

for i in range(1, len(pb["phaseData"])):    
    com = (coms[i-i]  ).tolist()
    #~ if i > 2:
        #~ com = (coms[i-1] + (coms[i] - coms[i-1]) *0.8).tolist()
    snew = gen_state(sprev, i , com, num_max_sample = 200, first = i == idfirst )
    all_states += [snew]   
    sprev = snew
    com2 = coms[i].tolist()
    snew2 = gen_state(snew, i , com2, num_max_sample = 2 )
    all_states += [snew2]
    sprev = snew2
    
all_states = all_states[4:-1]
configs = [ st.q() for st in all_states[:]]; i = 0

print "SID ", [s.sId for s in all_states]

beginId = 4

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
            

#~ raise ValueError

def nil():
    global paths
    #~ for j in range(2,len(all_states)-2,2):    
    for j in range(2,len(all_states)-2,1):    
        #~ print "state ", j
        pathId = fullBody.limbRRT(all_states[j].sId, all_states[j+1].sId)
        paths += [(pathId, ps.pathLength(pathId))]
        play_int()
        
from hpp.gepetto import PathPlayer
pp = PathPlayer (v)
