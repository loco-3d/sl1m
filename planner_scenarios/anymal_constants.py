from anymal_rbprm.anymal import Robot
from hpp.gepetto import Viewer
from hpp.corbaserver import Client
from tools.display_tools import *
import time

from hpp.corbaserver.rbprm import  rbprmstate
from hpp.corbaserver.rbprm import  state_alg

from hpp.corbaserver.rbprm.rbprmstate import StateHelper

from numpy import array, zeros, ones
from hpp.corbaserver.problem_solver import ProblemSolver
from sl1m.stand_alone_scenarios.constraints_anymal import limbNames, effectorNames


packageName = 'hpp-rbprm-corba'
meshPackageName = 'hpp-rbprm-corba'

fullBody = Robot ()
fullBody.setJointBounds ("root_joint", [-20,20, -10, 10, 0, 2.2])
fullBody.client.robot.setDimensionExtraConfigSpace(0)


from hpp.corbaserver.problem_solver import ProblemSolver
ps = ProblemSolver( fullBody )
r = Viewer (ps)
v = Viewer (ps,viewerClient=r.client, displayCoM = True)

# ~ fullBody.setConstrainedJointsBounds()
fullBody.setJointBounds("LF_KFE",[-1.4,0.])
fullBody.setJointBounds("RF_KFE",[-1.4,0.])
fullBody.setJointBounds("LH_KFE",[0.,1.4])
fullBody.setJointBounds("RH_KFE",[0.,1.4])
fullBody.setJointBounds ("root_joint", [-20,20, -20, 20, -20, 20])
dict_heuristic = {fullBody.rLegId:"static", fullBody.lLegId:"static", fullBody.rArmId:"fixedStep04", fullBody.lArmId:"fixedStep04"}
fullBody.loadAllLimbs(dict_heuristic,"ReferenceConfiguration",nbSamples=12)

q_ref = fullBody.referenceConfig[::]
q_init = q_ref[::] 
fullBody.setReferenceConfig(q_ref)
fullBody.setCurrentConfig (q_init)
fullBody.setPostureWeights(fullBody.postureWeights[::])
footOffset = 0.03
