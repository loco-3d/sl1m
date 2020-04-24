from sl1m.planner_scenarios.state_methods import  *
from sl1m.planner_scenarios.anymal_constants import *

from sl1m.stand_alone_scenarios.anymal.flat_ground import solve

#load scene
from hpp.corbaserver.affordance.affordance import AffordanceTool
afftool = AffordanceTool ()
afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.00005])
afftool.loadObstacleModel ("hpp_environments", "multicontact/ori_plinth", "planning", v,reduceSizes=[0.2,0,0])
afftool.visualiseAffordances('Support', v, v.color.lightBrown)
v.addLandmark(v.sceneName,1)


#retrieve surfaces from scene for sl1m
from hpp.corbaserver.rbprm.tools import surfaces_from_path

#starting position
q_init [0:3] = [-1.7,0,0.465]
initContacts, initPos, initCom = initConstraintsFrom_q_init(fullBody, q_init, limbNames)
initPos = [array([-1.3, 0.2, -0.0]), array([-1.3, -0.2, -0.0]), array([-2.1, 0.2, -0.0]), array([-2.1, -0.2, -0.0])]
v(q_init)



from sl1m.stand_alone_scenarios.anymal.plinth import solve, overrideSurfaces

#in case reduceSize changed
overrideSurfaces(surfaces_from_path.getAllSurfaces(afftool))

#compute contact sequence
print ("initPos ", initPos)
pb, coms, footpos, allfeetpos, res = solve(initPos=initPos)
# ~ pb, coms, footpos, allfeetpos, res = solve()

#create init state
s = rbprmstate.State(fullBody, q = q_init, limbsIncontact = limbNames[:])
states = [s]

#compute whole-body states from contact sequences
run(fullBody, states, allfeetpos, limbNames, coms, pb)

#play motion
play(states,v)

#export contact sequence
cs = exportCS(fullBody, q_init, states, pb, allfeetpos, limbNames, effectorNames)

#save file
cs.saveAsBinary("anymal_flatGround.cs")

#display footsteps
displaySteppingStones(cs, v.client.gui, v.sceneName, fullBody)


