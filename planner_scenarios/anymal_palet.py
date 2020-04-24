from sl1m.planner_scenarios.state_methods import  *
from sl1m.planner_scenarios.anymal_constants import *

from sl1m.stand_alone_scenarios.anymal.flat_ground import solve

#load scene
from hpp.corbaserver.affordance.affordance import AffordanceTool
afftool = AffordanceTool ()
afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.05])
afftool.loadObstacleModel ("hpp_environments", "ori/modular_palet_flat", "planning", v,reduceSizes=[0.1,0,0])
afftool.visualiseAffordances('Support', v, v.color.lightBrown)
v.addLandmark(v.sceneName,1)


#retrieve surfaces from scene for sl1m
from hpp.corbaserver.rbprm.tools import surfaces_from_path

#starting position
q_init [0:3] = [-1.2,0,0.6]

initContacts, initPos, initCom = initConstraintsFrom_q_init(fullBody, q_init, limbNames)

initPos = [array([-0.8, 0.2, 0.13]),
 array([-0.8, -0.2, 0.13]),
 array([-1.6, 0.2, 0.13]),
 array([-1.6, -0.2, 0.13])]
endCom =  [10., 0., 0.4]

from sl1m.stand_alone_scenarios.anymal.palet import solve, overrideSurfaces

#in case reduceSize changed
overrideSurfaces(surfaces_from_path.getAllSurfaces(afftool))

#compute contact sequence
pb, coms, footpos, allfeetpos, res = solve(initPos = initPos, endCom = endCom)
# ~ pb, coms, footpos, allfeetpos, res = solve(initPos = initPos)


# ~ footOffset = 0.

q_init[2] += footOffset
v(q_init)

#create init state
s = rbprmstate.State(fullBody, q = q_init, limbsIncontact = limbNames[:])
states = [s]

#compute whole-body states from contact sequences
run(fullBody, states, allfeetpos, limbNames, coms, pb, footOffset)

#play motion

#export contact sequence
cs = exportCS(fullBody, q_init, states, pb, allfeetpos, limbNames, effectorNames)

#save file
cs.saveAsBinary("anymal_palet.cs")
play(states,v)

#display footsteps
displaySteppingStones(cs, v.client.gui, v.sceneName, fullBody)


