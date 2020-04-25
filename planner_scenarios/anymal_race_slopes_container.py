from sl1m.planner_scenarios.state_methods import  *
from sl1m.planner_scenarios.anymal_constants import *

from sl1m.stand_alone_scenarios.anymal.flat_ground import solve

#load scene
from hpp.corbaserver.affordance.affordance import AffordanceTool
afftool = AffordanceTool ()
afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.05])
afftool.loadObstacleModel ("hpp_environments", "ori/race_slopes_container", "planning", v)
afftool.visualiseAffordances('Support', v, v.color.lightBrown)
v.addLandmark(v.sceneName,1)


#retrieve surfaces from scene for sl1m
from hpp.corbaserver.rbprm.tools import surfaces_from_path

# ~ initPos = [array([-0.8, 0.2, 0.13]),
 # ~ array([-0.8, -0.2, 0.13]),
 # ~ array([-1.6, 0.2, 0.13]),
 # ~ array([-1.6, -0.2, 0.13])]
# ~ endCom =  [10., 0., 0.4]


endCom =  [10., 0., 0.4]

from sl1m.stand_alone_scenarios.anymal.race_slopes_container import solve, overrideSurfaces


#starting position
q_init [0:3] = [0, 0, 0.58]
v(q_init)
initContacts, initPos, initCom = initConstraintsFrom_q_init(fullBody, q_init, limbNames)

initPos = [array([0.38437629, 0.18974121, 0.13850503]),
 array([ 0.38437629, -0.18974121,  0.1326605 ]),
 array([-0.38437629,  0.18974121,  0.11856727]),
 array([-0.38437629, -0.18974121,  0.05008679])]


#in case reduceSize changed
# ~ overrideSurfaces(surfaces_from_path.getAllSurfaces(afftool))

endCom = [0, 1.25, 0.4]

# ~ #compute contact sequence
pb, coms, footpos, allfeetpos, res = solve(initPos = initPos, endCom = endCom)

# ~ q_init[2] += footOffset
# ~ v(q_init)

# ~ #create init state
# ~ s = rbprmstate.State(fullBody, q = q_init, limbsIncontact = limbNames[:])
s = gen_init_state(initPos, fullBody, q_init, limbNames, footOffset, normal = z)
states = [s]
q_init = s.q()[:]
v(q_init)

# ~ #compute whole-body states from contact sequences
run(fullBody, states, allfeetpos, limbNames, coms, pb, footOffset)

# ~ #play motion

# ~ #export contact sequence
cs = exportCS(fullBody, q_init, states, pb, allfeetpos, limbNames, effectorNames, squeeze = True)

# ~ #save file
cs.saveAsBinary("anymal_race_slopes_container.cs")
play(states,v)

# ~ #display footsteps
# ~ displaySteppingStones(cs, v.client.gui, v.sceneName, fullBody)


