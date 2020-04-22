from sl1m.planner_scenarios.state_methods import  *
from sl1m.planner_scenarios.anymal_constants import *


from sl1m.stand_alone_scenarios.anymal.flat_ground import solve


initContacts, initPos, initCom = initConstraintsFrom_q_init(fullBody, q_init, limbNames)

# ~ pb, coms, footpos, allfeetpos, res = solve(initCom = None, initPos = None)

#compute contact sequence
pb, coms, footpos, allfeetpos, res = solve(initCom = initCom, initPos = initPos)


v(q_init)
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
