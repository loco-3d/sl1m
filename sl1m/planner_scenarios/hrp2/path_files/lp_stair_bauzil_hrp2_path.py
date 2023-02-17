from hpp.corbaserver.rbprm.rbprmbuilder import Builder
from hpp.gepetto import Viewer
from hpp.corbaserver import Client
from hpp.corbaserver.robot import Robot as Parent


from hpp.corbaserver.rbprm.hrp2_abstract import Robot

# ~ rootJointType = 'freeflyer'
packageName = "hpp-rbprm-corba"
meshPackageName = "hpp-rbprm-corba"
# ~ urdfName = 'hrp2_trunk_flexible'
# ~ urdfNameRoms =  ['hrp2_larm_rom','hrp2_rarm_rom','hrp2_lleg_rom','hrp2_rleg_rom']
# ~ urdfSuffix = ""
# ~ srdfSuffix = ""

rbprmBuilder = Robot()
rbprmBuilder.setJointBounds("root_joint", [-2, 2, -1, 1, 0, 2.2])
rbprmBuilder.setFilter([Robot.rLegId, Robot.lLegId])

rbprmBuilder.setAffordanceFilter(
    Robot.rLegId,
    [
        "Support",
    ],
)
rbprmBuilder.setAffordanceFilter(Robot.lLegId, ["Support"])
rbprmBuilder.boundSO3([-0.0, 0, -1, 1, -1, 1])

# ~ from hpp.corbaserver.rbprm. import ProblemSolver
from hpp.corbaserver.problem_solver import ProblemSolver

ps = ProblemSolver(rbprmBuilder)
r = Viewer(ps)


q_init = rbprmBuilder.getCurrentConfig()
q_init[0:3] = [-1, -0.82, 0.5]
rbprmBuilder.setCurrentConfig(q_init)
r(q_init)

q_goal = q_init[::]
q_goal[3:7] = [0.0, 0.0, 0.14943813, 0.98877108]
q_goal[0:3] = [1.49, -0.65, 1.2]
r(q_goal)
ps.addPathOptimizer("RandomShortcut")
ps.setInitialConfig(q_init)
ps.addGoalConfig(q_goal)

from hpp.corbaserver.affordance.affordance import AffordanceTool

afftool = AffordanceTool()
afftool.setAffordanceConfig("Support", [0.5, 0.03, 0.00005])
afftool.loadObstacleModel(packageName, "stair_bauzil", "planning", r)
afftool.visualiseAffordances("Support", r, [0.25, 0.5, 0.5])

ps.client.problem.selectConfigurationShooter("RbprmShooter")
ps.client.problem.selectPathValidation("RbprmPathValidation", 0.05)
t = ps.solve()


print(t)
# ~ if isinstance(t, list):
# ~ t = t[0]* 3600000 + t[1] * 60000 + t[2] * 1000 + t[3]
# ~ f = open('log.txt', 'a')
# ~ f.write("path computation " + str(t) + "\n")
# ~ f.close()


from hpp.gepetto import PathPlayer

pp = PathPlayer(r)
# ~ pp.fromFile("/home/stonneau/dev/hpp/src/hpp-rbprm-corba/script/paths/stair.path")
# ~

for i in range(1, 10):
    rbprmBuilder.client.problem.optimizePath(i)

# ~ pp (10)
# ~ pp (0)

# ~ pp (1)
# ~ pp.toFile(1, "/home/stonneau/dev/hpp/src/hpp-rbprm-corba/script/paths/stair.path")
# ~ rbprmBuilder.exportPath (r, ps.client.problem, 1, 0.01, "stair_bauzil_hrp2_path.txt")

# ~ cl = Client()
# ~ cl.problem.selectProblem("rbprm_path")
# ~ rbprmBuilder2 = Robot ("toto")
# ~ ps2 = ProblemSolver( rbprmBuilder2 )
# ~ cl.problem.selectProblem("default")
# ~ cl.problem.movePathToProblem(1,"rbprm_path",rbprmBuilder.getAllJointNames())
# ~ r2 = Viewer (ps2, viewerClient=r.client)
# ~ r.client.gui.setVisibility("toto", "OFF")
# ~ r.client.gui.setVisibility("hrp2_trunk_flexible", "OFF")
# ~ r2(q_far)
