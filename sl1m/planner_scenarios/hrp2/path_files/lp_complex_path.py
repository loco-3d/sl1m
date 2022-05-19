from hpp.gepetto import PathPlayer
from hpp.corbaserver.affordance.affordance import AffordanceTool
from hpp.corbaserver.problem_solver import ProblemSolver
from hpp.corbaserver.rbprm.rbprmbuilder import Builder
from hpp.gepetto import Viewer
from hpp.corbaserver import Client
from hpp.corbaserver.robot import Robot as Parent

from hpp.corbaserver.rbprm.hrp2_abstract import Robot

packageName = "hpp_environments"
meshPackageName = "hpp_environments"

rbprmBuilder = Robot()
rbprmBuilder.setJointBounds("root_joint", [-20, 20, -10, 10, 0, 2.2])
rbprmBuilder.setFilter([Robot.rLegId, Robot.lLegId])

rbprmBuilder.setAffordanceFilter(
    Robot.rLegId,
    [
        "Support",
    ],
)
rbprmBuilder.setAffordanceFilter(Robot.lLegId, ["Support"])
rbprmBuilder.boundSO3([-0.0, 0, -1, 1, -1, 1])

ps = ProblemSolver(rbprmBuilder)
r = Viewer(ps)

q_init = rbprmBuilder.getCurrentConfig()
q_init[0:3] = [-1, -0.82, 0.5]
rbprmBuilder.setCurrentConfig(q_init)
r(q_init)

q_goal = q_init[::]
q_goal[3:7] = [0.0, 0.0, 0.14943813, 0.98877108]
q_goal[0:3] = [100.49, -0.65, 1.2]
r(q_goal)
ps.addPathOptimizer("RandomShortcut")
ps.setInitialConfig(q_init)
ps.addGoalConfig(q_goal)
r(q_goal)

afftool = AffordanceTool()
afftool.setAffordanceConfig("Support", [0.5, 0.03, 0.00005])
afftool.loadObstacleModel(
    packageName, "multicontact/bauzil_ramp", "planning", r
)
afftool.visualiseAffordances("Support", r, [0.25, 0.5, 0.5])

ps.client.problem.selectConfigurationShooter("RbprmShooter")
ps.client.problem.selectPathValidation("RbprmPathValidation", 0.05)

pp = PathPlayer(r)
