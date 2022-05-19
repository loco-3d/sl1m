from hpp.gepetto import PathPlayer, Viewer, ViewerFactory
from hpp.corbaserver import Client
from hpp.corbaserver.affordance.affordance import AffordanceTool
from hpp.corbaserver.problem_solver import ProblemSolver
from talos_rbprm.talos_abstract import Robot as TalosAbstract
import time

TalosAbstract.urdfName += "_large"


def compute_path():
    talos_abstract = TalosAbstract()
    talos_abstract.setJointBounds(
        "root_joint", [-3.2, 1.8, 0.19, 0.21, 0.95, 1.7]
    )
    # As this scenario only consider walking, we fix the DOF of the torso :
    talos_abstract.setJointBounds("torso_1_joint", [0, 0])
    talos_abstract.setJointBounds("torso_2_joint", [0.0, 0.0])
    vMax = 1.0  # linear velocity bound for the root
    aMax = 2.0  # linear acceleration bound for the root
    extraDof = 6
    mu = 0.5  # coefficient of friction
    talos_abstract.setFilter([TalosAbstract.rLegId, TalosAbstract.lLegId])

    talos_abstract.setAffordanceFilter(
        TalosAbstract.rLegId,
        [
            "Support",
        ],
    )
    talos_abstract.setAffordanceFilter(TalosAbstract.lLegId, ["Support"])
    talos_abstract.boundSO3([-4.0, 4.0, -0.1, 0.1, -0.1, 0.1])
    # Add 6 extraDOF to the problem, used to store the linear velocity and acceleration of the root
    talos_abstract.client.robot.setDimensionExtraConfigSpace(extraDof)
    # We set the bounds of this extraDof with velocity and acceleration bounds (expect on z axis)
    extraDofBounds = [
        -vMax,
        vMax,
        -vMax,
        vMax,
        -10.0,
        10.0,
        -aMax,
        aMax,
        -aMax,
        aMax,
        -10.0,
        10.0,
    ]
    talos_abstract.client.robot.setExtraConfigSpaceBounds(extraDofBounds)
    (
        talos_abstract.getConfigSize()
        - talos_abstract.client.robot.getDimensionExtraConfigSpace()
    )

    ps = ProblemSolver(talos_abstract)
    vf = ViewerFactory(ps)

    afftool = AffordanceTool()
    afftool.setAffordanceConfig("Support", [0.5, 0.03, 0.00005])
    afftool.loadObstacleModel(
        "package://hpp_environments/urdf/multicontact/bauzil_ramp_simplified.urdf",
        "planning",
        vf,
        reduceSizes=[0.07, 0.0, 0.0],
    )
    v = vf.createViewer(displayArrows=True)
    afftool.visualiseAffordances("Support", v, [0.25, 0.5, 0.5])
    v.addLandmark(v.sceneName, 1)

    ps.setParameter("Kinodynamic/velocityBound", vMax)
    ps.setParameter("Kinodynamic/accelerationBound", aMax)
    # force the orientation of the trunk to match the direction of the motion
    ps.setParameter("Kinodynamic/forceYawOrientation", True)
    ps.setParameter("Kinodynamic/synchronizeVerticalAxis", True)
    ps.setParameter("Kinodynamic/verticalAccelerationBound", 10.0)
    ps.setParameter("DynamicPlanner/sizeFootX", 0.2)
    ps.setParameter("DynamicPlanner/sizeFootY", 0.12)
    ps.setParameter("DynamicPlanner/friction", mu)
    # sample only configuration with null velocity and acceleration :
    ps.setParameter("ConfigurationShooter/sampleExtraDOF", False)
    ps.setParameter("PathOptimization/RandomShortcut/NumberOfLoops", 100)
    # Choosing RBPRM shooter and path validation methods.
    ps.selectConfigurationShooter("RbprmShooter")
    ps.addPathOptimizer("RandomShortcutDynamic")
    ps.selectPathValidation("RbprmPathValidation", 0.05)
    # Choosing kinodynamic methods :
    ps.selectSteeringMethod("RBPRMKinodynamic")
    ps.selectDistance("Kinodynamic")
    ps.selectPathPlanner("DynamicPlanner")

    ### BEGIN climb the stairs #####
    ps.setParameter("Kinodynamic/velocityBound", 0.3)
    ps.setParameter("Kinodynamic/accelerationBound", 0.1)
    q_init = talos_abstract.getCurrentConfig()
    q_init[0:3] = [-0.3, 0.2, 0.98]
    v(q_init)  # between rubble and stairs
    q_goal = q_init[::]
    q_goal[0:3] = [1.7, 0.2, 1.58]
    v(q_goal)  # top of stairs
    q_goal[-6:-3] = [0, 0, 0]
    ps.setInitialConfig(q_init)
    q_init_0 = q_init[::]
    ps.addGoalConfig(q_goal)
    v(q_goal)

    ps.solve()
    print("done planning, optimize path ...")
    for _ in range(5):
        ps.optimizePath(ps.numberPaths() - 1)

    pId_stairs = ps.numberPaths() - 1
    ### END climb the stairs #####
    talos_abstract.setJointBounds(
        "root_joint", [-3.2, 2.5, -0.8, 0.3, 1.4, 2.0]
    )
    ps.resetGoalConfigs()
    ### BEGIN turn around on the platform #####
    ps.setParameter("Kinodynamic/velocityBound", 0.2)
    ps.setParameter("Kinodynamic/accelerationBound", 0.07)
    q_init = talos_abstract.getCurrentConfig()
    q_init = q_goal[::]
    v(q_init)  # top of stairs
    q_init[-6:-3] = [0.2, 0, 0]
    q_goal = q_init[::]
    q_goal[0:3] = [1.7, -0.6, 1.58]
    v(q_goal)  # before bridge
    q_goal[3:7] = [0, 0, 1, 0]
    q_goal[-6:-3] = [-0.2, 0, 0]
    ps.setInitialConfig(q_init)
    ps.addGoalConfig(q_goal)
    v(q_goal)

    ps.solve()
    print("done planning, optimize path ...")
    for _ in range(5):
        ps.optimizePath(ps.numberPaths() - 1)

    pId_platform = ps.numberPaths() - 1
    ### END turn around on the platform #####
    ps.resetGoalConfigs()
    ### BEGIN bridge cross #####
    ps.setParameter("Kinodynamic/velocityBound", 0.3)
    ps.setParameter("Kinodynamic/accelerationBound", 0.2)
    q_init = talos_abstract.getCurrentConfig()
    q_init = q_goal[::]
    v(q_init)  # top of stairs
    q_goal[0:3] = [-1.7, -0.6, 1.58]
    v(q_goal)  # after bridge
    q_goal[3:7] = [0, 0, 1, 0]
    q_goal[-6:-3] = [0.0, 0, 0]
    ps.setInitialConfig(q_init)
    ps.addGoalConfig(q_goal)
    v(q_goal)

    ps.solve()
    print("done planning, optimize path ...")
    for _ in range(5):
        ps.optimizePath(ps.numberPaths() - 1)

    pId_bridge = ps.numberPaths() - 1
    pathId = pId_stairs
    ps.concatenatePath(pathId, pId_platform)
    ps.concatenatePath(pathId, pId_bridge)

    print("done optimizing.")
    pp = PathPlayer(v)
    pp.dt = 0.1
    pp.displayVelocityPath(pathId)
    v.client.gui.setVisibility(
        "path_" + str(pathId) + "_root", "ALWAYS_ON_TOP"
    )
    pp.dt = 0.01

    q_far = q_goal[::]
    q_far[2] = -5
    v(q_far)
    q_init = q_init_0[::]

    return talos_abstract, ps, afftool, pathId, v, q_init
