import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

import sl1m.tools.plot_tools as plot
from sl1m.tools.heightmap_tools import Heightmap

from solo_rbprm.robot import Robot

from hpp.corbaserver.affordance.affordance import AffordanceTool
from hpp.corbaserver.problem_solver import ProblemSolver
from hpp.gepetto import ViewerFactory

# --------------------------------- PROBLEM DEFINITION ---------------------------------------------------------------

ENV_URDF = "/opt/openrobots/share/hpp_environments/urdf/multicontact/bauzil_stairs.urdf"

HEIGHTMAP_FILENAME = "heightmap.pickle"
N_X = 100
N_Y = 50
X_BOUNDS = [-0.5, 2.5]
Y_BOUNDS = [0.5, 2.]

LIMBS = ['solo_RHleg_rom', 'solo_LHleg_rom', 'solo_LFleg_rom', 'solo_RFleg_rom']

# --------------------------------- METHODS ---------------------------------------------------------------


def init_affordance():
    """
    Initialize the affordance tool and return the solo abstract rbprm builder, the surface
    dictionary and all the affordance points
    """
    robot = Robot()
    robot.setJointBounds("root_joint", [-5., 5., -5., 5., 0.241, 1.5])
    robot.boundSO3([-3.14, 3.14, -0.01, 0.01, -0.01, 0.01])
    robot.setFilter(LIMBS)
    for limb in LIMBS:
        robot.setAffordanceFilter(limb, ['Support'])
    ps = ProblemSolver(robot)
    vf = ViewerFactory(ps)
    afftool = AffordanceTool()
    afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.00005])
    afftool.loadObstacleModel(ENV_URDF, "environment", vf)
    ps.selectPathValidation("RbprmPathValidation", 0.05)

    return afftool.getAffordancePoints('Support')


# --------------------------------- MAIN ---------------------------------------------------------------
if __name__ == "__main__":
    affordances = init_affordance()

    heightmap = Heightmap(N_X, N_Y, X_BOUNDS, Y_BOUNDS)
    heightmap.build(affordances)
    heightmap.save_pickle(HEIGHTMAP_FILENAME)

    ax_heightmap = plot.plot_heightmap(heightmap)
    plt.show()
