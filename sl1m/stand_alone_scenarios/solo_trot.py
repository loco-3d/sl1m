
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.flat_ground import scene

import sl1m.tools.plot_tools as plot

USE_SL1M = False
USE_COM = True
gait = "trot"

GAITS = {}
GAITS["walk"] = [np.array([1, 0, 1, 1]), np.array([1, 1, 0, 1]), np.array([0, 1, 1, 1]), np.array([1, 1, 1, 0])]
GAITS["trot"] = [np.array([0, 1, 1, 0]), np.array([1, 0, 0, 1])]
GAITS["jumping_trot"] = [np.array([1, 0, 1, 0]), np.zeros(4), np.array([0, 1, 0, 1]), np.zeros(4)]
GAIT = GAITS[gait]

if gait == "jumping_trot":
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import jumping_trot_surfaces as surfaces
elif gait == "trot":
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import trot_surfaces as surfaces
else:
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import walk_surfaces as surfaces

STEP_LENGTH = [0.2, 0.0]


# paths = [os.environ["INSTALL_HPP_DIR"] + "/share/solo-rbprm/com_inequalities/feet_quasi_flat/",
#          os.environ["INSTALL_HPP_DIR"] + "/share/solo-rbprm/relative_effector_positions/"]
paths = [os.environ["SOLO3D_ENV_DIR"] + "/com_inequalities/",
         os.environ["INSTALL_HPP_DIR"] + "/solo-rbprm/relative_effector_positions/"]
suffix_com = "_effector_frame.obj"
others = ['HR_FOOT', 'HL_FOOT', 'FL_FOOT', 'FR_FOOT']
limbs = ['HRleg', 'HLleg', 'FLleg', 'FRleg']
offsets = {'FRleg':  [0.1946, -0.0875, -0.241], 'FLleg': [0.1946, 0.0875, -0.241],
           'HRleg': [-0.1946, -0.0875, -0.241], 'HLleg': [-0.1946, 0.0875, -0.241]}

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    q_init = [1., 0., 0.241]
    initial_contacts = []
    for limb in limbs:
        contact = np.array(q_init) + offsets[limb]
        contact[2] = 0.
        initial_contacts.append(contact)
    t_1 = clock()

    pb = Problem(limb_names=limbs, other_names=others, constraint_paths=paths, suffix_com=suffix_com)
    # pb = Problem(limb_names=limbs, other_names=others, constraint_paths=paths)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init[:3])
    t_2 = clock()

    coms = [np.array(q_init) + [0. + i * STEP_LENGTH[0]/2, 0., 1.] for i in range(len(surfaces))]
    print(coms)

    COSTS = {"step_size": [10.0, STEP_LENGTH]}   

    if USE_SL1M:
        result = solve_L1_combinatorial(pb, costs=COSTS, com=USE_COM)
    else:
        result = solve_MIP(pb, costs=COSTS, com=USE_COM)
    t_end = clock()

    print(result)
    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Initializing the problem takes          ", 1000. * (t_1 - t_init))
    print("Generating the problem dictionary takes ", 1000. * (t_2 - t_1))
    print("Solving the problem takes               ", 1000. * (t_end - t_2))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(scene)
    if(result.success):
        plot.plot_planner_result(result.all_feet_pos, coms=result.coms, step_size=STEP_LENGTH, ax=ax, show=True)
    else:
        plt.show()
