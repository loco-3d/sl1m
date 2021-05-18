
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from anymal_rbprm.anymal import Robot as Anymal

from sl1m.generic_solver import solve_L1_combinatorial_gait, solve_MIP_gait
from sl1m.problem_definition_gait import Problem
# ~ from sl1m.stand_alone_scenarios.surfaces.flat_ground import scene
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import surfaces as scene
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import walk_surfaces as surfaces

import sl1m.tools.plot_tools as plot

USE_SL1M = False

WALK = [np.array([1, 0, 1, 1]), np.array([1, 1, 0, 1]), np.array([0, 1, 1, 1]), np.array([1, 1, 1, 0])]
TROT = [np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1])]
GAIT = WALK


STEP_LENGTH = [0.3, 0.0]
COSTS = {"step_size": [1.0, STEP_LENGTH], "posture": [1.0]}

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    anymal = Anymal()
    q_init = anymal.referenceConfig.copy()
    q_init[:2] = [-0.5, 0.0]

    q_floor = q_init.copy()
    q_floor[2] = 0.
    initial_contacts = [(np.array(q_floor[:3]) + anymal.dict_limb_offset[foot])
                        for foot in anymal.limbs_names]

    t_2 = clock()

    anymal.kinematic_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
        "/share/anymal-rbprm/com_inequalities/feet_quasi_flat/anymal_"
    anymal.relative_feet_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
        "/share/anymal-rbprm/relative_effector_positions/anymal_"

    pb = Problem(anymal, suffix_com="_effector_frame_quasi_static_upscale.obj")
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init[:3])
    t_3 = clock()

    if USE_SL1M:
        result = solve_L1_combinatorial_gait(pb, surfaces, costs=COSTS)
    else:
        print ("MIP")
        result = solve_MIP_gait(pb, surfaces, costs=COSTS, com=False)
    t_end = clock()

    print(result)

    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Computing the surfaces takes            ", 1000. * (t_1 - t_init))
    print("Computing the initial contacts takes    ", 1000. * (t_2 - t_1))
    print("Generating the problem dictionary takes ", 1000. * (t_3 - t_2))
    print("Solving the problem takes               ", 1000. * (t_end - t_3))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(scene)
    plot.plot_initial_contacts(initial_contacts, ax=ax)
    if(result.success):
        plot.plot_planner_result(result.coms, result.all_feet_pos, ax, True)
    else:
        plt.show()
