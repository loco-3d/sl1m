
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from anymal_rbprm.anymal import Robot as Anymal

from sl1m.generic_solver import solve_MIP_gait
from sl1m.problem_definition_gait import Problem
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import quadruped_surfaces_gait as surfaces
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import scene

import sl1m.tools.plot_tools as plot

GAIT = [np.array([1, 1, 1, 0]), np.array([1, 1, 0, 1]), np.array([1, 0, 1, 1]), np.array([0, 1, 1, 1])]

USE_COM = False

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

    result = solve_MIP_gait(pb, com=USE_COM)
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
    if(result.success):
        plot.plot_planner_result(result.coms, result.all_feet_pos, ax, True)
    else:
        plt.show()
