
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from solo_rbprm.solo import Robot as Solo

from sl1m.solver import Solvers
from sl1m.generic_solver import solve_L1_combinatorial_gait, solve_MIP_gait
from sl1m.problem_definition_gait import Problem
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import solo_surfaces_gait as surfaces
# from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import solo_surfaces_gait_small as surfaces
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import solo_scene as scene

import sl1m.tools.plot_tools as plot

GAIT = [np.array([1, 0, 1, 1]), np.array([1, 1, 0, 1]), np.array([0, 1, 1, 1]), np.array([1, 1, 1, 0])]
COSTS = {"posture": [1.0]}

USE_SL1M = False
USE_COM = True
TEST_CBC = True

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    solo = Solo()
    q_init = solo.referenceConfig.copy()
    q_init[:2] = [-0.25, 0.1]

    initial_contacts = [np.array(q_init[:3]) +
                        solo.dict_ref_effector_from_root[limb_name] +
                        solo.dict_offset[solo.dict_limb_joint[limb_name]].translation
                        for limb_name in solo.limbs_names]

    t_2 = clock()

    solo.kinematic_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
        "/share/solo-rbprm/com_inequalities/feet_quasi_flat/"
    solo.relative_feet_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
        "/share/solo-rbprm/relative_effector_positions/"

    pb = Problem(solo)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init[:3])
    t_3 = clock()

    if USE_SL1M:
        result = solve_L1_combinatorial_gait(pb, surfaces, costs=COSTS, com=USE_COM)
    else:
        if TEST_CBC:
            result = solve_MIP_gait(pb, costs=COSTS, com=USE_COM)
        else:
            result = solve_MIP_gait(pb, costs=COSTS, com=USE_COM)
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
        plt.show(block=False)
        
    if TEST_CBC:
        print ("CBC results")
        t_3 = clock()
        result = solve_MIP_gait(pb, costs=COSTS, com=USE_COM, solver=Solvers.CVXPY)
        t_end = clock()
        
        print(result)

        print("Optimized number of steps:              ", pb.n_phases)
        print("Solving the problem takes               ", 1000. * (t_end - t_3))
        print("The LP and QP optimizations take        ", result.time)
        
        ax = plot.draw_scene(scene)
        if(result.success):
            plot.plot_planner_result(result.coms, result.all_feet_pos, ax, True)
        else:
            plt.show(block=False)
