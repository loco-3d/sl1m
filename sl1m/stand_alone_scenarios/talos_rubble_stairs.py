import os
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock

from sl1m.generic_solver import solve_L1_combinatorial_gait, solve_MIP_gait, solve_MIP, solve_L1_combinatorial
from sl1m.problem_definition_gait import Problem
from talos_rbprm.talos import Robot as Talos
from sl1m.stand_alone_scenarios.surfaces.rubble_stair_surfaces import gait_surfaces as surfaces
from sl1m.stand_alone_scenarios.surfaces.rubble_stair_surfaces import scene
import sl1m.tools.plot_tools as plot

GAIT = [0, 1]

USE_COM = True
GAIT = [np.array([1, 0]), np.array([0, 1])]

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    initial_contacts = [np.array([-2.7805096486250154, 0.335, 0.]),
                        np.array([-2.7805096486250154, 0.145, 0.])]
    t_2 = clock()

    talos = Talos()
    talos.kinematic_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
            "/share/talos-rbprm/com_inequalities/feet_quasi_flat/talos_"
    talos.relative_feet_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
            "/share/talos-rbprm/relative_effector_positions/talos_"
    pb = Problem(talos, suffix_com="_effector_frame_REDUCED.obj", suffix_feet="_quasi_flat_REDUCED.obj", limb_names= ["LF", "RF"])
    pb.generate_problem(R, surfaces, GAIT, initial_contacts)
    t_3 = clock()

    result = solve_L1_combinatorial_gait(pb, surfaces, com=USE_COM)
    # result = solve_MIP_gait(pb, com=USE_COM)
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

