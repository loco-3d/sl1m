import numpy as np
import os
import matplotlib.pyplot as plt
from time import perf_counter as clock

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.ramp_surfaces import surfaces_gait as surfaces
from sl1m.stand_alone_scenarios.surfaces.ramp_surfaces import scene

import sl1m.tools.plot_tools as plot

USE_COM = False
GAIT = [np.array([1, 0]), np.array([0, 1])]

paths = [os.environ["INSTALL_HPP_DIR"] + "/share/talos-rbprm/com_inequalities/feet_quasi_flat/talos_",
         os.environ["INSTALL_HPP_DIR"] + "/share/talos-rbprm/relative_effector_positions/talos_"]
limbs = ["LF", "RF"]
suffix_com = "_effector_frame_REDUCED.obj"
suffix_feet = "_quasi_flat_REDUCED.obj"

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    initial_contacts = [np.array([1.7,  0.285,  0.6]), np.array([1.7,  0.115,  0.6])]
    t_2 = clock()

    pb = Problem(limb_names=limbs, constraint_paths=paths, suffix_com=suffix_com, suffix_feet=suffix_feet)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts)
    t_3 = clock()

    result = solve_L1_combinatorial(pb, surfaces, com=USE_COM)
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
