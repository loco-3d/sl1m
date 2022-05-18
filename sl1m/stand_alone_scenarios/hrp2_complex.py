import numpy as np
import os
import matplotlib.pyplot as plt
from time import perf_counter as clock

import sl1m.stand_alone_scenarios
from sl1m.generic_solver import solve_L1_combinatorial
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.complex_surfaces import surfaces_gait as surfaces
from sl1m.stand_alone_scenarios.surfaces.complex_surfaces import scene

import sl1m.tools.plot_tools as plot

GAIT = [np.array([1, 0]), np.array([0, 1])]
USE_COM = False

if __name__ == '__main__':
    paths = [os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/" ,
             os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/"]
    limbs = ["LF", "RF"]    

    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    initial_contacts = [np.array([-2.7805096486250154, 0.33499999999999996, 0.]),
                        np.array([-2.7805096486250154, 0.145, 0.])]
    t_2 = clock()

    pb = Problem(limb_names=limbs, constraint_paths=paths)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts)
    t_3 = clock()

    # result = solve_MIP(pb, com=USE_COM)
    result = solve_L1_combinatorial(pb, com=USE_COM)
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
        plot.plot_planner_result(result.all_feet_pos, coms=result.coms, ax=ax, show=True)
    else:
        plt.show()