
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock

from sl1m.generic_solver import solve_L1_combinatorial_biped, solve_MIP_biped
from sl1m.stand_alone_scenarios.problem_definition_hrp2 import Problem
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import surfaces

import sl1m.tools.plot_tools as plot

GAIT = [0, 1]

if __name__ == '__main__':
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    initial_contacts = [np.array([-0.0, 0.33499999999999996, 0.]),
                        np.array([-0.0, 0.145, 0.])]
    t_2 = clock()

    pb = Problem()
    pb.generate_problem(R, surfaces, GAIT, initial_contacts)
    t_3 = clock()

    result = solve_L1_combinatorial_biped(pb, surfaces)
    t_end = clock()

    print(result)

    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Computing the surfaces takes            ", 1000. * (t_1 - t_init))
    print("Computing the initial contacts takes    ", 1000. * (t_2 - t_1))
    print("Generating the problem dictionary takes ", 1000. * (t_3 - t_2))
    print("Solving the problem takes               ", 1000. * (t_end - t_3))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(surfaces, GAIT)
    plot.plot_initial_contacts(initial_contacts, ax=ax)
    if(result.success):
        plot.plot_planner_result(result.coms, result.moving_foot_pos, result.all_feet_pos, ax, True)
    else:
        plt.show()
