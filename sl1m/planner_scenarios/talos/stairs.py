import numpy as np
import os
import sl1m.tools.plot_tools as plot
import matplotlib.pyplot as plt

from sl1m.rbprm.surfaces_from_planning import getSurfacesFromGuideContinuous
from sl1m.stand_alone_scenarios.problem_definition_talos import Problem as TalosProblem
from sl1m.problem_definition import Problem
from sl1m.generic_solver import *

from time import perf_counter as clock

GAIT = [np.array([1, 0]), np.array([0, 1])]

USE_BIPED_PLANNER = False
USE_MIP = False
USE_COM = True

paths = [os.environ["INSTALL_HPP_DIR"] + "/share/talos-rbprm/com_inequalities/feet_quasi_flat/talos_",
         os.environ["INSTALL_HPP_DIR"] + "/share/talos-rbprm/relative_effector_positions/talos_"]
limbs = ["LF", "RF"]
suffix_com = "_effector_frame_REDUCED.obj"
suffix_feet = "_quasi_flat_REDUCED.obj"

if __name__ == '__main__':
    t_init = clock()

    from sl1m.planner_scenarios.talos import lp_stairs_path as tp
    t_1 = clock()

    R, surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder, tp.ps, tp.afftool, tp.pathId, tp.v, 0.7, False)
    t_2 = clock()

    p0 = [np.array(tp.q_init[:3]) + [0, 0.085, -0.98], np.array(tp.q_init[:3]) + [0, -0.085, -0.98]]
    t_3 = clock()

    if USE_BIPED_PLANNER:
        pb = TalosProblem()
        pb.generate_problem(R, surfaces, [0, 1], p0)
        t_4 = clock()
        if USE_MIP:
            result = solve_MIP_biped(pb)
        else:
            result = solve_L1_combinatorial_biped(pb, surfaces)
    else:
        surfaces_gait = [[surface] for surface in surfaces]

        pb = Problem(limb_names=limbs, constraint_paths=paths, suffix_com=suffix_com, suffix_feet=suffix_feet)
        pb.generate_problem(R, surfaces_gait, GAIT, p0, tp.q_init[:3])
        t_4 = clock()

        if USE_MIP:
            result = solve_MIP(pb, com=USE_COM)
        else:
            result = solve_L1_combinatorial(pb, surfaces_gait, com=USE_COM)

    t_end = clock()

    print(result)

    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Computing the path takes                ", 1000. * (t_1 - t_init))
    print("Computing the surfaces takes            ", 1000. * (t_2 - t_1))
    print("Computing the initial contacts takes    ", 1000. * (t_3 - t_2))
    print("Generating the problem dictionary takes ", 1000. * (t_4 - t_3))
    print("Solving the problem takes               ", 1000. * (t_end - t_4))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(surfaces)
    if(result.success):
        plot.plot_planner_result(result.coms, result.all_feet_pos, ax=ax, show=True)
    else:
        plt.show()
