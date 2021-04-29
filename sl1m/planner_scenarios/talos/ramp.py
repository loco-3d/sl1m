import numpy as np
import os
import sl1m.tools.plot_tools as plot

from talos_rbprm.talos import Robot as Talos
from sl1m.rbprm.surfaces_from_planning import getSurfacesFromGuideContinuous
from sl1m.stand_alone_scenarios.problem_definition_talos import Problem as TalosProblem
from sl1m.problem_definition import Problem
from sl1m.generic_solver import *

from time import perf_counter as clock

GAIT = [0, 1]
LIMB_NAMES = ["LF", "RF"]

USE_BIPED_PLANNER = True
USE_MIP = False

if __name__ == '__main__':
    t_init = clock()

    from sl1m.planner_scenarios.talos import lp_ramp_path as tp
    t_1 = clock()

    R, surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder, tp.ps, tp.afftool, tp.pathId, tp.v, 0.7, False)
    t_2 = clock()

    talos = Talos()
    lf_0 = np.array(tp.q_init[:3]) + np.array([0, 0.085, -0.98])
    rf_0 = np.array(tp.q_init[:3]) + np.array([0, -0.085, -0.98])
    initial_contacts = [lf_0, rf_0]
    t_3 = clock()

    if USE_BIPED_PLANNER:
        pb = TalosProblem()
        pb.generate_problem(R, surfaces, GAIT, initial_contacts)
        t_4 = clock()
        if USE_MIP:
            result = solve_MIP_biped(pb, surfaces)
        else:
            result = solve_L1_combinatorial_biped(pb, surfaces)
    else:
        talos.kinematic_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
            "/share/talos-rbprm/com_inequalities/feet_quasi_flat/talos_"
        talos.relative_feet_constraints_path = os.environ["INSTALL_HPP_DIR"] + \
            "/share/talos-rbprm/relative_effector_positions/talos_"
        pb = Problem(talos, suffix_com="_effector_frame_REDUCED.obj",
                     suffix_feet="_quasi_flat_REDUCED.obj", limb_names=LIMB_NAMES)
        pb.generate_problem(R, surfaces, GAIT, initial_contacts, tp.q_init[:3])
        t_4 = clock()
        
        if USE_MIP:
            result = solve_MIP(pb, surfaces)
        else:
            result = solve_L1_combinatorial(pb, surfaces)

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
