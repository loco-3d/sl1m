import numpy as np
import sl1m.tools.plot_tools as plot

from talos_rbprm.talos import Robot as Talos
from sl1m.rbprm.surfaces_from_planning import getSurfacesFromGuideContinuous
from sl1m.planner_scenarios.talos.problem_definition_talos import generate_problem
from sl1m.generic_solver import solve_L1_combinatorial_biped, solve_MIP_cost_biped

from time import perf_counter as clock

COSTS = {"step_size": None, "final_com": None, "effector_positions": None, "coms": None, "posture": True}

GAIT = [0, 1]

if __name__ == '__main__':
    t_init = clock()

    from sl1m.planner_scenarios.talos import lp_stairs_path as tp
    t_1 = clock()

    R, surfaces = getSurfacesFromGuideContinuous(tp.rbprmBuilder, tp.ps, tp.afftool, tp.pathId, tp.v, 0.7, True)
    t_2 = clock()

    talos = Talos()
    lf_0 = np.array(tp.q_init[:3]) + np.array([0, 0.085, -0.98])
    rf_0 = np.array(tp.q_init[:3]) + np.array([0, -0.085, -0.98])
    initial_contacts = [lf_0, rf_0]
    t_3 = clock()

    pb = generate_problem(talos, R, surfaces, GAIT, initial_contacts, eq_as_ineq=False)
    t_4 = clock()

    result = solve_MIP_cost_biped(pb, surfaces, costs=None)
    t_end = clock()

    print("Optimized number of steps:              ", pb["n_phases"])
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Computing the path takes                ", 1000. * (t_1 - t_init))
    print("Computing the surfaces takes            ", 1000. * (t_2 - t_1))
    print("Computing the initial contacts takes    ", 1000. * (t_3 - t_2))
    print("Generating the problem dictionary takes ", 1000. * (t_4 - t_3))
    print("Solving the problem takes               ", 1000. * (t_end - t_4))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(surfaces, GAIT)
    plot.plot_initial_contacts(initial_contacts, ax=ax)
    plot.plot_planner_result(result.coms, result.moving_foot_pos, result.all_feet_pos, ax, True)