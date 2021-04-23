from time import perf_counter as clock
from sl1m.generic_solver import solve_L1_combinatorial_biped
from sl1m.planner_scenarios.talos.problem_definition_talos import generate_problem
from talos_rbprm.talos import Robot as Talos
import sl1m.tools.plot_tools as plot
import numpy as np

GAIT = [0, 1]

begin = np.array([[1.75, 1.65, 1.65, 1.75],
                  [0.3, 0.3, 0.1, 0.1],
                  [0.6, 0.6, 0.6, 0.6]])

platform = np.array([[2.5, 1.5, 1.5, 2.5],
                     [0.9, 0.9, -1.1, -1.1],
                     [0.6, 0.6, 0.6, 0.6]])

bridge = np.array([[-1.5, -1.5, 1.5, 1.5],
                   [-0.5, -0.8, -0.8, -0.5],
                   [0.6, 0.6, 0.6, 0.6]])

end = np.array([[-1.5, -2.0, -2.0, -1.5],
                [-0.4, -0.4, -1.1, -1.1],
                [0.6, 0.6, 0.6, 0.6]])


if __name__ == '__main__':
    t_init = clock()

    surfaces = []
    surfaces += [[begin]]+[[begin]]
    for i in range(20):
        surfaces += [[platform, bridge]]
    surfaces += [[end]]+[[end]]

    R = []
    for i in range(len(surfaces)):
        R.append(np.identity(3))
    t_1 = clock()

    initial_contacts = [np.array([1.7,  0.285,  0.6]), np.array([1.7,  0.115,  0.6])]
    t_2 = clock()

    talos = Talos()
    pb = generate_problem(talos, R, surfaces, GAIT, initial_contacts, eq_as_ineq=False)
    t_3 = clock()

    pb_data = solve_L1_combinatorial_biped(pb, surfaces)
    t_end = clock()

    print("Optimized number of steps:              ", pb["n_phases"])
    print("Total time is:                          ", 1000. * (t_end-t_init))
    print("Computing the surfaces takes            ", 1000. * (t_1 - t_init))
    print("Computing the initial contacts takes    ", 1000. * (t_2 - t_1))
    print("Generating the problem dictionary takes ", 1000. * (t_3 - t_2))
    print("Solving the problem takes               ", 1000. * (t_end - t_3))
    print("The LP and QP optimizations take        ", pb_data.time)

    ax = plot.draw_scene(surfaces, GAIT)
    plot.plot_initial_contacts(initial_contacts, ax=ax)
    plot.plot_planner_result(result.coms, result.moving_foot_pos, result.all_feet_pos, ax, True)
