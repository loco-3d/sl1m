import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
from time import perf_counter as clock
import talos_rbprm

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.ramp_surfaces import surfaces_gait as surfaces
from sl1m.stand_alone_scenarios.surfaces.ramp_surfaces import scene

import sl1m.tools.plot_tools as plot

USE_SL1M = True
USE_COM = False

GAIT = [np.array([1, 0]), np.array([0, 1])]

talos_rbprm_path = (
    Path(talos_rbprm.__file__).resolve().parent.parent.parent.parent.parent
    / "share"
    / "talos-rbprm"
)
paths = [
    str(talos_rbprm_path / "com_inequalities" / "feet_quasi_flat") + os.sep,
    str(talos_rbprm_path / "relative_effector_positions") + os.sep,
]
limbs = ["LF", "RF"]
suffix_com = "_effector_frame_REDUCED.obj"
suffix_feet = "_quasi_flat_REDUCED.obj"

if __name__ == "__main__":
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    t_1 = clock()

    initial_contacts = [np.array([1.7, 0.285, 0.6]), np.array([1.7, 0.115, 0.6])]
    t_2 = clock()

    pb = Problem(
        limb_names=limbs,
        constraint_paths=paths,
        suffix_com=suffix_com,
        suffix_feet=suffix_feet,
    )
    pb.generate_problem(R, surfaces, GAIT, initial_contacts)
    t_3 = clock()

    COSTS = {
        "step_size": [10.0, [0.1, 0]],
        "posture": [1.0],
        # "final_com": [10.0, [2.5, 0., 0.241]]
    }

    if USE_SL1M:
        result = solve_L1_combinatorial(pb, com=USE_COM, costs=COSTS)
    else:
        result = solve_MIP(pb, com=USE_COM, costs=COSTS)
    t_end = clock()

    print(result)
    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000.0 * (t_end - t_init))
    print("Computing the surfaces takes            ", 1000.0 * (t_1 - t_init))
    print("Computing the initial contacts takes    ", 1000.0 * (t_2 - t_1))
    print("Generating the problem dictionary takes ", 1000.0 * (t_3 - t_2))
    print("Solving the problem takes               ", 1000.0 * (t_end - t_3))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(scene)
    plot.plot_initial_contacts(initial_contacts, ax=ax)
    if result.success:
        plot.plot_planner_result(
            result.all_feet_pos, coms=result.coms, ax=ax, show=True
        )
    else:
        plt.show()
