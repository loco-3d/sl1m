import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os
from pathlib import Path
import solo_rbprm

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.flat_ground import scene

import sl1m.tools.plot_tools as plot

"""
Validity of settings:

For `trot`:

| USE_SL1M | USE_COM | Result                           |
| -------- | ------- | -------------------------------- |
| True     | True    | FAIL                             |
| False    | False   | Result returned, but implausible |
| True     | False   | OK                               |
| False    | True    | FAIL                             |
"""

USE_SL1M = True
USE_COM = False
gait = "trot"

GAITS = {}
GAITS["walk"] = [
    np.array([1, 0, 1, 1]),
    np.array([1, 1, 0, 1]),
    np.array([0, 1, 1, 1]),
    np.array([1, 1, 1, 0]),
]
GAITS["trot"] = [np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1])]
GAITS["jumping_trot"] = [
    np.array([1, 0, 1, 0]),
    np.zeros(4),
    np.array([0, 1, 0, 1]),
    np.zeros(4),
]
GAIT = GAITS[gait]

if gait == "jumping_trot":
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import (
        jumping_trot_surfaces as surfaces,
    )
elif gait == "trot":
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import (
        trot_surfaces as surfaces,
    )
else:
    from sl1m.stand_alone_scenarios.surfaces.flat_ground import (
        walk_surfaces as surfaces,
    )

STEP_LENGTH = [0.1, 0.0]
COSTS = {
    "step_size": [10.0, STEP_LENGTH],
    "posture": [1.0],
    # "final_com": [10.0, [2.5, 0., 0.241]]
}

solo_rbprm_path = (
    Path(solo_rbprm.__file__).resolve().parent.parent.parent.parent.parent
    / "share"
    / "solo-rbprm"
)
paths = [
    str(solo_rbprm_path / "com_inequalities" / "feet_quasi_flat") + os.sep,
    str(solo_rbprm_path / "relative_effector_positions") + os.sep,
]
others = ["HR_FOOT", "HL_FOOT", "FL_FOOT", "FR_FOOT"]
limbs = ["HRleg", "HLleg", "FLleg", "FRleg"]
offsets = {
    "FRleg": [0.1946, -0.0875, -0.241],
    "FLleg": [0.1946, 0.0875, -0.241],
    "HRleg": [-0.1946, -0.0875, -0.241],
    "HLleg": [-0.1946, 0.0875, -0.241],
}

if __name__ == "__main__":
    t_init = clock()
    R = [np.identity(3)] * len(surfaces)
    q_init = [1.0, 0.0, 0.241]
    initial_contacts = [
        np.array(q_init) + np.array(offsets[limb]) for limb in limbs
    ]
    t_1 = clock()

    pb = Problem(limb_names=limbs, other_names=others, constraint_paths=paths)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init[:3])
    t_2 = clock()

    if USE_SL1M:
        result = solve_L1_combinatorial(pb, costs=COSTS, com=USE_COM)
    else:
        result = solve_MIP(pb, costs=COSTS, com=USE_COM)
    t_end = clock()

    print(result)
    print("Optimized number of steps:              ", pb.n_phases)
    print(
        "Total time is:                          ", 1000.0 * (t_end - t_init)
    )
    print("Initializing the problem takes          ", 1000.0 * (t_1 - t_init))
    print("Generating the problem dictionary takes ", 1000.0 * (t_2 - t_1))
    print("Solving the problem takes               ", 1000.0 * (t_end - t_2))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(scene)
    if result.success:
        plot.plot_planner_result(
            result.all_feet_pos,
            coms=result.coms,
            step_size=STEP_LENGTH,
            ax=ax,
            show=True,
        )
    else:
        plt.show()
