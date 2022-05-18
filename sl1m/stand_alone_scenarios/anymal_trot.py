import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from time import perf_counter as clock
import os
import anymal_rbprm

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem

import sl1m.tools.plot_tools as plot

USE_SL1M = False  # Cf. `anymal_stairs.py`
USE_COM = True
gait = "trot"

GAITS = {}
# Caracal limb swing order: LH, LF, RH, RF
GAITS["walk"] = [np.array([1, 0, 1, 1]), np.array([0, 1, 1, 1]), np.array([1, 1, 1, 0]), np.array([1, 1, 0, 1])]

# Trot pairs: LF-RH, RF-LH -- does not work
GAITS["trot"] = [np.array([1, 0, 0, 1]), np.array([0, 1, 1, 0])]

# This pattern outputs something with SL1M and CoM optimisation [Pinocchio order]
GAITS["trot"] = [np.array([1, 1, 0, 0]),
                 np.array([0, 0, 1, 1])]

# # HH, FF: Does not work [works with SL1M order + CoM optimisation]
# GAITS["trot"] = [np.array([1, 0, 1, 0]),
#                  np.array([0, 1, 0, 1])]

GAIT = GAITS[gait]

# 0.2m step length
STEP_LENGTH = [0.2, 0.0]
COSTS = {
    "step_size": [1.0, STEP_LENGTH],      # Maintain step size as close as possible to desired. Does not guarantee exceeding.
    "posture": [0.01],                      # Keep relative feet positions as close as possible to initial
    "final_com": [10.0, [0.0, 0.0, 0.47]]  # Achieve desired CoM at end of plan
}

paths = [os.path.join(os.path.dirname(anymal_rbprm.__file__), "../../../..", "share/anymal-rbprm/com_inequalities/feet_quasi_flat/anymal_"),
         os.path.join(os.path.dirname(anymal_rbprm.__file__), "../../../..", "share/anymal-rbprm/relative_effector_positions/anymal_")]
# ## SL1M order
# limbs = ['RHleg', 'RFleg', 'LHleg', 'LFleg']
# others = ['RH_ADAPTER_TO_FOOT', 'RF_ADAPTER_TO_FOOT', 'LH_ADAPTER_TO_FOOT', 'LF_ADAPTER_TO_FOOT']
## Pinocchio order
limbs = ['LFleg', 'LHleg', 'RFleg', 'RHleg']
others = ['LF_ADAPTER_TO_FOOT', 'LH_ADAPTER_TO_FOOT', 'RF_ADAPTER_TO_FOOT', 'RH_ADAPTER_TO_FOOT']
suffix = "_effector_frame_quasi_static_upscale.obj"

# The nominal configuration would be 0.37, 0.2 
#  --> i.e. this is 6.4cm wider per leg than usual --> 12.8cm total
# offsets = {'RFleg': [0.373, -0.264, -0.47], 'LFleg': [0.373, 0.264, -0.47],
#            'RHleg': [-0.373, -0.264, -0.47], 'LHleg': [-0.373, 0.264, -0.47]}
# Usual nominal config:
offsets = {'RFleg': [0.373, -0.2, -0.47], 'LFleg': [0.373, 0.2, -0.47],
           'RHleg': [-0.373, -0.2, -0.47], 'LHleg': [-0.373, 0.2, -0.47]}

if __name__ == '__main__':
    # Change this to -3 to make it fail...
    q_init = [-2.0, 0., 0.47]
    
    if "final_com" in COSTS:
        # Figure out distance
        dist = np.linalg.norm(np.asarray(COSTS["final_com"][1])-q_init)
        n_phases = 2*ceil(dist / STEP_LENGTH[0])
        print("Distance to traverse", dist, "# of phases", n_phases)
    else:
        n_phases = 20 # dummy value

    surf_size = 4
    surface = np.array([[-surf_size, -surf_size, surf_size, surf_size],
                        [surf_size, -surf_size, -surf_size, surf_size],
                        [0., 0., 0., 0.]])
    scene = [[surface], [surface]]
    surfaces = [scene] * n_phases

    t_init = clock()
    R = [np.identity(3)] * n_phases
    
    initial_contacts = [np.array(q_init) + offsets[limb] for limb in limbs]
    t_1 = clock()

    pb = Problem(limb_names=limbs, other_names=others, constraint_paths=paths, suffix_com=suffix)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init[:3])
    t_2 = clock()

    if USE_SL1M:
        result = solve_L1_combinatorial(pb, costs=COSTS,  com=USE_COM)
    else:
        result = solve_MIP(pb, costs=COSTS, com=USE_COM)
    t_end = clock()

    print(result)
    print("Optimized number of steps:              {0:>8}".format(pb.n_phases))
    print("Total time is:                          {0:>8.3f}ms".format(1000. * (t_end-t_init)))
    print("Initializing the problem takes          {0:>8.3f}ms".format(1000. * (t_1 - t_init)))
    print("Generating the problem dictionary takes {0:>8.3f}ms".format(1000. * (t_2 - t_1)))
    print("Solving the problem takes               {0:>8.3f}ms".format(1000. * (t_end - t_2)))
    print("The LP and QP optimizations take        {0:>8.3f}ms".format(result.time))

    # Analyse result:
    if result.success:
        print("\n\n\n")
        # Check final-CoM
        if "final_com" in COSTS:
            print("Desired Final-CoM = ({0:.2f}, {1:.2f}, {2:.2f})"
                  " vs Actual Final-CoM = ({3:.2f}, {4:.2f}, {5:.2f})".format(
                  *COSTS['final_com'][1], *result.coms[-1]))

        # Check footstep length
        for i in range(4):
            for j in range(1, len(result.all_feet_pos)):
                if result.all_feet_pos[i][j] is None: continue
                if result.all_feet_pos[i][j-1] is None: continue
                print("Foot #{0}, Step #{1} Dist {2:.2f}m".format(i, j, np.linalg.norm(result.all_feet_pos[i][j]-result.all_feet_pos[i][j-1])))

    if result.success:
        ax = plot.draw_scene(scene)
        plot.plot_planner_result(result.all_feet_pos, coms=result.coms, step_size=STEP_LENGTH, ax=ax, show=True)
        plt.show()
