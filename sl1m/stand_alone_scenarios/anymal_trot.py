import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import scene

import sl1m.tools.plot_tools as plot
import pinocchio as pin
import anymal_rbprm

# Problem definition
paths = [os.path.join(os.path.dirname(anymal_rbprm.__file__), "../../../..", "share/anymal-rbprm/com_inequalities/feet_quasi_flat/anymal_"),
         os.path.join(os.path.dirname(anymal_rbprm.__file__), "../../../..", "share/anymal-rbprm/relative_effector_positions/anymal_")]
suffix_com = "_effector_frame_quasi_static_reduced.obj"
limbs = ['LFleg', 'RFleg', 'LHleg', 'RHleg']
others = ['LF_ADAPTER_TO_FOOT', 'RF_ADAPTER_TO_FOOT', 'LH_ADAPTER_TO_FOOT', 'RH_ADAPTER_TO_FOOT']
suffix_feet = "_reduced.obj"

# Order : [LF, RF, LH, RH]
GAITS = {}
GAITS["walk"] = np.array([[0., 1., 1., 1.], [1., 0., 1., 1.], [1., 1., 0., 1.], [1., 1., 1., 0.]])
GAITS["trot"] = np.array([[1., 0., 1., 0.], [0., 1., 0., 1.]])
# Number of phase in one gait
N_GAIT = {}
N_GAIT["walk"] = 4
N_GAIT["trot"] = 2

# Parameters of the optimisation
USE_COM = True
gait = "trot"  # type of gait chosen
n_gait = N_GAIT[gait]  # Number of phase in the gait
GAIT = GAITS[gait]
T_gait = 1. # Period of the gait
N_phase = 8 * n_gait  # Global number of phase
q_init = np.array([-0.5, 0.2, 0.47])  # Initial Position
bvref = np.array([0.2, 0., 0.])  # Reference velocity

# List of surfaces following the format :
# vertices = array([[x0, ... , xn],
#                  [y0, ... , yn],
#                  [z0, ... , zn]])
all_surfaces = []
for sf in scene:
    all_surfaces.append(sf[0])

# Order : [LF, RF, LH, RH]
offsets_feet = np.array([[0.37, 0.37, -0.37, -0.37], [0.2, -0.2, 0.2, -0.2], [0., 0., 0., 0.]])


def get_potential_surfaces(configs, gait, all_surfaces):
    """
    Get the potential surfaces condidates for each moving foot.
    Here, all the surfaces are allocated for each foot, could be decreases using a guide path
    such as RBPRM.

    Args:
         - configs (list): the list of configurations.
         - gait (array nx4): gait matrix.
         - all_surfaces (list): list of surfaces.

    Returns:
        - (list): List of surface candidates for each moving foot.
    """
    surfaces_list = []
    for id, config in enumerate(configs):
        foot_surfaces = []
        stance_feet = np.nonzero(gait[id % len(gait)] == 1)[0]
        previous_swing_feet = np.nonzero(gait[(id - 1) % len(gait)] == 0)[0]
        moving_feet = stance_feet[np.in1d(stance_feet, previous_swing_feet, assume_unique=True)]
        for elt in moving_feet:
            foot_surfaces.append(all_surfaces)
        surfaces_list.append(foot_surfaces)

    return surfaces_list


def compute_com_positions(configs, pb):
    """
    Compute the CoM positions.

    Args:
        - configs (list): the list of configurations.
        - pb (Problem): sl1m problem.

    Returns:
        - (list): the list of CoM.
    """
    com_positions = []
    for phase in pb.phaseData:
        com = configs[phase.id][:3]
        com_positions.append(configs[phase.id][:3])

    return com_positions


def compute_effector_positions(configs, bvref, pb):
    """
    Compute the desired effector positions.

    Args:
        - configs (list): the list of configurations.
        - bvref (Array x3): the desired velocity in base frame.
        - pb (Problem): sl1m problem.

    Returns:
        - (list): Effector positions.
    """
    t_stance = T_gait / n_gait
    effector_positions = np.zeros((4, pb.n_phases, 2))

    for phase in pb.phaseData:
        for foot in phase.moving:
            rpy = pin.rpy.matrixToRpy(pin.Quaternion(np.array(configs[phase.id][3:7])).toRotationMatrix())
            yaw = rpy[2]  # Get yaw for the predicted configuration
            shoulders = np.zeros(2)
            # Compute heuristic position in horizontal frame
            rpy[2] = 0.  # Yaw = 0. in horizontal frame
            Rp = pin.rpy.rpyToMatrix(rpy)[:2, :2]
            heuristic = 0.5 * t_stance * Rp @ bvref[:2] + Rp @ offsets_feet[:2, foot]

            # Compute heuristic in world frame, rotation
            shoulders[0] = heuristic[0] * np.cos(yaw) - heuristic[1] * np.sin(yaw)
            shoulders[1] = heuristic[0] * np.sin(yaw) + heuristic[1] * np.cos(yaw)
            effector_positions[foot][phase.id] = np.array(configs[phase.id][:2] + shoulders)

    return effector_positions


def get_height_terrain(x, y):
    """ Get height and pitch angle of the terrain given x,y position. Hardcoded heightmap.

    Args:
        - x (float): x position.
        - y (float): y position.

    Returns:
        - (float): Z position.
        - (float): Pitch (rad).
    """
    slope_terrain = 0.4 / 2
    if x < 0.:
        return 0., 0.
    elif x > 2:
        return 0.4, 0.
    else:
        return slope_terrain * x, -np.arctan2(slope_terrain * x, x)


if __name__ == '__main__':
    t_init = clock()
    t_1 = clock()

    initial_config = np.zeros(7)  # Initial config
    initial_config[-1] = 1.
    initial_config[:3] = q_init

    # List of configurations in planned horizon, using the reference velocity.
    configs = []
    configs.append(initial_config.tolist())
    for i in range(1, N_phase):
        config = np.zeros(7)
        config[:3] = bvref * (T_gait / n_gait) * i + initial_config[:3]
        height, pitch = get_height_terrain(config[0], config[1])
        config[2] += height
        rpy = np.array([0., pitch, 0.])
        config[3:] = pin.Quaternion(pin.rpy.rpyToMatrix(rpy)).coeffs()
        configs.append(config.tolist())

    # List of rotation matrix
    R = [pin.XYZQUATToSE3(np.array(config)).rotation for config in configs]

    # Initian contacts
    q_tmp = np.zeros(3)
    q_tmp[:2] = q_init[:2] # Tmp vector, z = 0
    initial_contacts = [q_tmp + offsets_feet[:, foot] for foot in range(4)]

    surfaces = get_potential_surfaces(configs, GAIT, all_surfaces)
    t_2 = clock()

    pb = Problem(limb_names=limbs,
                 other_names=others,
                 constraint_paths=paths,
                 suffix_com=suffix_com,
                 suffix_feet=suffix_feet)
    # pb.generate_problem(R, surfaces, GAIT, initial_contacts, q_init)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, configs[0][:3], com=USE_COM)

    # Generate costs
    com_positions = compute_com_positions(configs, pb)
    effector_positions = compute_effector_positions(configs, bvref, pb)
    costs = {
        "effector_positions": [1.0, effector_positions],
        "coms_xy": [0.5, com_positions],
        "coms_z": [0.5, com_positions]
    }
    t_3 = clock()

    result = solve_MIP(pb, costs=costs, com=USE_COM)
    t_end = clock()

    print(result)
    print("Optimized number of steps:                              ", pb.n_phases)
    print("Total time is:                                          ", 1000. * (t_end - t_init))
    print("Initialization (Configs, Rot, contacts, surfaces) takes:", 1000. * (t_2 - t_1))
    print("Generating the problem dictionary takes                 ", 1000. * (t_3 - t_2))
    print("Solving the problem takes                               ", 1000. * (t_end - t_3))
    print("The LP and QP optimizations take                        ", result.time)

    ax = plot.draw_scene(scene)
    if (result.success):
        plot.plot_planner_result(result.all_feet_pos, coms=result.coms, ax=ax, show=True)
    else:
        plt.show()