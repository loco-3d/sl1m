import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter as clock
import os

from sl1m.generic_solver import solve_L1_combinatorial, solve_MIP
from sl1m.problem_definition import Problem
# from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import quadruped_surfaces_gait as surfaces
from sl1m.stand_alone_scenarios.surfaces.stair_surfaces import scene

import sl1m.tools.plot_tools as plot
from example_robot_data.robots_loader import ANYmalLoader
import pinocchio as pin

# Problem definition
paths = [os.environ["INSTALL_HPP_DIR"] + "/share/anymal-rbprm/com_inequalities/feet_quasi_flat/anymal_",
         os.environ["INSTALL_HPP_DIR"] + "/share/anymal-rbprm/relative_effector_positions/anymal_"]
suffix_com = "_effector_frame_quasi_static_reduced.obj"
limbs = ['LFleg', 'RFleg', 'LHleg', 'RHleg']
others = ['LF_ADAPTER_TO_FOOT', 'RF_ADAPTER_TO_FOOT', 'LH_ADAPTER_TO_FOOT', 'RH_ADAPTER_TO_FOOT']
suffix_feet = "_reduced.obj"

# Parameters of the optimisation
USE_COM = False
# GAIT = np.array([[0.,1.,1.,1.],
                #  [1.,0.,1.,1.],
                #  [1.,1.,0.,1.],
                #  [1.,1.,1.,0.]])
GAIT = np.array([[1.,0.,1.,0.],
                 [0.,1.,0.,1.]])
n_gait = 2  # Number of phase in the gait
N_phase = 2 * n_gait  # Global number of phase

all_surfaces = []
for sf in scene:
    all_surfaces.append(sf[0])

offsets_feet = np.array([[0.37, 0.37, -0.37, -0.37],
                         [1.98572559e-01, -1.98572559e-01, 1.98572559e-01, -1.98572559e-01],
                         [2.13273153e-06, 2.13273153e-06, 2.13273153e-06, 2.13273153e-06]])


def get_potential_surfaces(configs, gait, all_surfaces):
    """
        Get the rotation matrix and surface condidates for each configuration in configs
        :param configs: a list of successive configurations of the robot
        :param gait: a gait matrix
        :return: a list of surface candidates
        """
    surfaces_list = []
    empty_list = False
    for id, config in enumerate(configs):
        foot_surfaces = []
        stance_feet = np.nonzero(gait[id % len(gait)] == 1)[0]
        previous_swing_feet = np.nonzero(gait[(id - 1) % len(gait)] == 0)[0]
        moving_feet = stance_feet[np.in1d(stance_feet, previous_swing_feet, assume_unique=True)]

        foot_surfaces.append(all_surfaces)
        surfaces_list.append(foot_surfaces)

    return surfaces_list, empty_list


def compute_com_positions(configs, pb):
    """
    Compute the com positions
    :param configs the list of configurations
    """
    com_positions = []
    for phase in pb.phaseData:
        com = configs[phase.id][:3]
        com_positions.append(configs[phase.id][:3])

    return com_positions


def compute_effector_positions(configs, bvref, pb):
    """
    Compute the desired effector positions
    :param configs the list of configurations
    :param bvref, Array (x3) the desired velocity in base frame
    """
    t_stance = T_gait / n_gait
    effector_positions = np.zeros((4, pb.n_phases, 2))

    for phase in pb.phaseData:
        for foot in phase.moving:
            rpy = pin.rpy.matrixToRpy(pin.Quaternion(np.array(configs[phase.id][3:7])).toRotationMatrix())
            yaw = rpy[2]  # Get yaw for the predicted configuration
            shoulders = np.zeros(2)
            rpy[2] = 0.  # Yaw = 0. in horizontal frame
            Rp = pin.rpy.rpyToMatrix(rpy)[:2, :2]
            heuristic = 0.5 * t_stance * Rp @ bvref[:2] + Rp @ offsets_feet[:2, foot]

            # Compute heuristic in world frame, rotation
            shoulders[0] = heuristic[0] * np.cos(yaw) - heuristic[1] * np.sin(yaw)
            shoulders[1] = heuristic[0] * np.sin(yaw) + heuristic[1] * np.cos(yaw)
            effector_positions[foot][phase.id] = np.array(configs[phase.id][:2] + shoulders)

    return effector_positions


if __name__ == '__main__':
    t_init = clock()
    t_1 = clock()

    q_init = [-0.5, 0.2, 0.47]
    initial_config = np.zeros(7)  # Initial config
    initial_config[-1] = 1.
    initial_config[:3] = np.array(q_init)

    T_gait = 2.
    bvref = np.array([0.3, 0., 0.])  # Reference velocity

    # List of configurations in planned horizon, using the reference velocity.
    configs = []
    configs.append(initial_config.tolist())
    for i in range(1, N_phase):
        config = np.zeros(7)
        config[:3] = bvref * (T_gait / n_gait) * i + initial_config[:3]
        rpy = np.array([0., 0., 0.])
        config[3:] = pin.Quaternion(pin.rpy.rpyToMatrix(rpy)).coeffs()
        configs.append(config.tolist())

    R = [pin.XYZQUATToSE3(np.array(config)).rotation for config in configs]

    q_tmp = np.zeros(3)
    q_tmp[:2] = q_init[:2]
    initial_contacts = [q_tmp + offsets_feet[:, foot] for foot in range(4)]

    surfaces, empty_list = get_potential_surfaces(configs, GAIT, all_surfaces)
    t_2 = clock()

    pb = Problem(limb_names=limbs,
                 other_names=others,
                 constraint_paths=paths,
                 suffix_com=suffix_com,
                 suffix_feet=suffix_feet)
    pb.generate_problem(R, surfaces, GAIT, initial_contacts, configs[0][:3], com=USE_COM)

    # Generate costs
    com_positions = compute_com_positions(configs, pb)
    effector_positions = compute_effector_positions(configs, bvref, pb)
    final_COM = np.array([1.5,0.,0.4+0.47])
    costs = {"effector_positions": [0.1, effector_positions]}
    t_3 = clock()

    result = solve_MIP(pb, costs=costs, com=USE_COM)
    t_end = clock()

    print(result)
    print("Optimized number of steps:              ", pb.n_phases)
    print("Total time is:                          ", 1000. * (t_end - t_init))
    print("Computing the surfaces takes            ", 1000. * (t_1 - t_init))
    print("Computing the initial contacts takes    ", 1000. * (t_2 - t_1))
    print("Generating the problem dictionary takes ", 1000. * (t_3 - t_2))
    print("Solving the problem takes               ", 1000. * (t_end - t_3))
    print("The LP and QP optimizations take        ", result.time)

    ax = plot.draw_scene(scene)
    if (result.success):
        plot.plot_planner_result(result.all_feet_pos, coms=result.coms, ax=ax, show=True)
    else:
        plt.show()