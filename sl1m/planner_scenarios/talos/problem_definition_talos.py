from sl1m.constants_and_tools import *
import sl1m.stand_alone_scenarios
from sl1m.constants_and_tools import default_transform_from_pos_normal, convert_surface_to_inequality
from sl1m.tools.obj_to_constraints import load_obj, as_inequalities, rotate_inequalities
import numpy as np
import os

# General sl1m problem definition
#
# pb["n_effectors"]: number of effectors
# pb["p0"]:      initial feet positions
# pb["c0"]:      initial com positions
# pb["nphases"]: number of phases
# pb["phaseData"][i]["Moving"]: moving effector in phase i
# pb["phaseData"][i]["K"]: Com constraints for phase i, for each limb and each surface
# pb["phaseData"][i]["allRelativeK"]: Relative constraints for phase i for each limb and each surface
# pb["phaseData"][i]["rootOrientation"]: root orientation for phase i
# pb["phaseData"][i]["S"]: surfaces of phase i

LIMB_NAMES = ["LF", "RF"]
DIR = os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/"

Z = array([0., 0., 1.])
LF = 0
RF = 1


def generate_problem(Robot, R, surfaces, gait, p0, c0=None, eq_as_ineq=True):
    """
    Build a SL1M problem for the Mixed Integer formulation,
    with all the kinematics and foot relative position constraints required
    :param Robot: an rbprm robot
    :param R: a list of rotation matrix for the base of the robot (must be the same size as surfaces)
    :param surfaces: A list of surfaces candidates, with one set of surface candidates for each phase
    :param gait: The gait of the robot (list of id of the moving foot)
    :param p0: The initial positions of the limbs
    :param c0: The initial position of the com
    :return: a "res" dictionnary with the format required by SL1M
    """
    n_phases = len(surfaces)
    n_effectors = len(LIMB_NAMES)
    res = {"n_effectors": n_effectors, "p0": p0, "c0": c0, "n_phases": n_phases}
    res["phaseData"] = [{"Moving": gait[i % n_effectors],
                         "K": com_constraint(i, R, min_height=0.3),
                         "allRelativeK": relative_constraint(i, R),
                         "S": surfaces[i]} for i in range(n_phases)]

    return res


def com_constraint(index, rotations, min_height=None):
    """
    Generate the constraints on the CoM position for all the effectors
    :param Robot:
    :param rotations: the rotation to apply to the constraints
    :param normals: the default contact normals of each effectors
    :return: a list of [A,b] inequalities, in the form Ax <= b
    """
    if index == 0:
        trLF = default_transform_from_pos_normal(zeros(3), Z, rotations[index])
        trRF = trLF.copy()
    else:
        trLF = default_transform_from_pos_normal(zeros(3), Z, rotations[index-(index % 2)])
        trRF = default_transform_from_pos_normal(zeros(3), Z, rotations[index-((index + 1) % 2)])

    KLF = com_in_effector_frame_constraint(trLF, LF)
    KRF = com_in_effector_frame_constraint(trRF, RF)
    if min_height is not None:
        KLF = [np.vstack([KLF[0], -Z]), np.concatenate([KLF[1], -
                                                        np.ones(1) * min_height]).reshape((-1,))]
        KRF = [np.vstack([KRF[0], -Z]), np.concatenate([KRF[1], -
                                                        np.ones(1) * min_height]).reshape((-1,))]
    return [KLF, KRF]


def com_in_effector_frame_constraint(transform, foot):
    """
    Generate the inequalities constraints for the CoM position given a contact position for one limb
    :param Robot:
    :param transform: Transformation to apply to the constraints
    :param foot: the Id of the limb used (see Robot.limbs_names list)
    :return: [A, b] the inequalities, in the form Ax <= b
    """
    filekin = DIR + "COM_constraints_in_" + LIMB_NAMES[foot] + "_effector_frame.obj"
    obj = load_obj(filekin)
    ineq_right_foot = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2, 3] += 0.105
    ine = rotate_inequalities(ineq_right_foot, transform2)
    return (ine.A, ine.b)


def relative_constraint(index, rotation):
    """
    Generate all the relative position constraints for all limbs
    :param Robot:
    :param rotation: the rotation to apply to the constraints
    :param normals: the default contact normals of each effectors
    :return: a list of [A,b] inequalities, in the form Ax <= b
    """
    if index == 0:
        trLF = default_transform_from_pos_normal(zeros(3), Z, rotation[index])
        trRF = trLF.copy()
    else:
        trLF = default_transform_from_pos_normal(zeros(3), Z, rotation[index-(index % 2)])
        trRF = default_transform_from_pos_normal(zeros(3), Z, rotation[index-((index + 1) % 2)])
    res = [None, None]
    res[LF] = [(RF, foot_in_other_foot_frame_constraint(trRF, LF, RF))]
    res[RF] = [(LF, foot_in_other_foot_frame_constraint(trLF, RF, LF))]
    return res


def foot_in_other_foot_frame_constraint(transform, foot, other):
    """
    Generate the constraints for the position of a given effector, given another effector position
    :param transform: The transform to apply to the constraints
    :param foot: 
    :param other: 
    :return: [A, b] the inequalities, in the form Ax <= b
    """
    filekin = DIR + LIMB_NAMES[foot] + "_constraints_in_" + LIMB_NAMES[other] + ".obj"
    obj = load_obj(filekin)
    ineq_rf_in_rl = as_inequalities(obj)
    ine = rotate_inequalities(ineq_rf_in_rl, transform)
    return (ine.A, ine.b)
