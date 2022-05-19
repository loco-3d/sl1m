#  General sl1m problem definition
#
# pb.n_effectors = number of effectors
# pb.p0 = initial feet positions
# pb.c0 = initial com positions
# pb.nphases = number of phases
# pb.phaseData[i].moving =  moving effector in phase i
# pb.phaseData[i].K =  Com constraints for phase i, for each limb and each surface
# pb.phaseData[i].allRelativeK =  Relative constraints for phase i for each limb and each surface
# pb.phaseData[i].root_orientation =  root orientation for phase i
# pb.phaseData[i].S =  surfaces of phase i

from sl1m.constants_and_tools import default_transform_from_pos_normal
from sl1m.tools.obj_to_constraints import load_obj, as_inequalities, rotate_inequalities
import numpy as np
import os
import sl1m.stand_alone_scenarios

LIMB_NAMES = ["LF", "RF"]
Z = np.array([0.0, 0.0, 1.0])
LF = 0
RF = 1
DIR = os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/"


class PhaseData:
    def __init__(
        self, i, R, surfaces, moving_foot, normal, n_effectors, com_obj, foot_obj
    ):
        self.moving = moving_foot
        self.root_orientation = R[i]
        self.S = surfaces
        self.n_surfaces = len(self.S)
        self.generate_K(i, R, com_obj)
        self.generate_relative_K(i, R, foot_obj)

    def generate_K(self, i, R, com_obj):
        """
        Generate the constraints on the CoM position for all the effectors as a list of [A,b]
        inequalities, in the form Ax <= b
        :param n_effectors:
        :param obj: com constraint objects
        """
        if i == 0:
            trLF = default_transform_from_pos_normal(np.zeros(3), Z, R[i])
            trRF = trLF.copy()
        else:
            trLF = default_transform_from_pos_normal(np.zeros(3), Z, R[i - (i % 2)])
            trRF = default_transform_from_pos_normal(
                np.zeros(3), Z, R[i - ((i + 1) % 2)]
            )

        trLF[2, 3] += 0.105
        ine = rotate_inequalities(com_obj[LF], trLF)
        KLF = (ine.A, ine.b)

        trRF[2, 3] += 0.105
        ine = rotate_inequalities(com_obj[RF], trRF)
        KRF = (ine.A, ine.b)

        KLF = [
            np.vstack([KLF[0], -Z]),
            np.concatenate([KLF[1], -np.ones(1) * 0.3]).reshape((-1,)),
        ]
        KRF = [
            np.vstack([KRF[0], -Z]),
            np.concatenate([KRF[1], -np.ones(1) * 0.3]).reshape((-1,)),
        ]

        self.K = [KLF, KRF]

    def generate_relative_K(self, i, R, foot_obj):
        """
        Generate all the relative position constraints for all limbs as a list of [A,b]
        inequalities, in the form Ax <= b
        :param n_effectors:
        :param obj: foot relative constraints
        """
        if i == 0:
            trLF = default_transform_from_pos_normal(np.zeros(3), Z, R[i])
            trRF = trLF.copy()
        else:
            trLF = default_transform_from_pos_normal(np.zeros(3), Z, R[i - (i % 2)])
            trRF = default_transform_from_pos_normal(
                np.zeros(3), Z, R[i - ((i + 1) % 2)]
            )

        ineRF = rotate_inequalities(foot_obj[LF], trLF)
        ineLF = rotate_inequalities(foot_obj[RF], trRF)

        self.allRelativeK = [None, None]
        self.allRelativeK[LF] = [(RF, (ineLF.A, ineLF.b))]
        self.allRelativeK[RF] = [(LF, (ineRF.A, ineRF.b))]


class Problem:
    def __init__(self):
        self.n_effectors = 2

        self.com_objects = []
        self.foot_objects = []
        for foot in range(self.n_effectors):
            f = (
                DIR
                + "COM_constraints_in_"
                + LIMB_NAMES[foot]
                + "_effector_frame_quasi_static_reduced.obj"
            )
            self.com_objects.append(as_inequalities(load_obj(f)))

            f = (
                DIR
                + LIMB_NAMES[(foot + 1) % 2]
                + "_constraints_in_"
                + LIMB_NAMES[foot]
                + "_reduced.obj"
            )
            self.foot_objects.append(as_inequalities(load_obj(f)))

    def generate_problem(self, R, surfaces, gait, p0, c0=None):
        """
        Build a SL1M problem for the Mixed Integer formulation,
        with all the kinematics and foot relative position constraints required
        :param Robot: an rbprm robot
        :param R: a list of rotation matrix for the base of the robot (must be the same size as surfaces)
        :param surfaces: A list of surfaces candidates, with one set of surface candidates for each phase
        :param gait: The gait of the robot (list of id of the moving foot)
        :param p0: The initial positions of the limbs
        :param c0: The initial position of the com
        :return: a "res" dictionary with the format required by SL1M
        """
        normal = np.array([0, 0, 1])
        self.p0 = p0
        self.c0 = c0
        self.n_phases = len(surfaces)
        self.phaseData = []
        for i in range(self.n_phases):
            self.phaseData.append(
                PhaseData(
                    i,
                    R,
                    surfaces[i],
                    gait[i % self.n_effectors],
                    normal,
                    self.n_effectors,
                    self.com_objects,
                    self.foot_objects,
                )
            )
