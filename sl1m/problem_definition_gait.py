from sl1m.constants_and_tools import default_transform_from_pos_normal, convert_surface_to_inequality
from sl1m.tools.obj_to_constraints import load_obj, as_inequalities, rotate_inequalities
import numpy as np


class PhaseData:
    """
    phaseData.moving =  moving feet in phase i
    phaseData.root_orientation =  root orientation for phase i
    phaseData.S =  surfaces of phase i
    phaseData.n_surfaces =  number of surfaces
    phaseData.K =  Com constraints for the phase, for each foot and each surface
    phaseData.allRelativeK =  Relative constraints for the phase for each foot and each surface
    """

    def __init__(self, i, R, surfaces, gait, normal,  n_effectors, com_obj, foot_obj):
        self.id = i
        previous_swing_feet = np.nonzero(gait[(i-1) % len(gait)] == 0)[0]
        self.stance = np.nonzero(gait[i % len(gait)] == 1)[0]
        self.moving = self.stance[np.in1d(self.stance, previous_swing_feet, assume_unique=True)]
        self.root_orientation = R
        self.S = [[convert_surface_to_inequality(s, True) for s in foot_surfaces] for foot_surfaces in surfaces]
        self.n_surfaces = [len(s) for s in self.S]
        self.transform = default_transform_from_pos_normal(np.zeros(3), normal, R)
        self.generate_K(n_effectors, com_obj)
        self.generate_relative_K(n_effectors, foot_obj)

    def generate_K(self, n_effectors, obj):
        """
        Generate the constraints on the CoM position for all the effectors as a list of [A,b] 
        inequalities, in the form Ax <= b
        :param n_effectors:
        :param obj: com constraint objects
        """
        self.K = []
        for foot in range(n_effectors):
            ine = rotate_inequalities(obj[foot], self.transform.copy())
            self.K.append((ine.A, ine.b))

    def generate_relative_K(self, n_effectors, obj):
        """
        Generate all the relative position constraints for all limbs as a list of [A,b] 
        inequalities, in the form Ax <= b
        :param n_effectors:
        :param obj: foot relative constraints
        """
        self.allRelativeK = []
        for foot in range(n_effectors):
            foot_res = []
            for other in range(n_effectors):
                if other != foot:
                    ine = rotate_inequalities(obj[foot][other], self.transform.copy())
                    foot_res.append((other, (ine.A, ine.b)))
            self.allRelativeK += [foot_res]


class Problem:
    """
    General sl1m problem definition

    pb.n_effectors = number of effectors
    pb.p0 = initial feet positions
    pb.c0 = initial com positions
    pb.nphases = number of phases
    pb.phaseData list of Phase data objects
    """

    def __init__(self, Robot, suffix_com="_effector_frame_quasi_static_reduced.obj", suffix_feet="_reduced.obj", limb_names=None):
        if limb_names:
            self.n_effectors = len(limb_names)
        else:
            self.n_effectors = len(Robot.limbs_names)

        self.com_objects = []
        self.foot_objects = []
        for foot in range(self.n_effectors):
            if limb_names != None:
                foot_name = limb_names[foot]
            else:
                foot_name = Robot.limbs_names[foot]
            filekin = Robot.kinematic_constraints_path + "COM_constraints_in_" + foot_name + suffix_com
            self.com_objects.append(as_inequalities(load_obj(filekin)))

            foot_object = []
            for other in range(self.n_effectors):
                if other != foot:
                    if limb_names != None:
                        other_name = limb_names[other]
                    else:
                        other_name = Robot.dict_limb_joint[Robot.limbs_names[other]]
                    filekin = Robot.relative_feet_constraints_path + \
                        other_name + "_constraints_in_" + foot_name + suffix_feet
                    foot_object.append(as_inequalities(load_obj(filekin)))
                else:
                    foot_object.append(None)

            self.foot_objects.append(foot_object)

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
        :return: a "res" dictionnary with the format required by SL1M
        """
        normal = np.array([0, 0, 1])
        self.p0 = p0
        self.c0 = c0
        self.n_phases = len(surfaces)
        self.phaseData = []
        for i in range(self.n_phases):
            self.phaseData.append(PhaseData(i, R[i], surfaces[i], gait, normal,
                                            self.n_effectors, self.com_objects, self.foot_objects))

    def __str__(self):
        string = "Problem: "
        string += "\n \t n_effectors: " + str(self.n_effectors)
        string += "\n \t n_phases: " + str(self.n_phases)
        string += "\n \t p0: " + str(self.p0)
        string += "\n \t c0: " + str(self.c0)
        for i in range(self.n_phases):
            string += "\n \t \t Phase: " + str(i)
            string += "\n \t \t moving: " + str(self.phaseData[i].moving)
            string += "\n \t \t n_surfaces: " + str(self.phaseData[i].n_surfaces)
        return string
