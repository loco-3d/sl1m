from sl1m.constants_and_tools import *

import os

if "INSTALL_HPP_DIR" in os.environ:
    insdir = os.environ["INSTALL_HPP_DIR"] + "/share/"
elif "DEVEL_HPP_DIR" in os.environ:
    insdir = os.environ["DEVEL_HPP_DIR"] + "/install/share/"

# load hrp2 constraints
__ineq_right_foot_hrp2 = None
__ineq_left_foot_hrp2 = None
__ineq_right_foot_hrp2_reduced = None
__ineq_left_foot_hrp2_reduced = None

# add foot offset
def right_foot_hrp2_constraints(transform):
    global __ineq_right_foot_hrp2
    if __ineq_right_foot_hrp2 is None:
        filekin = (
            insdir
            + "hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_RF_effector_frame_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_right_foot_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2, 3] += 0.105
    ine = rotate_inequalities(__ineq_right_foot_hrp2, transform2)
    return (ine.A, ine.b)


def left_foot_hrp2_constraints(transform):
    global __ineq_left_foot_hrp2
    if __ineq_left_foot_hrp2 is None:
        filekin = (
            insdir
            + "hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_LF_effector_frame_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_left_foot_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2, 3] += 0.105
    ine = rotate_inequalities(__ineq_left_foot_hrp2, transform2)
    return (ine.A, ine.b)


def right_foot_talos_constraints(transform):
    global __ineq_right_foot_hrp2
    if __ineq_right_foot_hrp2 is None:
        filekin = (
            insdir
            + "talos-rbprm/com_inequalities/feet_quasi_flat/talos_COM_constraints_in_RF_effector_frame_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_right_foot_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_right_foot_hrp2, transform.copy())
    return (ine.A, ine.b)


def left_foot_talos_constraints(transform):
    global __ineq_left_foot_hrp2
    if __ineq_left_foot_hrp2 is None:
        filekin = (
            insdir
            + "talos-rbprm/com_inequalities/feet_quasi_flat/talos_COM_constraints_in_LF_effector_frame_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_left_foot_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_left_foot_hrp2, transform.copy())
    return (ine.A, ine.b)


__ineq_rf_in_rl_hrp2 = None
__ineq_lf_in_rf_hrp2 = None

# add foot offset
def right_foot_in_lf_frame_hrp2_constraints(transform):
    global __ineq_rf_in_rl_hrp2
    if __ineq_rf_in_rl_hrp2 is None:
        filekin = (
            insdir
            + "hrp2-rbprm/relative_effector_positions/hrp2_RF_constraints_in_LF_quasi_flat_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_rf_in_rl_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    ine = rotate_inequalities(__ineq_rf_in_rl_hrp2, transform2)
    return (ine.A, ine.b)


def left_foot_in_rf_frame_hrp2_constraints(transform):
    global __ineq_lf_in_rf_hrp2
    if __ineq_lf_in_rf_hrp2 is None:
        filekin = (
            insdir
            + "hrp2-rbprm/relative_effector_positions/hrp2_LF_constraints_in_RF_quasi_flat_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_lf_in_rf_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    ine = rotate_inequalities(__ineq_lf_in_rf_hrp2, transform2)
    return (ine.A, ine.b)


# add foot offset
def right_foot_in_lf_frame_talos_constraints(transform):
    global __ineq_rf_in_rl_hrp2
    if __ineq_rf_in_rl_hrp2 is None:
        filekin = (
            insdir
            + "talos-rbprm/relative_effector_positions/talos_RF_constraints_in_LF_quasi_flat_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_rf_in_rl_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_rf_in_rl_hrp2, transform.copy())
    return (ine.A, ine.b)


def left_foot_in_rf_frame_talos_constraints(transform):
    global __ineq_lf_in_rf_hrp2
    if __ineq_lf_in_rf_hrp2 is None:
        filekin = (
            insdir
            + "talos-rbprm/relative_effector_positions/talos_LF_constraints_in_RF_quasi_flat_REDUCED.obj"
        )
        obj = load_obj(filekin)
        __ineq_lf_in_rf_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_lf_in_rf_hrp2, transform.copy())
    return (ine.A, ine.b)
