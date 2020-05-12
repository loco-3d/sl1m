from sl1m.constants_and_tools import * 
import sl1m.stand_alone_scenarios

import os

insdir = os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/"

# load  constraints
__ineq_right_foot = None
__ineq_left_foot  = None
__ineq_right_foot_reduced = None
__ineq_left_foot_reduced  = None

# add foot offset
def right_foot_constraints(transform):
    global __ineq_right_foot
    if __ineq_right_foot is None:
        filekin = insdir +"COM_constraints_in_RF_effector_frame.obj"
        obj = load_obj(filekin)
        __ineq_right_foot = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_right_foot, transform2)
    return (ine.A, ine.b)
            
def left_foot_constraints(transform):
    global __ineq_left_foot
    if __ineq_left_foot is None:
        filekin = insdir +"COM_constraints_in_LF_effector_frame.obj"
        obj = load_obj(filekin)
        __ineq_left_foot = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_left_foot, transform2)
    return (ine.A, ine.b)
    

__ineq_rf_in_rl = None
__ineq_lf_in_rf  = None
    
# add foot offset
def right_foot_in_lf_frame_constraints(transform):
    global __ineq_rf_in_rl
    if __ineq_rf_in_rl is None:
        filekin = insdir +"RF_constraints_in_LF.obj"
        obj = load_obj(filekin)
        __ineq_rf_in_rl = as_inequalities(obj)
    transform2 = transform.copy()
    ine = rotate_inequalities(__ineq_rf_in_rl, transform2)
    return (ine.A, ine.b)
        
def left_foot_in_rf_frame_constraints(transform):
    global __ineq_lf_in_rf
    if __ineq_lf_in_rf is None:
        filekin = insdir +"LF_constraints_in_RF.obj"
        obj = load_obj(filekin)
        __ineq_lf_in_rf = as_inequalities(obj)
    transform2 = transform.copy()
    ine = rotate_inequalities(__ineq_lf_in_rf, transform2)
    return (ine.A, ine.b)
    
