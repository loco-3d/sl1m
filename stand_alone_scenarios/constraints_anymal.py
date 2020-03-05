from sl1m.constants_and_tools import * 
import sl1m.stand_alone_scenarios

import os

insdir = os.path.dirname(sl1m.stand_alone_scenarios.__file__) + "/constraints_files/anymal/"

limbNames = ['LFleg', 'RFleg', 'LHleg', 'RHleg']

effectorNames = ['LF_ADAPTER_TO_FOOT', 'RF_ADAPTER_TO_FOOT', 'LH_ADAPTER_TO_FOOT', 'RH_ADAPTER_TO_FOOT']


__ineq_com = [None for _ in range(len(limbNames))]
__ineq_relative = [ [None for _ in range(len(limbNames)-1)] for __ in range(len(limbNames))]


# add foot offset
def com_in_limb_effector_frame_constraint(transform, limbId):
    global __ineq_com
    assert (limbId <= len(effectorNames))
    if __ineq_com[limbId] is None:
        # ~ filekin = insdir +"anymal_LF_ADAPTER_TO_FOOT_constraints_in_LHleg_reduced.obj"
        filekin = insdir +"anymal_COM_constraints_in_"+limbNames[limbId]+ "_effector_frame_quasi_static_reduced.obj"
        obj = load_obj(filekin)
        __ineq_com[limbId] = as_inequalities(obj)
    transform2 = transform.copy()
    # ~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_com[limbId], transform2)
    return (ine.A, ine.b)
        
    
def realIdxFootId(limbId, footId):
    if footId > limbId:
        return footId -1 
    return footId
    
# add foot offset
def foot_in_limb_effector_frame_constraint(transform, limbId, footId):
    global __ineq_relative
    assert (limbId != footId)
    if __ineq_relative[limbId][realIdxFootId(limbId,footId)] is None:
        # ~ filekin = insdir +"RF_constraints_in_LF.obj"
        filekin = insdir +"anymal_"+effectorNames[footId]+"_constraints_in_"+limbNames[limbId]+"_reduced.obj"
        obj = load_obj(filekin)
        __ineq_relative[limbId][realIdxFootId(limbId,footId)] = as_inequalities(obj)
    transform2 = transform.copy()
    ine = rotate_inequalities(__ineq_relative[limbId][realIdxFootId(limbId,footId)], transform2)
    return (ine.A, ine.b)

def genCOMConstraints(rotation = [Id for _ in range(len(limbNames))], normals = [z for _ in range(len(limbNames))]):
    return [com_in_limb_effector_frame_constraint(default_transform_from_pos_normal_(rotation[idx], zero3, normals[idx]),idx) for idx in range(len(limbNames))]
    
def genRelativeConstraints(rotation = [Id for _ in range(len(limbNames))], normals = [z for _ in range(len(limbNames))]):
    transforms = [default_transform_from_pos_normal_(rotation[idx], zero3, normals[idx]) for idx in range(len(limbNames))]
    res = []
    for limbId, transform in enumerate(transforms):
        res += [[(footId, foot_in_limb_effector_frame_constraint(transform, limbId, footId)) for footId in range(len(limbNames)) if footId != limbId]]
    return res
