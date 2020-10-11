from multicontact_api import ContactSequence, ContactPhase, ContactPatch
from pinocchio import SE3, Quaternion
import numpy as np
from numpy.linalg import norm
import mlp.utils.util as utils

# Hardcoded data for solo !
rLegId = 'talos_rleg_rom'
rleg = 'leg_right_1_joint'
rfoot = 'leg_right_sole_fix_joint'

lLegId = 'talos_lleg_rom'
lleg = 'leg_left_1_joint'
lfoot = 'leg_left_sole_fix_joint'

rArmId = 'talos_rarm_rom'
rarm = 'arm_right_1_joint'
rhand = 'arm_right_7_joint'

lArmId = 'talos_larm_rom'
larm = 'arm_left_1_joint'
lhand = 'arm_left_7_joint'

limbs_names = [rLegId,lLegId]#[rLegId, lLegId]#[rArmId, rLegId, lArmId, lLegId]  # List of effector used to create contact
dict_limb_joint = {rLegId: rfoot, lLegId: lfoot}#, rArmId: rhand, lArmId: lhand}
rLegOffset = [0., 0., 0.0]
lLegOffset = [0., 0., 0.0]
rArmOffset = [-0.01, 0., -0.154]
lArmOffset = [-0.01, 0., -0.154]
# various offset used by scripts
MRsole_offset = SE3.Identity()
MRsole_offset.translation = np.matrix(rLegOffset).T
MLsole_offset = SE3.Identity()
MLsole_offset.translation = np.matrix(lLegOffset).T
MRhand_offset = SE3.Identity()
MRhand_offset.translation = np.matrix(rArmOffset).T
MLhand_offset = SE3.Identity()
MLhand_offset.translation = np.matrix(lArmOffset).T
dict_offset = {rfoot: MRsole_offset, lfoot: MLsole_offset}#, rhand: MRhand_offset, lhand: MLhand_offset}
DEFAULT_COM_HEIGHT = 0.86


def normal_from_ineq(s_ineq):
    return s_ineq[2]


def rotationFromNormal(n):
    """
    return a rotation matrix corresponding to the given contact normal
    :param n: the normal of the surface (as a numpy array)
    :return: a rotation matrix
    """
    z_up = np.array([0., 0., 1.])
    q = Quaternion.FromTwoVectors(z_up, n)
    return q.matrix()


def computeCenterOfSupportPolygonFromPhase(phase, DEFAULT_HEIGHT):
    """
    Compute 3D point in the center of the support polygon and at a given height
    :param phase: the ContactPhase used to compute the support polygon
    :param DEFAULT_HEIGHT: the height of the returned point
    :return: a numpy array of size 3
    """
    com = np.zeros(3)
    for patch in phase.contactPatches().values():
        com += patch.placement.translation
    com /= phase.numContacts()
    com[2] += DEFAULT_HEIGHT
    return com

def generateConfigFromPhase(fb, phase, projectCOM=False):
    """
    Compute a wholebody configuration corresponding to the contact placements defined in the given phase
    :param fb: an rbprm.FullBody object
    :param phase: the phase used to get the contact placements
    :param projectCOM: if True, the CoM position of the configuration is projected toward the center of the support polygon
    :return: a list of joint position
    """
    fb.usePosturalTaskContactCreation(False)
    effectorsInContact = phase.effectorsInContact()
    contacts = []  # contacts should contains the limb names, not the effector names
    list_effector = list(fb.dict_limb_joint.values())
    for eeName in effectorsInContact:
        contacts += [list(fb.dict_limb_joint.keys())[list_effector.index(eeName)]]
    #q = phase.q_init.tolist() # should be the correct config for the previous phase, if used only from high level helper methods
    q = fb.referenceConfig[::] + [0] * 6  # FIXME : more generic !
    root = computeCenterOfSupportPolygonFromPhase(phase, fb.DEFAULT_COM_HEIGHT).tolist()
    q[0:2] = root[0:2]
    q[2] += root[2] - fb.DEFAULT_COM_HEIGHT
    quat = Quaternion(phase.root_t.evaluateAsSE3(phase.timeInitial).rotation)
    q[3:7] = [quat.x, quat.y, quat.z, quat.w]
    # create state in fullBody :
    state = State(fb, q=q, limbsIncontact=contacts)
    # check if q is consistent with the contact placement in the phase :
    fb.setCurrentConfig(q)
    for limbId in contacts:
        eeName = fb.dict_limb_joint[limbId]
        placement_fb = SE3FromConfig(fb.getJointPosition(eeName))
        placement_phase = phase.contactPatch(eeName).placement
        if placement_fb != placement_phase:  # add a threshold instead of 0 ? how ?
            # need to project the new contact :
            placement = phase.contactPatch(eeName).placement
            p = placement.translation.tolist()
            n = utils.computeContactNormal(placement).tolist()
            state, success = StateHelper.addNewContact(state, limbId, p, n, 1000)
            if not success:
                print("Cannot project the configuration to contact, for effector : ", eeName)
                return state.q()
            if projectCOM:
                success = utils.projectCoMInSupportPolygon(state)
                if not success:
                    print("cannot project com to the middle of the support polygon.")
    phase.q_init = np.array(state.q())

    return state.q()

def placement_from_sl1m(ee_name, pos, phase_data):
    #pos[2] += EPS_Z # FIXME: apply epsz along the normal
    pos = dict_offset[ee_name].actInv(pos)
    print("Move effector ", ee_name)
    print("To position ", pos)
    placement = SE3.Identity()
    placement.translation = pos
    # compute orientation of the contact from the surface normal:
    # n = normal_from_ineq(phase_data["S"][0])
    n = np.array([0,0,1])
    placement.rotation = rotationFromNormal(n)
    print("new contact placement : ", placement)
    # TODO add yaw rotation from guide here !
    return placement


def build_cs_from_sl1m_mip(pb, allfeetpos, fb, q_init):
    # init contact sequence with first phase : q_ref move at the right root pose and with both feet in contact
    # FIXME : allow to customize that first phase
    num_steps = len(pb["phaseData"]) - 1 # number of contact repositionning
    num_effectors = len(limbs_names)
    print(" limbs names : ", limbs_names)
    cs = ContactSequence(0)
    # create the first contact phase :
    cp_init = ContactPhase()
    # cp_init.q_init = np.array(q_init)
    for k, pos in enumerate(allfeetpos[0]):
        phase_data = pb["phaseData"][0]
        ee_name = dict_limb_joint[limbs_names[k]]
        cp_init.addContact(ee_name, ContactPatch(placement_from_sl1m(ee_name, pos, phase_data)))
    cs.append(cp_init)
    # cp_init.q_init = generateConfigFromPhase(fb, cs.contactPhases[0])
    # print("Initial phase added, contacts : ", cs.contactPhases[0].effectorsInContact())
    # loop for all effector placements, and create the required contact phases
    # previous_eff_placements = allfeetpos[0]
    # if len(previous_eff_placements) != num_effectors:
    #     raise NotImplementedError("A phase in the output of SL1M do not have all the effectors in contact.")
    for pid, eff_placements in enumerate(allfeetpos[1:]):
        print("Loop allfeetpos, id = ", pid)
        # if len(eff_placements) != num_effectors:
        #     raise NotImplementedError("A phase in the output of SL1M do not have all the effectors in contact.")
        # switch = False # True if a repostionning have been detected
        # for k, pos in enumerate(eff_placements):
        #     if norm(pos - previous_eff_placements[k]) > 1e-3:
        #         if switch:
        #             raise NotImplementedError("Several contact changes between two adjacent phases in SL1M output")
        # switch = True
        pos = eff_placements#[(pid)%2]
        ee_name = dict_limb_joint[limbs_names[(pid)%2]]
        phase_data = pb["phaseData"][pid+1] # +1 because the for loop start at id = 1
        placement = placement_from_sl1m(ee_name, pos, phase_data)
        cs.moveEffectorToPlacement(ee_name, placement)
        # if not switch:
        #    raise RuntimeError("No contact changes between two adjacent phases in SL1M output")
        # assign com position to the last two phases :
        # swinging phase, same init and final position
        """
        cs.contactPhases[-2].c_init = coms[pid * 2]
        cs.contactPhases[-2].c_final = coms[pid * 2 + 1]
        # phase with all the contacts:
        cs.contactPhases[-1].c_init = coms[pid * 2 + 1]
        if pid * 2 + 2 < len(coms):
            cs.contactPhases[-1].c_final = coms[pid * 2 + 2]
        else:
            cs.contactPhases[-1].c_final = cs.contactPhases[-1].c_init
        """
        # previous_eff_placements = eff_placements

    #final
    # pos = eff_placements[(pid+1)%2]
    if pid%2 == 0: # right
        pos = pos + np.array([0,0.085*2,0])
    else : # left
        pos = pos - np.array([0,0.085*2,0])
    ee_name = dict_limb_joint[limbs_names[(pid+1)%2]]
    phase_data = pb["phaseData"][pid] # +1 because the for loop start at id = 1
    placement = placement_from_sl1m(ee_name, pos, phase_data)
    cs.moveEffectorToPlacement(ee_name, placement)

    p_final = cs.contactPhases[-1]
    p_final.q_end = np.array(q_end)
    p_final.c_final = computeCenterOfSupportPolygonFromPhase(p_final, DEFAULT_COM_HEIGHT)
    p_final.c_init = p_final.c_final
    # cs.contactPhases[-1]=p_final
    return cs

