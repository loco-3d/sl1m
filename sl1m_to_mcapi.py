from multicontact_api import ContactSequence, ContactPhase, ContactPatch
from pinocchio import SE3, Quaternion
import numpy as np
from numpy.linalg import norm

# Hardcoded data for solo !
rLegId = 'FRleg'
rleg = 'FR_HAA'
rfoot = 'FR_FOOT'
lLegId = 'FLleg'
lleg = 'FL_HAA'
lfoot = 'FL_FOOT'
lArmId = 'HLleg'
larm = 'HL_HAA'
lhand = 'HL_FOOT'
rArmId = 'HRleg'
rarm = 'HR_HAA'
rhand = 'HR_FOOT'
limbs_names = [rArmId, rLegId, lArmId, lLegId]  # List of effector used to create contact
dict_limb_joint = {rLegId: rfoot, lLegId: lfoot, rArmId: rhand, lArmId: lhand}
offset = [0., 0., -0.018]  # Contact position in the effector frame
# various offset used by scripts
MRsole_offset = SE3.Identity()
MRsole_offset.translation = np.matrix(offset).T
MLsole_offset = MRsole_offset.copy()
MRhand_offset = MRsole_offset.copy()
MLhand_offset = MRsole_offset.copy()
dict_offset = {rfoot: MRsole_offset, lfoot: MLsole_offset, rhand: MRhand_offset, lhand: MLhand_offset}
DEFAULT_COM_HEIGHT = 0.2

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

def placement_from_sl1m(ee_name, pos, phase_data):
    #pos[2] += EPS_Z # FIXME: apply epsz along the normal
    pos = dict_offset[ee_name].actInv(pos)
    print("Move effector ", ee_name)
    print("To position ", pos)
    placement = SE3.Identity()
    placement.translation = pos
    # compute orientation of the contact from the surface normal:
    n = normal_from_ineq(phase_data["S"][phase_data["id_surface"]])
    placement.rotation = rotationFromNormal(n)
    print("new contact placement : ", placement)
    # TODO add yaw rotation from guide here !
    return placement


def build_cs_from_sl1m_mip(pb, allfeetpos):
    # init contact sequence with first phase : q_ref move at the right root pose and with both feet in contact
    # FIXME : allow to customize that first phase
    num_steps = len(pb["phaseData"]) - 1 # number of contact repositionning
    num_effectors = len(limbs_names)
    print(" limbs names : ", limbs_names)
    cs = ContactSequence(0)
    # create the first contact phase :
    cp_init = ContactPhase()
    for k, pos in enumerate(allfeetpos[0]):
        phase_data = pb["phaseData"][0]
        ee_name = dict_limb_joint[limbs_names[k]]
        cp_init.addContact(ee_name, ContactPatch(placement_from_sl1m(ee_name, pos, phase_data)))
    cs.append(cp_init)
    print("Initial phase added, contacts : ", cs.contactPhases[0].effectorsInContact())
    # loop for all effector placements, and create the required contact phases
    previous_eff_placements = allfeetpos[0]
    if len(previous_eff_placements) != num_effectors:
        raise NotImplementedError("A phase in the output of SL1M do not have all the effectors in contact.")
    for pid, eff_placements in enumerate(allfeetpos[1:]):
        print("Loop allfeetpos, id = ", pid)
        if len(eff_placements) != num_effectors:
            raise NotImplementedError("A phase in the output of SL1M do not have all the effectors in contact.")
        switch = False # True if a repostionning have been detected
        for k, pos in enumerate(eff_placements):
            if norm(pos - previous_eff_placements[k]) > 1e-3:
                if switch:
                    raise NotImplementedError("Several contact changes between two adjacent phases in SL1M output")
                switch = True
                ee_name = dict_limb_joint[limbs_names[k]]
                phase_data = pb["phaseData"][pid+1] # +1 because the for loop start at id = 1
                placement = placement_from_sl1m(ee_name, pos, phase_data)
                cs.moveEffectorToPlacement(ee_name, placement)

        if not switch:
           raise RuntimeError("No contact changes between two adjacent phases in SL1M output")
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
        previous_eff_placements = eff_placements
    p_final = cs.contactPhases[-1]
    p_final.c_final = computeCenterOfSupportPolygonFromPhase(p_final, DEFAULT_COM_HEIGHT)
    p_final.c_init = p_final.c_final
    return cs

