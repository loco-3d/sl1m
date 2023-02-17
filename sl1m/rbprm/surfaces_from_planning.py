from numpy import arange, array
import numpy as np
from hpp.corbaserver.rbprm.tools.narrow_convex_hull import (
    getSurfaceExtremumPoints,
    removeDuplicates,
    normal,
)
from pinocchio import XYZQUATToSE3

MAX_SURFACE = 1.0  # if a contact surface is greater than this value, the intersection is used instead of the whole surface
LF = 0
RF = 1


def listToArray(seqs):
    nseq = []
    nseqs = []
    for seq in seqs:
        nseq = []
        for surface in seq:
            nseq.append(array(surface).T)
        nseqs.append(nseq)
    return nseqs


def area(s):
    area = 0
    for i in range(1, len(s) - 1):
        p0, p1, p2 = array(s[0]), array(s[i]), array(s[i + 1])
        area += np.linalg.norm(np.cross(p1 - p0, p2 - p0))
    return area * 0.5


def getRotationsFromConfigs(configs):
    R = []
    for config in configs:
        q_rot = config[3:7]
        R.append(XYZQUATToSE3([0, 0, 0] + q_rot).rotation)
    return R


def getContactsNames(rbprmBuilder, i, q):
    if i % 2 == LF:  # left leg
        step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(
            q, "talos_lleg_rom"
        )
    elif i % 2 == RF:  # right leg
        step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(
            q, "talos_rleg_rom"
        )
    return step_contacts


def getContactsIntersections(rbprmBuilder, i, q):
    if i % 2 == LF:  # left leg
        intersections = rbprmBuilder.getContactSurfacesAtConfig(q, "talos_lleg_rom")
    elif i % 2 == RF:  # right leg
        intersections = rbprmBuilder.getContactSurfacesAtConfig(q, "talos_rleg_rom")
    return intersections


# get all the contact surfaces (pts and normal)


def getAllSurfaces(afftool):
    l = afftool.getAffordancePoints("Support")
    return [(getSurfaceExtremumPoints(el), normal(el[0])) for el in l]


def getAllSurfacesDict(afftool):
    all_surfaces = getAllSurfaces(afftool)
    all_names = afftool.getAffRefObstacles("Support")  # id in names and surfaces match
    surfaces_dict = dict(
        zip(all_names, all_surfaces)
    )  # map surface names to surface points
    return surfaces_dict


def getSurfacesFromGuideContinuous(
    rbprmBuilder,
    ps,
    afftool,
    pathId,
    viewer=None,
    step=1.0,
    useIntersection=False,
):
    if viewer:
        from hpp.corbaserver.rbprm.tools.display_tools import (
            displaySurfaceFromPoints,
        )

    if pathId is None:
        pathId = ps.numberPaths() - 1
    pathLength = ps.pathLength(pathId)  # length of the path
    # get surface information
    surfaces_dict = getAllSurfacesDict(afftool)  # map surface names to surface points
    current_phase_end = step
    window_size = 0.5  # smaller step at which we check the colliding surfaces
    # Get surfaces for each phase along the path
    seqs = (
        []
    )  # list of list of surfaces : for each phase contain a list of surfaces. One phase is defined by moving of 'step' along the path
    i = 0
    t = 0.0
    end = False
    while not end:  # for all the path
        phase_contacts = []
        phase_contacts_names = []
        while (
            t < current_phase_end
        ):  # get the names of all the surfaces that the rom collide while moving from current_phase_end-step to current_phase_end
            q = ps.configAtParam(pathId, t)
            step_contacts = getContactsIntersections(rbprmBuilder, i, q)
            step_contacts_names = getContactsNames(
                rbprmBuilder, i, q
            )  # (FIX ME - RPBRM) step_contacts and step_contacts_names not in the same order
            # (FIX - PYTHON) Re-associate every intersections to their corresponding surfaces
            step_contacts_aux = []
            step_contacts_names_aux = []
            for index_intersection, intersection in enumerate(step_contacts):
                if len(intersection) > 0:  # Filter empty intersections
                    name = getNameSurfaceIntersected(
                        intersection, step_contacts_names, surfaces_dict
                    )
                    step_contacts_aux.append(intersection)
                    step_contacts_names_aux.append(name)
                    # print("  Name of surface ",step_contacts_names[index_intersection]," replaced by ",name)
            # Do not consider the duplicates yet
            phase_contacts += step_contacts_aux
            phase_contacts_names += step_contacts_names_aux
            t += window_size
            assert len(phase_contacts) == len(phase_contacts_names)
        # end current phase
        # Add surfaces or intersection to this phase.
        seq, seq_name = [], []
        if useIntersection:
            for index_phase, contact in enumerate(phase_contacts):
                # Add all intersections with ROMs => There can be two intersections on the same surface.
                # But in Sl1m paper : All surfaces must not intersect with one another (ERROR?)
                """
                if contact != [] and area(contact) > MAX_SURFACE:
                    seq.append(contact) # Add intersection
                    seq_name.append(phase_contacts_names[index_phase]) # Do something with it for solutions a,b or c
                """
                # Solutions : (a) Combine duplicates in one convex surface;
                #             (b) Keep the intersection with the biggest area ?
                #             (c) keep only the last intersection with this surface ?
                # Solution C implemented
                if area(contact) > MAX_SURFACE:
                    surface_name = phase_contacts_names[index_phase]
                    if surface_name in seq_name:
                        index_surface = seq_name.index(surface_name)
                        seq[index_surface] = contact  # Replace intersection
                    else:
                        seq.append(contact)  # Add intersection
                        seq_name.append(
                            phase_contacts_names[index_phase]
                        )  # Do something with it for solutions a,b or c
        else:
            for index_phase, name_contact in enumerate(phase_contacts_names):
                # Add surfaces intersecting with ROMs, only once per phase => No duplicate
                if not name_contact in seq_name:
                    seq.append(surfaces_dict[name_contact][0])  # Add entire surface
                    seq_name.append(name_contact)
        seq.sort()
        seqs.append(seq)
        # increase value for the next phase
        t = current_phase_end
        current_phase_end += step
        i += 1  # phase number
        if t == pathLength:
            current_phase_end = pathLength + 0.01
        if t > pathLength:
            end = True
    # end of all guide path
    # get rotation matrix of the root at each step
    seqs = listToArray(seqs)
    configs = []
    for t in arange(0, pathLength, step):
        configs.append(ps.configAtParam(pathId, t))
    R = getRotationsFromConfigs(configs)
    return R, seqs


def getSurfacesFromGuide(
    rbprmBuilder,
    ps,
    afftool,
    pathId,
    viewer=None,
    step=0.6,
    useIntersection=False,
):
    if viewer:
        from hpp.corbaserver.rbprm.tools.display_tools import (
            displaySurfaceFromPoints,
        )
    if pathId is None:
        pathId = ps.numberPaths() - 1
    pathLength = ps.pathLength(pathId)  # length of the path
    configs = []
    # get configuration along the path
    for s in arange(0, pathLength, step):
        configs.append(ps.configAtParam(pathId, s))

    # get surface information
    surface_dict = getAllSurfacesDict(afftool)

    # get surface candidate at each discretization step
    # suppose that we start with the left leg
    seqs = []
    for i, q in enumerate(configs):
        seq = []
        contacts = getContactsIntersections(rbprmBuilder, i, q)
        contact_names = getContactsNames(
            rbprmBuilder, i, q
        )  # Contacts and contact_names not in the same order
        # (FIX - PYTHON) Re-associate every intersections to their corresponding surfaces
        contacts_aux = []
        contact_names_aux = []
        for index_intersection, intersection in enumerate(contacts):
            if len(intersection) > 0:  # Filter empty intersections
                name = getNameSurfaceIntersected(
                    intersection, contact_names, surface_dict
                )
                contacts_aux.append(intersection)
                contact_names_aux.append(name)
                # print("  Name of surface ",contact_names[index_intersection]," replaced by ",name)
        contacts = contacts_aux
        contact_names = contact_names_aux
        # Add intersection or surface
        for j, contact in enumerate(contacts):
            if useIntersection:
                if area(contact) > MAX_SURFACE:
                    seq.append(contact)
            else:
                seq.append(surface_dict[contact_names[j]][0])
            if viewer:
                displaySurfaceFromPoints(viewer, contact, [0, 0, 1, 1])
        seq.sort()
        seqs.append(seq)

    # remove duplicates
    for i, seq in enumerate(seqs):
        seqs[i] = removeDuplicates(seq)

    seqs = listToArray(seqs)  # change the format from list to array
    R = getRotationsFromConfigs(configs)
    return R, seqs


# p : point tested
# p0, p1, p2 : points defining plane
def pointOnPlane(p, p0, p1, p2):
    n = np.cross(np.array(p1) - np.array(p0), np.array(p2) - np.array(p0))
    # n = n / np.linalg.norm(n)
    return np.dot(n, np.array(p) - np.array(p0)) < 1e-6


def sumAngles(p, surface):
    sum_angles = 0.0
    for i in range(len(surface) - 1):
        v0 = np.array(surface[i]) - np.array(p)
        v1 = np.array(surface[i + 1]) - np.array(p)
        v0, v1 = v0 / np.linalg.norm(v0), v1 / np.linalg.norm(v1)
        # print("  dot : ",np.dot(v0,v1))
        sum_angles += np.arccos(np.dot(v0, v1))
    # Last
    v0 = np.array(surface[-1]) - np.array(p)
    v1 = np.array(surface[0]) - np.array(p)
    v0, v1 = v0 / np.linalg.norm(v0), v1 / np.linalg.norm(v1)
    # print("  dot : ",np.dot(v0,v1))
    sum_angles += np.arccos(np.dot(v0, v1))
    return abs(sum_angles)


# Hypothesis :
# - intersection lies in on of the surface in nameSurfaces.
def getNameSurfaceIntersected(intersection, namesSurfaces, surfaces_dict):
    if len(intersection) > 0:
        for name in namesSurfaces:
            surface = surfaces_dict[name][0]
            # Check if a point of intersection is equal to one of this surface
            for p_intersection in intersection:
                for p_surface in surface:
                    if (
                        np.linalg.norm(np.array(p_intersection) - np.array(p_surface))
                        < 1e-6
                    ):
                        # Equal
                        return name
            # Check if the first point of intersection lies on plane of the surface
            p = intersection[0]
            if pointOnPlane(p, surface[0], surface[1], surface[2]):
                # Check if point inside or outside the surface
                # Sum of angles between p and all the vertex pairs.
                sum_angles = sumAngles(p, surface)
                if abs(sum_angles - np.pi * 2) < 1e-6:
                    return name
    return None
