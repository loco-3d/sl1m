import numpy as np
from numpy import (
    array,
    zeros,
    ones,
    vstack,
    hstack,
    identity,
    cross,
    concatenate,
)
from numpy.linalg import norm
from scipy.spatial import ConvexHull

from sl1m.tools.obj_to_constraints import (
    load_obj,
    as_inequalities,
    rotate_inequalities,
    inequalities_to_Inequalities_object,
)

# --------------------------------- CONSTANTS ---------------------------------------------------------------

EPSILON = 0.000001
EPSILON_EQ = 0.0

IDENTITY = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
GRAVITY = array([0.0, 0.0, -9.81])
GRAVITY_6 = array([0.0, 0.0, -9.81, 0.0, 0.0, 0.0])

X = array([1.0, 0.0, 0.0])
Y = array([0.0, 1.0, 0.0])
Z = array([0.0, 0.0, 1.0])

# --------------------------------- METHODS ---------------------------------------------------------------


def normalize(Ab):
    A = Ab[0]
    b = Ab[1]
    A_normalized = zeros(A.shape)
    b_normalized = zeros(b.shape)
    for i in range(A.shape[0]):
        n = norm(A[i, :])
        if n <= EPSILON:
            n = 1.0
        A_normalized[i, :] = A[i, :] / n
        b_normalized[i] = b[i] / n
    return A_normalized, b_normalized


def convert_surface_to_inequality(s, eq_as_ineq):
    # TODO does normal orientation matter ? It will for collisions
    n = cross(s[:, 1] - s[:, 0], s[:, 2] - s[:, 0])
    if n[2] <= 0.0:
        for i in range(3):
            n[i] = -n[i]
    norm_n = norm(n)
    if norm_n > 1e-6:
        n /= norm(n)
    else:
        # FIXME: to something better here
        n = array([0, 0, 1])
    return surface_points_to_inequalities(s, n, eq_as_ineq)


def normal_from_ineq(s_ineq):
    n = s_ineq[0][-1]
    if n[2] < 0:
        n = -n
    return n


def replace_surfaces_with_ineq_in_phaseData(phase, eq_as_ineq):
    phase.S = [convert_surface_to_inequality(S, eq_as_ineq) for S in phase.S]


def replace_surfaces_with_ineq_in_problem(pb, eq_as_ineq=False):
    [
        replace_surfaces_with_ineq_in_phaseData(phase, eq_as_ineq)
        for phase in pb.phaseData
    ]


def ineqQHull(hull):
    A = hull.equations[:, :-1]
    b = -hull.equations[:, -1]
    return A, b


def default_transform_from_pos_normal(pos, normal, transform=IDENTITY):
    """
    transform matrix is a rotation matrix from the root trajectory
    used to align the foot yaw (z-axis) orientation
    surface normal is used to align the foot roll and pitch (x- and y- axes) orientation
    """
    f = Z
    t = array(normal)
    t = t / norm(t)
    v = np.cross(f, t)
    c = np.dot(f, t)
    if abs(c) > 0.99:
        rot = identity(3)
    else:
        u = v / norm(v)
        h = (1.0 - c) / (1.0 - c**2)
        s = 1.0 - c**2
        assert s > 0
        s = np.sqrt(s)
        ux, uy, uz = u
        rot = array(
            [
                [c + h * ux**2, h * ux * uy - uz * s, h * ux * uz + uy * s],
                [h * ux * uy + uz * s, c + h * ux**2, h * uy * uz - ux * s],
                [h * ux * uz - uy * s, h * uy * uz + ux * s, c + h * uz**2],
            ]
        )

    rot = np.dot(transform, rot)
    return vstack([hstack([rot, pos.reshape((-1, 1))]), [0.0, 0.0, 0.0, 1.0]])


def surface_points_to_inequalities(S, normal, eq_as_ineq):
    # Homogenous tranformation matrix from world frame to surface frame.
    # Normal vector corresponds to z-axis.
    # 1st vertex is the Origin.
    n = array(normal)
    tr = default_transform_from_pos_normal(S.T[0], n)

    # Projection of vertices on the surface local frame --> z = 0. (2D)
    # P_world = R_ @ P_local + T_
    # --> P_local = R_-1 @ (P_world - T_)
    trpts = [np.dot(tr[:3, :3].T , (pos - tr[:-1,-1]))[:2] for pos in S.T]

    hull = ConvexHull(array(trpts))
    A, b = ineqQHull(hull)
    # Increase inequality matrix with z = 0 column
    A = hstack([A, zeros((A.shape[0], 1))])
    # In local frame : A_l @ P_local <= b_l
    # --> A_l @ R-1 @ P_world <= b_l + A_l @ R_-1 @ T_
    A_world = np.dot(A,tr[:3,:3].T)
    b_world = b + np.dot(np.dot(A,tr[:3,:3].T), tr[:-1,-1])
    ine = inequalities_to_Inequalities_object(A_world, b_world)

    d = array([n.dot(S[:, 0])])
    if eq_as_ineq:
        A = vstack([ine.A, n, -n])
        b = concatenate([ine.b, d + EPSILON_EQ, -d + EPSILON_EQ]).reshape((-1,))
    else:
        A = vstack([ine.A, n, -n])
        b = concatenate([ine.b, d, -d]).reshape((-1,))
        A = vstack([ine.A, n])
        b = concatenate([ine.b, d]).reshape((-1,))

    return A, b


def timMs(t1, t2):
    """
    Expresses the duration between t1 and t2 in milliseconds
    """
    return (t2 - t1) * 1000.0
