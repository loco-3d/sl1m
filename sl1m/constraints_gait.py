import numpy as np


class Constraints:
    """
    Implementation of the class constraints, which implements the constraints used by the generic planner

    The problem is a LP with these variables for each phase:
    [ com_x, com_y, com_z_1, com_z_2, {p_i_x, p_i_y, p_i_z}, {a_ij}                                                  ]
    [ 1    , 1    , 1      , 1      , 3 * phase_n_moving   , (0 if n_surfaces == 1, else n_surfaces) * phase_n_moving]

    Under the following constraints :
    - fixed_foot_com: Ensures the COM is 'above' the fixed feet
    - foot_relative_distance: Ensures the moving feet are close enough to the other feet
    - surface: Each foot belongs to one surface
    - slack_positivity: The slack variables are positive
    - com_weighted_equality: Fix the horizontal position of the COm at the barycenter of the contact points (fixed feet)
    """

    def __init__(self, n_effectors, com=True):
        self.n_effectors = n_effectors

        self.default_n_variables = 4 * int(com)

        self.M = 100

    def _default_n_variables(self, phase):
        """
        @param phase phase concerned
        @return the number of non slack variables in phase
        """
        return self.default_n_variables + 3 * len(phase.moving)

    def _expression_matrix(self, size, _default_n_variables, j):
        """
        Generate a selection matrix for a given variable 
        @param size number of rows of the variable
        @param _default_n_variables number of non slack variables in the phase
        @param x position of the variable in the phase variables
        @return a (size, number of default variables (without slacks)) matrix with identity at column j
        """
        M = np.zeros((size, _default_n_variables))
        M[:, j:j + size] = np.identity(size)
        return M

    def com_xy(self, phase):
        """
        Generate a selection matrix for the com x and y components
        @param phase phase data
        @return a (2, phase number of variables (without slacks)) matrix 
        """
        return self._expression_matrix(2, self._default_n_variables(phase), 0)

    def com_1(self, phase):
        """
        Generate a selection matrix for the com_1
        @param phase phase data
        @return a (3, phase number of variables (without slacks)) matrix 
        """
        return self._expression_matrix(3, self._default_n_variables(phase), 0)

    def com_2(self, phase):
        """
        Generate a selection matrix for the com_2
        @param phase phase data
        @return a (3, phase number of variables (without slacks)) matrix 
        """
        M = self._expression_matrix(3, self._default_n_variables(phase), 0)
        M[2, 2] = 0
        M[2, 3] = 1
        return M

    def com_1_z(self, phase):
        """
        Generate a selection matrix for the com_1 z component
        @param phase phase data
        @return a (1, phase number of variables (without slacks)) matrix 
        """
        return self._expression_matrix(1, self._default_n_variables(phase), 2)

    def com_2_z(self, phase):
        """
        Generate a selection matrix for the com_2 z component
        @param phase phase data
        @return a (1, phase number of variables (without slacks)) matrix 
        """
        return self._expression_matrix(3, self._default_n_variables(phase), 3)

    def foot(self, phase, foot=None, id=None):
        """
        Generate a selection matrix for a given foot 
        @param phase phase data
        @param foot foot to select
        @param id id of the foot in the moving feet list
        @return a (3, number of variables (without slacks)) selection matrix
        """
        if foot is not None:
            id = np.argmax(phase.moving == foot)
        elif id is None:
            print("Error in foot selection matrix: you must specify either foot or id")
        j = self.default_n_variables + 3 * id
        return self._expression_matrix(3, self._default_n_variables(phase), j)

    def foot_xy(self, phase, foot=None, id=None):
        """
        Generate a selection matrix for a given foot x and y coordinates
        @param phase phase data
        @param foot foot to select
        @param id id of the foot in the moving feet list
        @return a (3, number of variables (without slacks)) selection matrix
        """
        if foot is not None:
            id = np.argmax(phase.moving == foot)
        elif id is None:
            print("Error in foot selection matrix: you must specify either foot or id")
        j = self.default_n_variables + 3 * id
        return self._expression_matrix(2, self._default_n_variables(phase), j)

    def _fixed_foot_com_2_kinematic(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The COM 2 must belong to a reachable polytope above each stance foot of the phase
        For each effector id , K_id (c2 - p_id) <= k_id
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return  i_start + the number of rows used by the constraint
        """
        i = i_start
        j = js[-1]

        for foot, (K, k) in enumerate(phase.K):
            if foot in phase.stance:
                l = k.shape[0]
                G[i:i + l, j:j + self._default_n_variables(phase)] = K.dot(self.com_2(phase))
                h[i:i + l] = k
                if foot in phase.moving:
                    G[i:i + l, j:j + self._default_n_variables(phase)] -= K.dot(self.foot(phase, foot))
                elif feet_phase[foot] != -1:
                    j_f = js[feet_phase[foot]]
                    phase_f = pb.phaseData[feet_phase[foot]]
                    G[i:i + l, j_f:j_f + self._default_n_variables(phase_f)] = -K.dot(self.foot(phase_f, foot))
                else:
                    foot_pose = pb.p0[foot]
                    h[i:i + l] += K.dot(foot_pose)
                i += l
        return i

    def _fixed_foot_com_1_kinematic(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The COM 1 must belong to a reachable polytope above each stance foot of the previous phase
        For each effector id in stance phase, K_id (c1- p_id(t-1)) <= k_id
        This constraint should not be added at the first phase.
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return  i_start + the number of rows used by the constraint
        """
        i = i_start
        j = js[-1]
        for foot, (K, k) in enumerate(phase.K):
            if phase.id == 0:
                l = k.shape[0]
                foot_pose = pb.p0[foot]
                G[i:i + l, j:j + self._default_n_variables(phase)] = K.dot(self.com_1(phase))
                h[i:i + l] = k + K.dot(foot_pose)
                i += l
            elif foot in pb.phaseData[phase.id - 1].stance:
                l = k.shape[0]
                G[i:i + l, j:j + self._default_n_variables(phase)] = K.dot(self.com_1(phase))
                h[i:i + l] = k
                if feet_phase[foot] != -1:
                    j_f = js[feet_phase[foot]]
                    phase_f = pb.phaseData[feet_phase[foot]]
                    G[i:i + l, j_f:j_f + self._default_n_variables(phase_f)] -= K.dot(self.foot(phase_f, foot))
                else:
                    foot_pose = pb.p0[foot]
                    h[i:i + l] += K.dot(foot_pose)
                i += l
        return i

    def fixed_foot_com(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The COM must belong to a reachable polytope above each end-effector
        For each effector id , K_id (c(t) - p_id(t-1)) <= k_id
        This constraint should not be added at the first phase.
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase  List of the feet last moving phase, -1 if they haven't moved
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        if phase.id != 0:
            i = self._fixed_foot_com_1_kinematic(pb, phase, G, h, i, js, feet_phase)
        return self._fixed_foot_com_2_kinematic(pb, phase, G, h, i, js, feet_phase)

    def foot_relative_distance(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The distance between the moving effector and the other stance feet is limited
        For i in moving_foot, For j !=i, Ki (pj - pi) <= ki
        The distane between the moving effector and the previous phase stance feet (or initial 
        contacts for the fist phase) is also limited
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        j = js[-1]
        for foot in phase.moving:
            constraints = phase.allRelativeK[foot]
            for (other, (K, k)) in constraints:
                if other in phase.stance:
                    l = k.shape[0]
                    G[i:i + l, j:j + self._default_n_variables(phase)] = -K.dot(self.foot(phase, foot))
                    h[i:i + l] = k
                    if other in phase.moving:
                        G[i:i + l, j:j + self._default_n_variables(phase)] += K.dot(self.foot(phase, other))
                    elif feet_phase[other] != -1:
                        j_f = js[feet_phase[other]]
                        phase_f = pb.phaseData[feet_phase[other]]
                        G[i:i + l, j_f:j_f + self._default_n_variables(phase_f)] = K.dot(self.foot(phase_f, other))
                    else:
                        foot_pose = pb.p0[other]
                        h[i:i + l] -= K.dot(foot_pose)
                    i += l
                elif phase.id == 0:
                    l = k.shape[0]
                    G[i:i + l, j:j + self._default_n_variables(phase)] = -K.dot(self.foot(phase, foot))
                    h[i:i + l] = k
                    h[i:i + l] -= K.dot(pb.p0[other])
                    i += l
                elif other in pb.phaseData[phase.id - 1].stance:
                    l = k.shape[0]
                    G[i:i + l, j:j + self._default_n_variables(phase)] = -K.dot(self.foot(phase, foot))
                    h[i:i + l] = k
                    if feet_phase[other] != -1:
                        j_f = js[feet_phase[other]]
                        phase_f = pb.phaseData[feet_phase[other]]
                        G[i:i + l, j_f:j_f + self._default_n_variables(phase_f)] = K.dot(self.foot(phase_f, other))
                    else:
                        foot_pose = pb.p0[other]
                        h[i:i + l] -= K.dot(foot_pose)
                    i += l
        return i

    def slack_positivity(self, phase, G, h, i_start, j):
        """
        The slack variables (alpha) should be positive
        Sl for each surface s, -alpha_s <= 0
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param j           Column corresponding to this phase variables
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        j_alpha = j + self._default_n_variables(phase)
        for n_surface in phase.n_surfaces:
            if n_surface > 1:
                G[i:i + n_surface, j_alpha:j_alpha + n_surface] = -np.identity(n_surface)
                j_alpha += n_surface
                i += n_surface
        return i

    def surface_inequality(self, phase, G, h, i_start, j):
        """
        Each moving foot must belong to one surface
        For each surface l: Sl pi - alphal <= sl
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param j           Column corresponding to this phase variables
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        j_alpha = self._default_n_variables(phase)
        for id, surfaces in enumerate(phase.S):
            for S, s in surfaces:
                l = S.shape[0]
                G[i:i + l, j:j + self._default_n_variables(phase)] = S.dot(self.foot(phase, id=id))
                h[i:i + l] = s
                if phase.n_surfaces[id] > 1:
                    G[i:i + l, j + j_alpha] = -self.M * np.ones(l)
                    j_alpha += 1
                i += l
        return i

    def slack_equality(self, phase, C, d, i_start, j):
        """
        The slack variables (alpha) sum should be equal to the number of surfaces -1 
        Sl for each moving foot, sum(alpha_s) = n_surfaces - 1
        @param phase       The phase specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @param j           Column corresponding to this phase variables
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        j_alpha = j + self._default_n_variables(phase)
        for n_surface in phase.n_surfaces:
            if n_surface > 1:
                C[i, j_alpha:j_alpha + n_surface] = np.ones(n_surface)
                d[i] = n_surface - 1
                j_alpha += n_surface
                i += 1
        return i

    def _com_weighted_equality(self, pb, phase, C, d, i_start, js, feet_phase):
        """
        The 2D position of the com should be the baricenter of the fixed feet locations
        0 =  sum(fixed_foot i) WEIGHT * p_i_x_y - c_x_y
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        j = js[-1]

        weight = 1./len(phase.stance)

        for foot in range(self.n_effectors):
            if foot in phase.stance:
                if feet_phase[foot] != -1:
                    j_f = js[feet_phase[foot]]
                    phase_f = pb.phaseData[feet_phase[foot]]
                    C[i:i + 2, j_f:j_f + self._default_n_variables(phase_f)] = weight * self.foot_xy(phase_f, foot)
                else:
                    foot_pose = pb.p0[foot]
                    d[i:i + 2] -= weight * foot_pose[:2]
        C[i:i + 2, j:j + self._default_n_variables(phase)] = -self.com_xy(phase)
        i += 2
        return i

    def _com_equality_init(self, pb, phase, C, d, i_start):
        """
        The initial com position is defined in the problem
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        C[i:i + 2, :self._default_n_variables(phase)] = self.com_xy(phase)
        d[i:i + 2] = pb.c0[:2]
        i += 2
        return i

    def com(self, pb, phase, C, d, i_start, js, feet_phase):
        """
        The com position is defined by the initial data or as the barycenter of fixed feet positions
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return i_start + the number of rows used by the constraint
        """
        if phase.id == 0 and pb.c0 != None:
            return self._com_equality_init(pb, phase, C, d, i_start)
        else:
            return self._com_weighted_equality(pb, phase, C, d, i_start, js, feet_phase)
