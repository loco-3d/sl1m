import numpy as np


class Constraints:
    """
    Implementation of the class constraints, which implements the constraints used by the generic planner

    The problem is a LP with these variables for each phase:
    [ com_x, com_y, com_z_1, com_z_2, p_i_x, p_i_y, p_i_z, {a_i}                                ]
    [ 1    , 1    , 1      , 1      , 3                  , 0 if n_surfaces == 1, else n_surfaces]

    Under the following constraints :
    - fixed_foot_com: Ensures the COM is 'above' the fixed feet
    - foot_relative_distance: Ensures the moving feet is close enough to the other feet
    - surface: Each foot belongs to one surface
    - slack_positivity: The slack variables are positive
    - com_weighted_equality: Fix the horizontal position of the COm at the barycenter of the contact points (fixed feet)
    """

    def __init__(self, n_effectors):
        self.n_effectors = n_effectors
        self.default_n_variables = 7

        self.WEIGHTS = [1./float(n_effectors - 1)] * n_effectors
        
        self.com_xy = self._expression_matrix(2, 0)
        self.com_1 = self._expression_matrix(3, 0)
        self.com_2 = self._expression_matrix(3, 0)
        self.com_2[2, 2] = 0
        self.com_2[2, 3] = 1
        self.com_1_z = self._expression_matrix(1, 2)
        self.com_2_z = self._expression_matrix(1, 3)
        self.foot = self._expression_matrix(3, 4)
        self.foot_xy = self._expression_matrix(2, 4)

        self.M = 10.

    def _expression_matrix(self, size, j):
        """
        Generate a selection matrix for a given variable 
        @param size number of rows of the variable
        @param x position of the variable in the phase variables
        @return a (size, number of default variables (without slacks)) matrix with identity at column j
        """
        M = np.zeros((size, self.default_n_variables))
        M[:, j:j+size] = np.identity(size)
        return M

    def _fixed_foot_com_2_kinematic(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The COM 2 must belong to a reachable polytope above each end-effector of the phase
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
            if foot == phase.moving:
                l = k.shape[0]
                G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_2 - self.foot)
                h[i:i + l] = k
                i += l
            elif feet_phase[foot] != -1:
                l = k.shape[0]
                j_foot = js[feet_phase[foot]]
                G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_2)
                G[i:i + l, j_foot:j_foot + self.default_n_variables] = -K.dot(self.foot)
                h[i:i + l] = k
                i += l
            else:
                l = k.shape[0]
                G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_2)
                foot_pose = pb.p0[foot]
                h[i:i + l] = k + K.dot(foot_pose)
                i += l
        return i

    def _fixed_foot_com_1_kinematic(self, pb, phase,  G, h, i_start, js, feet_phase):
        """
        The COM 1 must belong to a reachable polytope above each end-effector of the previous phase
        For each effector id , K_id (c1- p_id(t-1)) <= k_id
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
            if foot == phase.moving:
                if feet_phase[foot] != -1:
                    l = k.shape[0]
                    j_foot = js[feet_phase[foot]]
                    G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_1)
                    G[i:i + l, j_foot:j_foot + self.default_n_variables] = -K.dot(self.foot)
                    h[i:i + l] = k
                    i += l
                else:
                    l = k.shape[0]
                    G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_1)
                    foot_pose = pb.p0[foot]
                    h[i:i + l] = k + K.dot(foot_pose)
                    i += l
            elif feet_phase[foot] != -1:
                l = k.shape[0]
                j_foot = js[feet_phase[foot]]
                G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_1)
                G[i:i + l, j_foot:j_foot + self.default_n_variables] = -K.dot(self.foot)
                h[i:i + l] = k
                i += l
            else:
                l = k.shape[0]
                G[i:i + l, j:j + self.default_n_variables] = K.dot(self.com_1)
                foot_pose = pb.p0[foot]
                h[i:i + l] = k + K.dot(foot_pose)
                i += l
        return i

    def fixed_foot_com(self, pb, phase,  G, h, i_start, js, phase_id, feet_phase):
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
        @param phase_id    Phase number
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return  i_start + the number of rows used by the constraint
        """
        i = i_start
        if phase_id != 0:
            i = self._fixed_foot_com_1_kinematic(pb, phase,  G, h, i_start, js, feet_phase)
        return self._fixed_foot_com_2_kinematic(pb, phase,  G, h, i, js, feet_phase)

    def foot_relative_distance(self, pb, phase, G, h, i_start, js, feet_phase):
        """
        The distance between the moving effector and the other ones is limited
        For i = moving_foot, For j !=i, Ki (pj - pi) <= ki
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
        constraints = phase.allRelativeK[phase.moving]
        for (foot, (K, k)) in constraints:
            l = k.shape[0]
            G[i:i + l, j:j + self.default_n_variables] = -K.dot(self.foot)
            if feet_phase[foot] != -1:
                j_foot = js[feet_phase[foot]]
                G[i:i + l, j_foot:j_foot + self.default_n_variables] = K.dot(self.foot)
                h[i:i + l] = k
            else:
                foot_pose = pb.p0[foot]
                h[i:i + l] = k - K.dot(foot_pose)
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
        n_surfaces = phase.n_surfaces
        if n_surfaces > 1:
            for s in range(n_surfaces):
                G[i+s, j + self.default_n_variables + s] = -1
            i += n_surfaces
        return i

    def surface_inequality(self, phase, G, h, i_start, j):
        """
        The moving foot must belong to one surface
        For each surface l: Sl pi - alphal <= sl
        @param phase       The phase specific data
        @param G           The inequality constraint matrix
        @param h           The inequality constraint vector
        @param i_start     Initial row to use
        @param j           Column corresponding to this phase variables
        @return i_start + the number of rows used by the constraint
        """
        i = i_start
        n_surfaces = phase.n_surfaces
        j_alpha = self.default_n_variables
        for S, s in phase.S:
            l = S.shape[0]
            G[i:i + l, j:j + self.default_n_variables] = S.dot(self.foot)
            h[i:i + l] = s
            if n_surfaces > 1:
                G[i:i + l, j + j_alpha] = -self.M * np.ones(l)
                j_alpha += 1
            i += l
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
        for foot in range(self.n_effectors):
            if foot == phase.moving:
                continue
            if feet_phase[foot] != -1:
                j_foot = js[feet_phase[foot]]
                C[i:i + 2, j_foot:j_foot + self.default_n_variables] = self.WEIGHTS[foot] * self.foot_xy
            else:
                foot_pose = pb.p0[foot]
                d[i:i + 2] -= self.WEIGHTS[foot] * foot_pose[:2]
        C[i:i + 2, j:j + self.default_n_variables] = -self.com_xy
        i += 2
        return i

    def _com_equality_init(self, pb, C, d, i_start):
        """
        The  initial com position is defined in the problem
        @param pb          The problem specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @return i_start + the number of rows used by the constraint
        """
        initial_position = pb.c0
        i = i_start
        C[i:i + 2, :self.default_n_variables] = self.com_xy
        d[i:i + 2] = initial_position[:2]
        i += 2
        return i

    def com(self, pb, phase, C, d, i_start, js, phase_id, feet_phase):
        """
        The com position is defined
        @param pb          The problem specific data
        @param phase       The phase specific data
        @param C           The equality constraint matrix
        @param d           The equality constraint vector
        @param i_start     Initial row to use
        @param js          List of column corresponding to the start of each phase
        @param phase_id    Phase number
        @param feet_phase List of the feet las moving phase, -1 if they haven't moved
        @return i_start + the number of rows used by the constraint
        """
        if phase_id != 0:
            return self._com_weighted_equality(pb, phase, C, d, i_start, js, feet_phase)
        else:
            return self._com_equality_init(pb, C, d, i_start)
