import numpy as np

# Implementation of the class constraints, which implements the constraints used by the generic planner
#
# The problem is a has these variables for each phase:
# [ p_i_x, p_i_y, p_i_z, com_z, {alpha_0_i, alpha_1_i}                   ]
# [ 1    , 1    , 1    , 1    , 0 if n_surfaces == 1, else 2 * n_surfaces]

# alpha_0_i, slack variable for inequality, alpha_1_i slack for equality constraint
# Extra constraints are alpha_1_i - alpha_0_i <= 0 ; -alpha_1_i - alpha_0_i <=  0;  alpha_0_i >= 0 is implied by the first 2
#
# Under the following constraints :
# - fixed_foot_com: Ensures the COM is 'above' the fixed feet
# - foot_relative_distance: Ensures the moving feet is close enough to the other feet
# - surface: Each foot belongs to one surface
# - slack_positivity: Inequality constraints on the slack variables


class BipedConstraints:
    def __init__(self):
        self.default_n_variables = 4
        self.slack_scale = 10.

        self.n_slack_per_surface = 2
        self.n_ineq_per_surface = 2

        self.foot = self.__expression_matrix(3, 0)
        self.foot_xy = self.foot[:2, :]
        self.foot_z = self.foot[2:, :]
        self.com_z = self.__expression_matrix(1, 3)
        self.com = self.__expression_matrix(3, 0)
        self.com[2, :] = self.com_z

    def __fixed(self, moving_foot):
        """
        Get the index of the fixed foot given the moving foot index
        @param moving_foot moving foot index
        @return fixed foot index
        """
        return (moving_foot + 1) % 2

    def __expression_matrix(self, size, j):
        """
        Generate a selection matrix for a given variable
        @param size number of rows of the variable
        @param x position of the variable in the phase variables
        @return a (size, number of default variables (without slacks)) matrix with identity at column j
        """
        M = np.zeros((size, self.default_n_variables))
        M[:, j:j + size] = np.identity(size)
        return M

    def fixed_foot_com(self, pb, id, phase, A, b, j_previous, j, i):
        """
        The COM must be in the reachable polytope of the fixed foot (whose position is stored in "p0")
        K_fixed (cz - p_fixed_z) <= k_fixed
        :param: id          Phase id
        :param: pb          The problem specific data
        :param: phase       The phase specific data
        :param: G           The inequality constraint matrix
        :param: h           The inequality constraint vector
        :param: i_start     Initial row to use
        :param: js          List of j_startumn corresponding to the start of each phase
        return  i_start + the number of rows used by the constraint
        """
        fixed_foot = self.__fixed(phase.moving)
        K, k = phase.K[fixed_foot]
        l = K.shape[0]

        if id == 0:
            if pb.p0 is not None:
                fixed_foot_z = pb.p0[fixed_foot][-1:]
                A[i:i + l, j:j + self.default_n_variables] = K[:, -1:].dot(self.com_z)
                b[i:i + l] = k + K[:, -1:].dot(fixed_foot_z)
        else:
            A[i:i + l, j:j + self.default_n_variables] = K[:, -1:].dot(self.com_z)
            A[i:i + l, j_previous:j_previous + self.default_n_variables] = - \
                K[:, -1:].dot(self.foot_z)
            b[i:i + l] = k

        return i + l

    def moving_foot_com(self, pb, id, phase, A, b, j_previous, j, i):
        """
        The COM must be in the reachable polytope of the moving foot (whose position is stored in "p0")
        K_moving_z * cz - K_moving * (p_moving - p_fixed) <= k_moving
        This contraint is only valid at the first phase
        :param: id          Phase id
        :param: pb          The problem specific data
        :param: phase       The phase specific data
        :param: G           The inequality constraint matrix
        :param: h           The inequality constraint vector
        :param: i_start     Initial row to use
        :param: js          List of j_startumn corresponding to the start of each phase
        return  i_start + the number of rows used by the constraint
        """
        moving_foot = phase.moving
        fixed_foot = self.__fixed(moving_foot)
        K, k = phase.K[moving_foot]
        l = K.shape[0]

        if id == 0:
            if pb.p0 is not None:
                com_xy = pb.p0[fixed_foot][:2]
                A[i:i + l, j:j + self.default_n_variables] = K[:, -
                                                               1:].dot(self.com_z) - K.dot(self.foot)
                b[i:i + l] = k - K[:, :2].dot(com_xy[:2])
        else:
            A[i:i+l, j:j + self.default_n_variables] = K[:, -1:].dot(self.com_z) - K.dot(self.foot)
            A[i:i+l, j_previous:j_previous + self.default_n_variables] = K[:, :2].dot(self.foot_xy)
            b[i:i+l] = k

        return i + l

    def feet_relative_distance(self, pb, id, phase, A, b, j_previous, j, i):
        """
        The distance between the fixed and moving effectors is limited
        K_moving (pj - pi) <= k_moving
        This contraint is only valid at the first phase
        :param: id          Phase id
        :param: pb          The problem specific data
        :param: phase       The phase specific data
        :param: G           The inequality constraint matrix
        :param: h           The inequality constraint vector
        :param: i_start     Initial row to use
        :param: js          List of j_startumn corresponding to the start of each phase
        return  i_start + the number of rows used by the constraint
        """
        moving_foot = phase.moving
        fixed_foot, (K, k) = phase.allRelativeK[moving_foot][0]
        l = K.shape[0]

        if id == 0:
            if pb.p0 is not None:
                fixed_foot_position = pb.p0[fixed_foot]
                A[i:i + l, j:j + self.default_n_variables] = K.dot(self.foot)
                b[i:i + l] = k + K.dot(fixed_foot_position)
                return i + l
        else:
            A[i:i + l, j:j + self.default_n_variables] = K.dot(self.foot)
            A[i:i + l, j_previous:j_previous + self.default_n_variables] = -K.dot(self.foot)
            b[i:i + l] = k

        return i + l


    def surface_inequality(self, phase, A, b, j, i_start):
        """
        The moving foot must belong to one surface
        For each surface l: Sl pi - alphal <= sl
        :param: phase       The phase specific data
        :param: C           The equality constraint matrix
        :param: d           The equality constraint vector
        :param: i     Initial row to use
        :param: j           Column corresponding to this phase variables
        return i + the number of rows used by the constraint
        """
        n_surfaces = len(phase.S)
        j_slack = self.default_n_variables
        i = i_start
        for (S, s) in phase.S:
            l = S.shape[0]-1
            A[i:i + l, j:j + self.default_n_variables] = S[:-1, :].dot(self.foot)
            b[i:i + l] = s[:-1]
            if n_surfaces > 1:
                A[i:i + l, j + j_slack] = - np.ones(l) * self.slack_scale
                j_slack += self.n_slack_per_surface
            i += l
        return i

    def surface_equality(self, phase, E, e, j, i_start):
        """
        The moving foot must belong to one surface
        For each surface l: di * p_i = e_i + beta_i
        :param: phase       The phase specific data
        :param: C           The equality constraint matrix
        :param: d           The equality constraint vector
        :param: i_start     Initial row to use
        :param: j           Column corresponding to this phase variables
        return i_start + the number of rows used by the constraint
        """
        n_surfaces = len(phase.S)
        i = i_start
        if n_surfaces == 1:
            E[i, j:j + self.default_n_variables] = phase.S[0][0][-1].dot(self.foot)
            e[i] = phase.S[0][1][-1]
            return i + 1
        else:
            j_slack = self.default_n_variables + 1
            for (S, s) in phase.S:
                E[i, j:j + self.default_n_variables] = S[-1, :].dot(self.foot)
                E[i, j + j_slack] = -1 * self.slack_scale
                e[i] = s[-1]
                j_slack += self.n_slack_per_surface
                i += 1
            return i

    
    def slack_positivity(self, phase, A, b, j_start, i_start):
        """
        The slack variables (alpha) should be positive
        for each surface s:
        -alpha_s + beta_s <= 0
        -alpha_s - beta_s <= 0
        :param: phase       The phase specific data
        :param: G           The inequality constraint matrix
        :param: h           The inequality constraint vector
        :param: i_start     Initial row to use
        :param: j           Column corresponding to this phase variables
        return i_start + the number of rows used by the constraint
        """
        n_surfaces = len(phase.S)
        i = i_start
        j = j_start
        if n_surfaces > 1:
            j += self.default_n_variables
            for j_slack in range(0, n_surfaces * self.n_slack_per_surface, self.n_slack_per_surface):
                A[i, j + j_slack:j + j_slack + 2] = [-1,  1]
                A[i + 1, j + j_slack:j + j_slack + 2] = [-1, -1]
                i += self.n_ineq_per_surface
        return i
