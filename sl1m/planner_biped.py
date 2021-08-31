import numpy as np
from sl1m.constants_and_tools import replace_surfaces_with_ineq_in_problem, normalize
from sl1m.constraints_biped import BipedConstraints


class BipedPlanner:
    """
    Implements an optimization problem to find the next surfaces to use for a biped robot

    The problem is a has these variables for each phase:
    [ p_i_x, p_i_y, p_i_z, com_z, {alpha_i, beta_i}                        ]
    [ 1    , 1    , 1    , 1    , 0 if n_surfaces == 1, else 2 * n_surfaces]

    Under the following constraints :
    - fixed_foot_com: Ensures the COM is 'above' the fixed feet
    - foot_relative_distance: Ensures the moving feet is close enough to the other feet
    - surface: Each foot belongs to one surface
    - slack_positivity: The slack variables are positive and beta is between -alpha and alpha
    """

    def __init__(self, mip=False):
        self.mip = mip
        self.default_n_variables = 4
        self.slack_scale = 10.

        self.n_slack_per_surface = 2
        self.n_ineq_per_surface = 2

        self.foot = self._expression_matrix(3, 0)
        self.com = self._expression_matrix(3, 0)
        self.com[2, 2] = 0
        self.com[2, 3] = 1
        self.com_xy = self.com[:2, :]
        self.foot_xy = self.foot[:2, :]

        self.cost_dict = {"final_com": self.end_com_cost,
                          "effector_positions": self.end_effectors_position_cost,
                          "coms": self.com_cost,
                          "posture": self.posture_cost,
                          "step_size": self.step_size_cost}

    def _expression_matrix(self, size, j):
        """
        Generate a selection matrix for a given variable
        @param size number of rows of the variable
        @param x position of the variable in the phase variables
        @return a(size, number of default variables(without slacks)) matrix with identity at column j
        """
        M = np.zeros((size, self.default_n_variables))
        M[:, j:j+size] = np.identity(size)
        return M

    def _slack_selection_vector(self):
        """
        Get the selection matrix of the slack variables
        @return the selection matrix
        """
        n_variables = self._total_n_variables()
        selection_vector = np.zeros(n_variables)
        i = 0
        for phase in self.pb.phaseData:
            phase_n_variables = self._phase_n_variables(phase)
            n_slacks = phase_n_variables - self.default_n_variables
            i_start = i + self.default_n_variables
            for i_slack in range(0, n_slacks, self.n_slack_per_surface):
                selection_vector[i_start + i_slack] = 1
            i += phase_n_variables
        return selection_vector

    def set_problem(self, pb):
        """
        Set the problem data
        @param pb new problem data
        """
        self.pb = pb
        self.alphas = self._slack_selection_vector()

    def _phase_n_variables(self, phase):
        """
        Counts the number of variables in a phase
        @param phase concerned phase
        @return number of variables in the phase
        """
        n_variables = self.default_n_variables
        n_surfaces = len(phase.S)
        if n_surfaces > 1:
            n_variables += n_surfaces * self.n_slack_per_surface
        return n_variables

    def _total_n_variables(self):
        """
        Counts the number of variables, inequalities constraints in the problem
        @return the number of variables of the problem
        """
        return sum([self._phase_n_variables(phase) for phase in self.pb.phaseData])

    def _phase_n_ineq(self, phase):
        """
        Counts the dimension of the inequalities in a phase
        - COM Kinematic constraints: summation over all effectors, times 2 because there are 2 height possible for the transition
        - Relative kinematic constraints between each effectors
        - Inequalities relative to each contact surface
        @param phase concerned phase
        """
        # COM kinematic constraints
        n_ineq = sum([k.shape[0] for (_, k) in phase.K])

        # Foot relative distance
        _, Ks = phase.allRelativeK[phase.moving][0]
        n_ineq += Ks[0].shape[0]

        # Surfaces
        n_ineq += sum([S[1].shape[0] - 1 for S in phase.S])

        # Slack positivity
        n_surfaces = len(phase.S)
        if n_surfaces > 1:
            n_ineq += n_surfaces * self.n_ineq_per_surface
        return n_ineq

    def _total_n_ineq(self):
        """
        Counts the number of inequality constraints
        @return the number of inequality constraints of the problem
        """
        return sum([self._phase_n_ineq(phase) for phase in self.pb.phaseData])

    def _phase_n_eq(self, phase):
        """
        Counts the dimension of the equality constraints of a phase
        @param phase concerned phase
        """
        n_eq = len(phase.S)
        if self.mip:
            for phase in self.pb.phaseData:
                n_surface = len(phase.S)
                if n_surface > 1:
                    n_eq += 1
        return n_eq

    def _total_n_eq(self):
        """
        Counts the number of equality constraints
        @return the number of equality constraints of the problem
        """
        return sum([self._phase_n_eq(phase) for phase in self.pb.phaseData])

    def convert_pb_to_LP(self, pb, convert_surfaces=True):
        """
        Compute the constraints:
        G x <= h
        C x = d

        @param pb problem data
        @return G, h, C, d
        """
        if convert_surfaces:
            replace_surfaces_with_ineq_in_problem(pb, True)

        self.set_problem(pb)

        n_variables = self._total_n_variables()
        n_ineq = self._total_n_ineq()
        n_eq = self._total_n_eq()

        A = np.zeros((n_ineq, n_variables))
        b = np.zeros(n_ineq)
        E = np.zeros((n_eq, n_variables))
        e = np.zeros(n_eq)

        i_start = 0
        i_start_eq = 0
        j = 0
        j_previous = 0

        constraints = BipedConstraints()
        for i, phase in enumerate(pb.phaseData):
            # inequality
            i_start = constraints.fixed_foot_com(pb, i, phase, A, b, j_previous, j, i_start)
            i_start = constraints.feet_relative_distance(pb, i, phase, A, b, j_previous, j, i_start)
            i_start = constraints.moving_foot_com(pb, i, phase, A, b, j_previous, j, i_start)
            i_start = constraints.surface_inequality(phase, A, b, j, i_start)
            i_start = constraints.slack_positivity(phase, A, b, j, i_start)

            # equality
            i_start_eq = constraints.surface_equality(phase, E, e, j, i_start_eq)

            if self.mip:
                i_start_eq = constraints.slack_equality(phase, E, e, i_start_eq, j)

            j_previous = j
            j += self._phase_n_variables(phase)

        A, b = normalize([A, b])
        E, e = normalize([E, e])
        return (A, b, E, e)

    def selected_surfaces(self, alphas):
        """
        : param: alphas the list of slack variables found by the planner
        Return the list of selected surfaces indices in the problem
        """
        indices = []
        for id, phase in enumerate(self.pb.phaseData):
            if phase.n_surfaces == 1:
                index = 0
            else:
                index = np.argmin(alphas[id])
            indices.append(index)
        return indices

    def get_alphas(self, result):
        """
        Retrieve the alphas from the result
        @param result vector of variables
        @return list of alpha values for each phase
        """
        alphas = []

        j = 0
        for phase in self.pb.phaseData:
            phase_alphas = []
            for j_slack in range(self.default_n_variables, self._phase_n_variables(phase), self.n_slack_per_surface):
                phase_alphas.append(result[j + j_slack])
            alphas.append(phase_alphas)

            j += self._phase_n_variables(phase)
        return alphas

    def get_result(self, result):
        """
        Retrieve the com and feet positions from the result
        @param result vector of variables
        @return com positions, moving foot positions and all feet positions
        """
        coms = []
        moving_foot_pos = []
        all_feet_pos = [[self.pb.p0[0]], [self.pb.p0[1]]]

        j = 0
        for i, phase in enumerate(self.pb.phaseData):
            moving_foot = phase.moving
            fixed_foot = (moving_foot + 1) % 2

            # CoM
            coms.append(self.com.dot(result[j:j + self.default_n_variables]))

            # Moving foot
            moving_foot_pos.append(self.foot.dot(result[j:j + self.default_n_variables]))

            # All feet
            all_feet_pos[moving_foot].append(moving_foot_pos[i])
            all_feet_pos[fixed_foot].append(all_feet_pos[fixed_foot][-1])

            j += self._phase_n_variables(phase)

        all_feet_pos[0].pop(0)
        all_feet_pos[1].pop(0)

        return coms, moving_foot_pos, all_feet_pos

    def com_cost(self, coms):
        """
        Compute a cost to keep the com close to a target com at each phase
        @param coms list of target positions for the com
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            A = np.zeros((2, n_variables))
            A[:, j: j + self.default_n_variables] = self.com
            b = coms[id][:2]

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self._phase_n_variables(phase)

        return P, q

    def end_com_cost(self, com):
        """
        Compute a cost to keep the final CoM position close to a target one
        @param com Target final com position
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        n_phases = self.pb["n_phases"]
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            if id == n_phases - 1:
                A = np.zeros((2, n_variables))
                A[:, j: j + self.default_n_variables] = self.com_xy
                b = com[:2]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
            j += self._phase_n_variables(phase)

        return P, q

    def end_effectors_position_cost(self, effector_positions):
        """
        Compute a cost to keep the final end effectors positions closed to target ones
        @param effector_positions list of each effector's final position
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        n_phases = self.pb["n_phases"]
        j = 0.
        for id, phase in enumerate(self.pb.phaseData):
            if id >= n_phases - 2:
                A = np.zeros((3, n_variables))
                A[:, j: j + self.default_n_variables] = self.foot
                b = effector_positions[phase.moving][:3]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
            j += self._phase_n_variables(phase)

        return P, q

    def posture_cost(self):
        """
        Compute a cost to keep the feet relative positions as close as possible from the initial ones
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        relative_positions = [self.pb.p0[1] - self.pb.p0[0]]

        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j_previous = 0
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            A = np.zeros((3, n_variables))
            b = np.zeros(3)
            if id == 0:
                if phase.moving == 0:
                    A[:, j:j + self.default_n_variables] = -self.foot
                    b = -self.pb.p0[1]
                else:
                    A[:, j:j + self.default_n_variables] = self.foot
                    b = self.pb.p0[0]
            else:
                if phase.moving == 0:
                    A[:, j:j + self.default_n_variables] = -self.foot
                    A[:, j_previous:j_previous + self.default_n_variables] = self.foot
                else:
                    A[:, j:j + self.default_n_variables] = self.foot
                    A[:, j_previous:j_previous + self.default_n_variables] = -self.foot
            b += relative_positions

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j_previous = j
            j += self._phase_n_variables(phase)

        return P, q

    def step_size_cost(self, step_size):
        """
        Compute a cost to keep the step sizes as close as possible of a target step size
        @param step_size desired size of the steps
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j_previous = 0
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            foot = phase.moving
            A = np.zeros((2, n_variables))
            b = step_size
            A[:, j:j + self.default_n_variables] = self.foot_xy
            if id == 0:
                b += self.pb.p0[foot][:2]
            else:
                A[:, j_previous:j_previous + self.default_n_variables] = -self.foot_xy

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j_previous = j
            j += self._phase_n_variables(phase)

        return P, q

    def compute_costs(self, costs):
        """
        This function computes the P and q cost matrix and vector given a cost dictionary

        @param costs the cost dictionary. The keys should match keys from self.cost_dict
        @return P, q the cost matrix and vector of the QP problem
        """
        n_variables = self._total_n_variables()

        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        if costs == {}:
            P += np.identity(n_variables)
        else:
            for key, val in [(k, v) for k, v in costs.items() if k in self.cost_dict]:
                if key in self.cost_dict.keys():
                    P_, q_ = self.cost_dict[key](*val[1:])
                    P += P_ * val[0]
                    q += q_ * val[0]
                else:
                    print("Unknown cost to add to the biped problem")

        return P, q
