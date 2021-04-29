import numpy as np
from sl1m.constants_and_tools import *
from sl1m.constraints import Constraints

# Implements an optimization problem to find the next surfaces to use
#
# The problem is a LP with these variables for each phase:
# [ com_x, com_y, com_z_1, com_z_2, p_i_x, p_i_y, p_i_z, {a_i}                                ]
# [ 1    , 1    , 1      , 1      , 1    , 1    , 1    , 0 if n_surfaces == 1, else n_surfaces]
#
# Under the following constraints :
# - fixed_foot_com: Ensures the COM is 'above' the fixed feet
# - foot_relative_distance: Ensures the moving feet is close enough to the other feet
# - surface: Each foot belongs to one surface
# - slack_positivity: The slack variables are positive
# - com_weighted_equality: Fix the horizontal position of the CoM at the barycenter of the contact points (fixed feet)


class Planner:
    def __init__(self):
        self.default_n_variables = 7
        self.default_n_equality_constraints = 2

        self.com_xy = self.__expression_matrix(2, 0)
        self.com_2 = self.__expression_matrix(3, 0)
        self.com_2[2, 2] = 0
        self.com_2[2, 3] = 1
        self.foot = self.__expression_matrix(3, 4)
        self.foot_xy = self.__expression_matrix(2, 4)

        self.cost_dict = {"final_com": self.end_com_cost,
                          "effector_positions": self.end_effectors_position_cost,
                          "coms": self.com_cost,
                          "posture": self.posture_cost,
                          "step_size": self.step_size_cost}

    def __expression_matrix(self, size, j):
        """
        Generate a selection matrix for a given variable 
        @param size number of rows of the variable
        @param x position of the variable in the phase variables
        @return a (size, number of default variables (without slacks)) matrix with identity at column j
        """
        M = np.zeros((size, self.default_n_variables))
        M[:, j:j+size] = np.identity(size)
        return M

    def __slack_selection_vector(self):
        """
        Get the selection vector of the slack variables
        @return the selection vector
        """
        n_variables = self.__total_n_variables()
        selection_vector = np.zeros(n_variables)
        i = 0
        for phase in self.pb.phaseData:
            n_variables_phase = self.__phase_n_variables(phase)
            n_slacks = n_variables_phase - self.default_n_variables
            i_start = i + self.default_n_variables
            selection_vector[i_start:i_start + n_slacks] = 1
            i += n_variables_phase
        return selection_vector

    def set_problem(self, pb):
        """
        Set the problem data
        @param pb new problem data
        """
        self.pb = pb
        self.n_effectors = pb.n_effectors
        self.alphas = self.__slack_selection_vector()

    def __feet_last_moving_phase(self, phase_id):
        """
        Ids of the last phase each foot moved at a given phase
        @param phase_id id of the phase
        @return a list of n_effetors terms corresponding to the last phase the feet moved, -1 if it hasn't moved yet
        """

        feet_last_moving_phase = [-1 for _ in range(self.n_effectors)]
        for id, phase in enumerate(self.pb.phaseData):
            if id >= phase_id:
                break
            moving_foot = phase.moving
            feet_last_moving_phase[moving_foot] = id
        return feet_last_moving_phase

    def __phase_n_variables(self, phase):
        """
        Counts the number of variables in a phase
        @param phase concerned phase
        @return number of variables in the phase
        """
        n_variables = self.default_n_variables
        n_surfaces = phase.n_surfaces
        if n_surfaces > 1:
            n_variables += n_surfaces
        return n_variables

    def __total_n_variables(self):
        """
        Counts the number of variables, inequalities constraints in the problem
        @return the number of variables of the problem
        """
        return sum([self.__phase_n_variables(phase) for phase in self.pb.phaseData])

    def __phase_n_ineq(self, id, phase):
        """
        Counts the dimension of the inequalities in a phase
        - COM Kinematic constraints: summation over all effectors, times 2 because there are 2 height possible for the transition
        - Relative kinematic constraints between each effectors
        - Inequalities relative to each contact surface
        @param id index of the phase in the problem
        @param phase concerned phase
        """
        n_ineq = 0
        # COM kinematic constraints
        for _, (K, _) in enumerate(phase.K):
            n_ineq += K.shape[0]
        if id != 0:
            n_ineq *= 2

        # Foot relative distance
        for (_, Ks) in phase.allRelativeK[phase.moving]:
            n_ineq += Ks[0].shape[0]

        # Surfaces
        n_surfaces = phase.n_surfaces
        n_ineq += sum([S[1].shape[0] for S in phase.S])

        # Slack positivity
        if n_surfaces > 1:
            n_ineq += n_surfaces
        return n_ineq

    def __total_n_ineq(self):
        """
        Counts the number of inequality constraints
        @return the number of inequality constraints of the problem
        """
        return sum([self.__phase_n_ineq(i, phase) for i, phase in enumerate(self.pb.phaseData)])

    def __total_n_eq(self):
        """
        Counts the number of equality constraints
        @return the number of equality constraints of the problem
        """
        return self.default_n_equality_constraints * self.pb.n_phases

    def convert_pb_to_LP(self, pb, convert_surfaces=False):
        """
        Compute the constraints:
        G x <= h
        C x  = d

        @param pb problem data
        @return G, h, C, d
        """
        if convert_surfaces:
            replace_surfaces_with_ineq_in_problem(pb, False)

        self.set_problem(pb)

        n_variables = self.__total_n_variables()
        n_ineq = self.__total_n_ineq()
        n_eq = self.__total_n_eq()

        G = np.zeros((n_ineq, n_variables))
        h = np.zeros(n_ineq)
        C = np.zeros((n_eq, n_variables))
        d = np.zeros(n_eq)

        i_start = 0
        i_start_eq = 0
        js = [0]
        cons = Constraints(self.n_effectors)
        for id, phase in enumerate(self.pb.phaseData):
            # inequalities
            j_next = js[-1] + self.__phase_n_variables(phase)
            feet_phase = self.__feet_last_moving_phase(id)
            i_start = cons.fixed_foot_com(self.pb, phase, G, h, i_start, js, id, feet_phase)
            i_start = cons.foot_relative_distance(self.pb, phase, G, h, i_start, js, feet_phase)
            i_start = cons.surface_inequality(phase, G, h, i_start, js[-1])
            i_start = cons.slack_positivity(phase, G, h, i_start, js[-1])

            # equalities
            i_start_eq = cons.com(self.pb, phase, C, d, i_start_eq, js, id, feet_phase)

            js.append(j_next)

        G, h = normalize([G, h])
        C, d = normalize([C, d])
        return G, h, C, d

    def selected_surfaces(self, alphas):
        """
        :param: alphas the list of slack variables found by the planner
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
            alphas.append(result[j + self.default_n_variables:j + self.__phase_n_variables(phase)])
            j += self.__phase_n_variables(phase)
        return alphas

    def get_result(self, result):
        """
        Retrieve the com and feet positions from the result
        @param result vector of variables
        @return com positions, moving foot positions and all feet positions
        """
        coms = []
        moving_foot_pos = []
        all_feet_pos = [[self.pb.p0[i]] for i in range(self.n_effectors)]

        j = 0
        for i, phase in enumerate(self.pb.phaseData):
            coms.append(self.com_2.dot(result[j:j + self.default_n_variables]))
            moving_foot_pos.append(self.foot.dot(result[j:j + self.default_n_variables]))

            for foot in range(self.n_effectors):
                if foot == phase.moving:
                    all_feet_pos[foot].append(moving_foot_pos[i])
                else:
                    all_feet_pos[foot].append(all_feet_pos[foot][-1])

            j += self.__phase_n_variables(phase)

        for foot in range(self.n_effectors):
            all_feet_pos[foot].pop(0)

        return coms, moving_foot_pos, all_feet_pos

    def com_cost(self, coms):
        """
        Compute a cost to keep the com close to a target com at each phase
        @param coms list of target positions for the com
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self.__total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            A = np.zeros((2, n_variables))
            A[:, j: j + self.default_n_variables] = self.com_xy
            b = coms[id][:2]

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self.__phase_n_variables(phase)

        return P, q

    def end_com_cost(self, com):
        """
        Compute a cost to keep the final CoM position close to a target one 
        @param com Target final com position
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self.__total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        n_phases = self.pb.n_phases
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            if id == n_phases - 1:
                A = np.zeros((2, n_variables))
                A[:, j: j + self.default_n_variables] = self.com_xy
                b = com[:2]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
            j += self.__phase_n_variables(phase)

        return P, q

    def end_effectors_position_cost(self, effector_positions):
        """
        Compute a cost to keep the final end effectors positions closed to target ones 
        @param effector_positions list of each effector's final position
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self.__total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        n_phases = self.pb.n_phases
        first_phase = max(n_phases-self.pb.n_effectors, 0)
        j = 0.
        for id, phase in enumerate(self.pb.phaseData):
            if id >= first_phase:
                A = np.zeros((3, n_variables))
                A[:, j: j + self.default_n_variables] = self.foot
                b = effector_positions[phase.moving][:3]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
            j += self.__phase_n_variables(phase)

        return P, q

    def posture_cost(self):
        """
        Compute a cost to keep the feet relative positions as close as possible from the initial ones
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        relative_positions = [np.array(self.pb.p0[i] - self.pb.p0[0])
                              for i in range(1, self.n_effectors)]

        n_variables = self.__total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        js = [0]
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            feet_phase = self.__feet_last_moving_phase(id)

            A_0 = np.zeros((3, n_variables))
            b_0 = np.zeros(3)
            if phase.moving == 0:
                A_0[:, j:j + self.default_n_variables] = -self.foot
            elif feet_phase[0] != -1:
                j0 = js[feet_phase[0]]
                A_0[:, j0:j0 + self.default_n_variables] = -self.foot
            else:
                b_0 = self.pb.p0[0]

            for foot in range(1, self.n_effectors):
                A = np.copy(A_0)
                b = b_0 + relative_positions[foot-1]
                if foot == phase.moving:
                    A[:, j:j + self.default_n_variables] = self.foot
                elif feet_phase[foot] != -1:
                    jf = js[feet_phase[foot]]
                    A[:, jf:jf + self.default_n_variables] = self.foot
                else:
                    b += -self.pb.p0[foot]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self.__phase_n_variables(phase)
            js.append(j)

        return P, q

    def step_size_cost(self, step_size):
        """
        Compute a cost to keep the step sizes as close as possible of a target step size
        @param step_size desired size of the steps
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self.__total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        js = [0]
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            feet_phase = self.__feet_last_moving_phase(id)

            foot = phase.moving

            A = np.zeros((2, n_variables))
            b = step_size
            A[:, j:j + self.default_n_variables] = self.foot_xy
            if feet_phase[foot] != -1:
                jf = js[feet_phase[0]]
                A[:, jf:jf + self.default_n_variables] = -self.foot_xy
            else:
                b += self.pb.p0[foot][:2]

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self.__phase_n_variables(phase)
            js.append(j)

        return P, q

    def compute_costs(self, costs):
        """
        This function computes the P and q cost matrix and vector given a cost dictionary
        @param costs the cost dictionary. The keys should match keys from self.cost_dict
        @return P, q the cost matrix and vector of the QP problem
        """
        n_variables = self.__total_n_variables()

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
