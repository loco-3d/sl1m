import numpy as np
from sl1m.constants_and_tools import *
from sl1m.constraints import Constraints


class Planner:
    """
    Implements an optimization problem to find the next surfaces to use

    The problem is a LP with these variables for each phase:
    [ com_x, com_y, com_z_1, com_z_2, p_i_x, p_i_y, p_i_z, {a_i}                                ]
    [ 1    , 1    , 1      , 1      , 1    , 1    , 1    , 0 if n_surfaces == 1, else n_surfaces]

    Under the following constraints :
    - fixed_foot_com: Ensures the COM is 'above' the fixed feet
    - foot_relative_distance: Ensures the moving feet is close enough to the other feet
    - surface: Each foot belongs to one surface
    - slack_positivity: The slack variables are positive
    - com_weighted_equality: Fix the horizontal position of the CoM at the barycenter of the contact points (fixed feet)
    """

    def __init__(self, mip=False, com=False):
        self.mip = mip
        self.com = com

        self.default_n_equality_constraints = 2 * int(self.com)
        self.default_n_variables = 4 * int(self.com)

        self.cost_dict = {
            "final_com": self.end_com_cost,
            "end_effector_positions": self.end_effectors_position_cost,
            "effector_positions": self.effector_position_cost,
            "coms_xy": self.com_cost_xy,
            "coms_z": self.com_cost_z,
            "coms_3D": self.com_cost_3D,
            "posture": self.posture_cost,
            "step_size": self.step_size_cost,
        }

        self.com_costs = ["final_com", "coms_xy", "coms_z", "coms_3D"]

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
        M[:, j : j + size] = np.identity(size)
        return M

    def com_xy(self, phase):
        """
        Generate a selection matrix for the com x and y components
        @param phase phase data
        @return a (2, phase number of variables (without slacks)) matrix
        """
        return self._expression_matrix(2, self._default_n_variables(phase), 0)

    def com_z(self, phase):
        """
        Generate a selection matrix for the com x and y components
        @param phase phase data
        @return a (2, phase number of variables (without slacks)) matrix
        """
        return self._expression_matrix(2, self._default_n_variables(phase), 2)

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

    def com_xyz(self, phase):
        """
        Generate a selection matrix for the height of com_2
        @param phase phase data
        @return a (1, phase number of variables (without slacks)) matrix
        """
        return self._expression_matrix(4, self._default_n_variables(phase), 0)

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

    def _slack_selection_vector(self):
        """
        Get the selection vector of the slack variables
        @return the selection vector
        """
        n_variables = self._total_n_variables()
        selection_vector = np.zeros(n_variables)
        i = 0
        for phase in self.pb.phaseData:
            n_variables_phase = self._phase_n_variables(phase)
            n_slacks = n_variables_phase - self._default_n_variables(phase)
            i_start = i + self._default_n_variables(phase)
            selection_vector[i_start : i_start + n_slacks] = np.ones(n_slacks)
            i += n_variables_phase
        return selection_vector

    def set_problem(self, pb):
        """
        Set the problem data
        @param pb new problem data
        """
        self.pb = pb
        self.n_effectors = pb.n_effectors
        self.alphas = self._slack_selection_vector()

    def _feet_last_moving_phase(self, phase_id):
        """
        Ids of the last phase each foot moved at a given phase
        @param phase_id id of the phase
        @return a list of n_effetors terms corresponding to the last phase the feet moved, -1 if it hasn't moved yet
        """
        feet_last_moving_phase = [-1] * self.n_effectors
        for phase in self.pb.phaseData:
            if phase.id >= phase_id:
                break
            for foot in phase.moving:
                feet_last_moving_phase[foot] = phase.id
        return feet_last_moving_phase

    def _phase_n_variables(self, phase):
        """
        Counts the number of variables in a phase
        @param phase concerned phase
        @return number of variables in the phase
        """
        n_variables = self._default_n_variables(phase)
        for n_surface in phase.n_surfaces:
            if n_surface > 1:
                n_variables += n_surface
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
        n_ineq = 0
        # COM kinematic constraints
        if self.com:
            for foot, (K, _) in enumerate(phase.K):
                if foot in phase.stance:
                    n_ineq += K.shape[0]
                if phase.id == 0:
                    n_ineq += K.shape[0]
                elif foot in self.pb.phaseData[phase.id - 1].stance:
                    n_ineq += K.shape[0]

        # Foot relative distance
        for foot in phase.moving:
            for (other, Ks) in phase.allRelativeK[foot]:
                if other in phase.stance:
                    n_ineq += Ks[0].shape[0]
                elif phase.id == 0:
                    n_ineq += Ks[0].shape[0]
                elif other in self.pb.phaseData[phase.id - 1].stance:
                    n_ineq += Ks[0].shape[0]

        # Surfaces
        n_ineq += sum(
            [sum([S[1].shape[0] for S in foot_surfaces]) for foot_surfaces in phase.S]
        )

        # Slack positivity
        for n_surface in phase.n_surfaces:
            if n_surface > 1:
                n_ineq += n_surface

        return n_ineq

    def _total_n_ineq(self):
        """
        Counts the number of inequality constraints
        @return the number of inequality constraints of the problem
        """
        return sum([self._phase_n_ineq(phase) for phase in self.pb.phaseData])

    def _total_n_eq(self):
        """
        Counts the number of equality constraints
        @return the number of equality constraints of the problem
        """
        n_eq = self.default_n_equality_constraints * self.pb.n_phases
        if self.mip:
            for phase in self.pb.phaseData:
                for n_surface in phase.n_surfaces:
                    if n_surface > 1:
                        n_eq += 1
        return n_eq

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

        n_variables = self._total_n_variables()
        n_ineq = self._total_n_ineq()
        n_eq = self._total_n_eq()

        G = np.zeros((n_ineq, n_variables))
        h = np.zeros(n_ineq)
        C = np.zeros((n_eq, n_variables))
        d = np.zeros(n_eq)

        i_start = 0
        i_start_eq = 0
        js = [0]
        cons = Constraints(self.n_effectors, self.com)
        for phase in self.pb.phaseData:
            feet_phase = self._feet_last_moving_phase(phase.id)

            i_start = cons.foot_relative_distance(
                self.pb, phase, G, h, i_start, js, feet_phase
            )
            i_start = cons.surface_inequality(phase, G, h, i_start, js[-1])
            i_start = cons.slack_positivity(phase, G, h, i_start, js[-1])

            if self.com:
                i_start = cons.fixed_foot_com(
                    self.pb, phase, G, h, i_start, js, feet_phase
                )
                # i_start_eq = cons.com(self.pb, phase, C, d, i_start_eq, js, feet_phase)

            if self.mip:
                i_start_eq = cons.slack_equality(phase, C, d, i_start_eq, js[-1])

            js.append(js[-1] + self._phase_n_variables(phase))

        G, h = normalize([G, h])
        C, d = normalize([C, d])
        return G, h, C, d

    def selected_surfaces(self, alphas):
        """
        :param: alphas nested list of slack variables for each foot in each phase
        Return the list of selected surfaces indices in the problem
        """
        indices = []
        for phase in self.pb.phaseData:
            phase_indices = []
            for i, n_surface in enumerate(phase.n_surfaces):
                if n_surface == 1:
                    phase_indices.append(0)
                else:
                    phase_indices.append(np.argmin(alphas[phase.id][i]))
            indices.append(phase_indices)
        return indices

    def get_alphas(self, result):
        """
        Retrieve the alphas from the result
        @param result vector of variables
        @return list of alpha values for each phase
        """
        alphas = []
        j_alpha = 0
        for phase in self.pb.phaseData:
            alpha_phase = []
            j_alpha += self._default_n_variables(phase)
            for n_surface in phase.n_surfaces:
                if n_surface == 1:
                    alpha_phase.append(None)
                else:
                    alpha_phase.append(result[j_alpha : j_alpha + n_surface])
                    j_alpha += n_surface
            alphas.append(alpha_phase)
        return alphas

    def get_result(self, result):
        """
        Retrieve the com and feet positions from the result
        @param result vector of variables
        @return com positions, moving foot positions and all feet positions
        """
        if self.com:
            coms = []
        else:
            coms = None
        moving_feet_pos = []
        all_feet_pos = [[self.pb.p0[i]] for i in range(self.n_effectors)]

        j = 0
        for i, phase in enumerate(self.pb.phaseData):
            if self.com:
                coms.append(
                    self.com_2(phase).dot(
                        result[j : j + self._default_n_variables(phase)]
                    )
                )
            phase_moving_feet = []
            for foot in phase.moving:
                phase_moving_feet.append(
                    self.foot(phase, foot).dot(
                        result[j : j + self._default_n_variables(phase)]
                    )
                )
            moving_feet_pos.append(phase_moving_feet)

            for foot in range(self.n_effectors):
                if foot in phase.moving:
                    id = np.argmax(phase.moving == foot)
                    all_feet_pos[foot].append(moving_feet_pos[i][id])
                elif foot in phase.stance:
                    all_feet_pos[foot].append(all_feet_pos[foot][-1])
                else:
                    all_feet_pos[foot].append(None)

            j += self._phase_n_variables(phase)

        return coms, moving_feet_pos, all_feet_pos

    def com_cost_xy(self, coms):
        """
        Compute a cost to keep the com in X,Y axis close to a target com at each phase
        @param coms list of target positions for the com
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for phase in self.pb.phaseData:
            A = np.zeros((2, n_variables))
            A[:, j : j + self._default_n_variables(phase)] = self.com_xy(phase)
            b = coms[phase.id][:2]

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self._phase_n_variables(phase)

        return P, q

    def com_cost_z(self, coms):
        """
        Compute a cost to keep the com in Z axis close to a target com at each phase
        @param coms list of target positions for the com
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for phase in self.pb.phaseData:
            A = np.zeros((2, n_variables))
            b = np.zeros(2)
            A[:, j : j + self._default_n_variables(phase)] = self.com_z(phase)
            b[0] = coms[phase.id][2]
            b[1] = coms[phase.id][2]

            P += np.dot(A.T, A)
            q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self._phase_n_variables(phase)

        return P, q

    def com_cost_3D(self, coms):
        """
        Compute a cost to keep the com close to a target com in X,Y,Z axis at each phase
        @param coms list of target positions for the com
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for phase in self.pb.phaseData:
            A = np.zeros((4, n_variables))
            b = np.zeros(4)
            A[:, j : j + self._default_n_variables(phase)] = self.com_xyz(phase)
            b[:3] = coms[phase.id][:]
            b[3] = coms[phase.id][2]

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

        j = 0
        for phase in self.pb.phaseData:
            if phase.id == self.pb.n_phases - 1:
                A = np.zeros((2, n_variables))
                A[:, j : j + self._default_n_variables(phase)] = self.com_xy(phase)
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

        js = [0]
        for phase in self.pb.phaseData:
            if phase.id < self.pb.n_phases - 1:
                js.append(js[-1] + self._phase_n_variables(phase))

        feet_phase = self._feet_last_moving_phase(self.pb.n_phases - 1)
        for foot in range(self.n_effectors):
            phase = self.pb.phaseData[-1]
            if foot in phase.moving:
                j = js[phase.id]
                A = np.zeros((3, n_variables))
                A[:, j : j + self._default_n_variables(phase)] = self.foot(phase, foot)
                b = effector_positions[foot][:3]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
            elif feet_phase[foot] != -1:
                phase = self.pb.phaseData[feet_phase[foot]]
                j = js[phase.id]
                A = np.zeros((3, n_variables))
                A[:, j : j + self._default_n_variables(phase)] = self.foot(phase, foot)
                b = effector_positions[foot][:3]

                P += np.dot(A.T, A)
                q += -np.dot(A.T, b).reshape(A.shape[1])
        return P, q

    def posture_cost(self):
        """
        Compute a cost to keep the feet relative positions as close as possible from the initial ones
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        relative_positions = [
            [None for _ in range(self.n_effectors)] for _ in range(self.n_effectors)
        ]
        for foot in range(self.n_effectors):
            for other in range(self.n_effectors):
                if (
                    foot != other
                    and self.pb.p0[foot] is not None
                    and self.pb.p0[other] is not None
                ):
                    relative_positions[foot][other] = np.array(
                        self.pb.p0[foot] - self.pb.p0[other]
                    )

        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        js = [0]
        j = 0
        for id, phase in enumerate(self.pb.phaseData):
            feet_phase = self._feet_last_moving_phase(id)

            for foot in range(self.n_effectors):
                for other in range(self.n_effectors):
                    if other > foot and relative_positions[foot][other] is not None:
                        A = np.zeros((3, n_variables))
                        b = relative_positions[foot][other]
                        if foot in phase.moving:
                            A[:, j : j + self._default_n_variables(phase)] = self.foot(
                                phase, foot
                            )
                        elif feet_phase[foot] != -1:
                            jf = js[feet_phase[foot]]
                            phase_f = self.pb.phaseData[feet_phase[foot]]
                            A[
                                :, jf : jf + self._default_n_variables(phase_f)
                            ] = self.foot(phase_f, foot)
                        else:
                            b -= self.pb.p0[foot]

                        if other in phase.moving:
                            A[:, j : j + self._default_n_variables(phase)] = -self.foot(
                                phase, other
                            )
                        elif feet_phase[other] != -1:
                            jf = js[feet_phase[other]]
                            phase_f = self.pb.phaseData[feet_phase[other]]
                            A[
                                :, jf : jf + self._default_n_variables(phase_f)
                            ] = -self.foot(phase_f, other)
                        else:
                            b += self.pb.p0[other]

                        P += np.dot(A.T, A)
                        q += -np.dot(A.T, b).reshape(A.shape[1])

            j += self._phase_n_variables(phase)
            js.append(j)

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

        js = [0]
        j = 0
        for phase in self.pb.phaseData:
            feet_phase = self._feet_last_moving_phase(phase.id)
            for foot in phase.moving:
                A = np.zeros((2, n_variables))
                b = step_size.copy()
                A[:, j : j + self._default_n_variables(phase)] = self.foot_xy(
                    phase, foot
                )
                if feet_phase[foot] != -1:
                    j_f = js[feet_phase[foot]]
                    phase_f = self.pb.phaseData[feet_phase[foot]]
                    A[
                        :, j_f : j_f + self._default_n_variables(phase_f)
                    ] = -self.foot_xy(phase_f, foot)
                elif self.pb.p0[foot] is not None:
                    b += self.pb.p0[foot][:2]
                else:
                    continue

                P += np.dot(A.T, A)
                q -= np.dot(A.T, b).reshape(A.shape[1])

            j += self._phase_n_variables(phase)
            js.append(j)

        return P, q

    def effector_position_cost(self, end_effector_positions):
        """
        Compute a cost to keep the effectors with a specified distance to the shoulders
        @param end_effector_positions desired effector positions
        @return P matrix and q vector s.t. we minimize x' P x + q' x
        """
        n_variables = self._total_n_variables()
        P = np.zeros((n_variables, n_variables))
        q = np.zeros(n_variables)

        j = 0
        for phase in self.pb.phaseData:
            # feet_phase = self._feet_last_moving_phase(phase.id)
            for foot in phase.moving:
                A = np.zeros((2, n_variables))
                A[:, j : j + self._default_n_variables(phase)] = self.foot_xy(
                    phase, foot
                )
                b = end_effector_positions[foot][phase.id]

                P += np.dot(A.T, A)
                q -= np.dot(A.T, b).reshape(A.shape[1])

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
                if not self.com and key in self.com_costs:
                    print(
                        "COM cost given, but no COM variable in planner. The cost is ignored"
                    )
                elif key in self.cost_dict.keys():
                    P_, q_ = self.cost_dict[key](*val[1:])
                    P += P_ * val[0]
                    q += q_ * val[0]
                else:
                    print("Unknown cost to add to the biped problem")
        return P, q
