
import numpy as np
from sl1m.planner_biped import BipedPlanner
from sl1m.planner_generic import Planner
from sl1m.solver import call_MIP_solver, Solvers, solve_MIP_gurobi_cost
from sl1m.fix_sparsity import fix_sparsity_combinatorial, fix_sparsity_combinatorial_gait, optimize_sparse_L1
from sl1m.problem_data import ProblemData


# ----------------------- L1 -----------------------------------------------------------------------

def solve_L1_combinatorial(pb, surfaces, lp_solver=Solvers.GUROBI, qp_solver=Solvers.GUROBI, costs={}, com=True):
    """
    Solve the problem by first chosing the surfaces with a L1 norm minimization problem handling the
    combinatorial if necesary, and then optimizing the feet positions with a QP
    @param pb problem to solve
    @surfaces surfaces to choose from
    @lp_solver solver to use for the LP
    @qp_solver solver to use for the QP
    @costs cost dictionary specifying the cost functions to use and their parameters
    @return ProblemData storing the result
    """
    planner = Planner(mip=False, com=com)
    sparsity_fixed, pb, surface_indices, t = fix_sparsity_combinatorial_gait(planner, pb, surfaces, lp_solver)
    if sparsity_fixed:
        pb_data = optimize_sparse_L1(planner, pb, costs, qp_solver, lp_solver)
        pb_data.surface_indices = surface_indices
    else:
        return ProblemData(False, t)
    pb_data.time += t
    return pb_data


def solve_L1_combinatorial_biped(pb, surfaces, lp_solver=Solvers.GUROBI, qp_solver=Solvers.GUROBI, costs={}):
    """
    Solve the problem for a biped by first chosing the surfaces with a L1 norm minimization problem
    handling the combinatorial if necesary, and then optimizing the feet positions with a QP
    @param pb problem to solve
    @surfaces surfaces to choose from
    @lp_solver solver to use for the LP
    @qp_solver solver to use for the QP
    @costs cost dictionary specifying the cost functions to use and their parameters
    @return ProblemData storing the result
    """
    planner = BipedPlanner()
    sparsity_fixed, pb, surface_indices, t = fix_sparsity_combinatorial(
        planner, pb, surfaces, lp_solver)
    if sparsity_fixed:
        pb_data = optimize_sparse_L1(planner, pb, costs, qp_solver, lp_solver)
        pb_data.surface_indices = surface_indices
    else:
        return ProblemData(False, t)
    pb_data.time += t
    return pb_data


# ----------------------- MIP -----------------------------------------------------------------------

def solve_MIP(pb, costs={}, solver=Solvers.GUROBI, com=False):
    """
    Solve the problem with a MIP solver
    @param pb problem to solve
    @surfaces surfaces to choose from
    @costs cost dictionary specifying the cost functions to use and their parameters
    @solver MIP solver to use
    @return ProblemData storing the result
    """
    planner = Planner(mip=True, com=com)
    G, h, C, d = planner.convert_pb_to_LP(pb)
    slack_selection_vector = planner.alphas
    P = None

    # If no combinatorial call directly a QP
    if solver == Solvers.CVXPY and np.linalg.norm(slack_selection_vector) < 1:
        return optimize_sparse_L1(planner, pb, costs)

    q = None
    if costs != None:
        P, q = planner.compute_costs(costs)
    result = call_MIP_solver(slack_selection_vector, P, q, G, h, C, d, solver=solver)

    if costs != None and solver == Solvers.CVXPY:
        alphas = planner.get_alphas(result.x)
        selected_surfaces = planner.selected_surfaces(alphas)
        for i, phase in enumerate(pb.phaseData):
            for j in range(len(phase.n_surfaces)):
                phase.S[j] = [phase.S[j][selected_surfaces[i][j]]]
                phase.n_surfaces[j] = 1
        return optimize_sparse_L1(planner, pb, costs, QP_SOLVER=Solvers.CVXPY, LP_SOLVER=Solvers.CVXPY)

    if result.success:
        alphas = planner.get_alphas(result.x)
        coms, moving_foot_pos, all_feet_pos = planner.get_result(result.x)
        surface_indices = planner.selected_surfaces(alphas)
        return ProblemData(True, result.time, coms, moving_foot_pos, all_feet_pos, surface_indices)
    return ProblemData(False, result.time)


def solve_MIP_biped(pb, costs={}, solver=Solvers.GUROBI):
    """
    Solve the problem with a MIP solver for a biped
    @param pb problem to solve
    @surfaces surfaces to choose from
    @costs cost dictionary specifying the cost functions to use and their parameters
    @solver MIP solver to use
    @return ProblemData storing the result
    """
    planner = BipedPlanner()
    G, h, C, d = planner.convert_pb_to_LP(pb)
    slack_selection_vector = planner.alphas

    if costs != None:
        P, q = planner.compute_costs(costs)
        result = solve_MIP_gurobi_cost(slack_selection_vector, P, q, G, h, C, d)
    else:
        result = call_MIP_solver(slack_selection_vector, G, h, C, d, solver=solver)

    if result.success:
        alphas = planner.get_alphas(result.x)
        coms, moving_foot_pos, all_feet_pos = planner.get_result(result.x)
        surface_indices = planner.selected_surfaces(alphas)
        return ProblemData(True, result.time, coms, moving_foot_pos, all_feet_pos, surface_indices)
    return ProblemData(False, result.time)
