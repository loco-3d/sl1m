import numpy as np
import itertools
import copy

from sl1m.solver import call_QP_solver, call_LP_solver, Solvers
from sl1m.tools.plot_tools import draw_scene, plot_planner_result

from sl1m.problem_data import ProblemData

try:
    from time import perf_counter as clock
except ImportError:
    from time import clock

ALPHA_THRESHOLD = 0.01


def optimize_sparse_L1(planner, pb, costs, QP_SOLVER, LP_SOLVER):
    """
    This solver is called when the sparsity is fixed.
    It assumes the only contact surface for each phase is the one used for contact creation.
    Solve the problem with a specific solver
    @param P, q, G, h, C, d problem datalse
    @return None if wrong SOLVER, else ResultData
    """
    G, h, C, d = planner.convert_pb_to_LP(pb, False)
    P, q = planner.compute_costs(costs)
    if costs != {}:
        result = call_QP_solver(P, q, G, h, C, d, QP_SOLVER)
    else:
        result = call_LP_solver(q, G, h, C, d, LP_SOLVER)

    if result.success:
        coms, moving_foot_pos, all_feet_pos = planner.get_result(result.x)
        return ProblemData(True, result.time, coms, moving_foot_pos, all_feet_pos)
    else:
        print("optimize_sparse_L1 failed to solve the QP")
        return ProblemData(False, result.time)


def fix_sparsity_combinatorial(planner, pb, surfaces, LP_SOLVER):
    """
    Calls the sl1m solver. Tries to solve non fixed sparsity by handling the combinatorial.
    Ultimately calls solve which provides the approriate cost function
    @param planner
    @param pb problem data
    @param surfaces potential surfaces
    @param SOLVER Solver choice
    @return true if the problem was solved, fixed surfaces problem, surface_indices and time
    """
    G, h, C, d = planner.convert_pb_to_LP(pb)
    q = 100. * planner.alphas

    result = call_LP_solver(q, G, h, C, d, LP_SOLVER)
    t = result.time
    if not result.success:
        print("Initial LP solver fails in fix_sparsity_combinatorial")
        return False, pb, [], t

    alphas = planner.get_alphas(result.x)
    if is_sparsity_fixed(pb, alphas):
        surface_indices = planner.selected_surfaces(alphas)
        for i, phase in enumerate(pb.phaseData):
            for j, surfaces in enumerate(phase.S):
                surfaces = [surfaces[surface_indices[i][j]]]
                phase.n_surfaces[j] = len(surfaces)
        return True, pb, surface_indices, t

    pbs = generate_fixed_sparsity_problems(pb, alphas)
    if pbs is None:
        print("No combinatorial problems was found")
        return False, pb, [], t

    # Handle the combinatorial
    sparsity_fixed = False
    i = 0
    for (fixed_pb, _, solution_indices) in pbs:
        G, h, C, d = planner.convert_pb_to_LP(fixed_pb, False)
        q = 100. * planner.alphas
        result = call_LP_solver(q, G, h, C, d, LP_SOLVER)
        t += result.time
        if result.success:
            alphas = planner.get_alphas(result.x)
            if is_sparsity_fixed(fixed_pb, alphas):
                sparsity_fixed = True
                break
        i += 1

    if not sparsity_fixed:
        print("Sparsity could not be fixed")
        return False, fixed_pb, [], t

    fixed_pb_surface_indices = planner.selected_surfaces(alphas)
    j = 0
    surface_indices = []
    for i, phase in enumerate(fixed_pb.phaseData):
        if len(phase.S) > 1:
            phase.S = [phase.S[fixed_pb_surface_indices[i]]]
            phase.n_surfaces = len(phase.S)
            surface_indices.append(fixed_pb_surface_indices[i])
        elif len(surfaces[i]) > 1:
            surface_indices.append(solution_indices[j])
            j += 1
        else:
            surface_indices.append(0)



    return sparsity_fixed, fixed_pb, surface_indices, t


def get_undecided_surfaces(pb, alphas):
    """
    Get the surfaces and indices of all the undecided surfaces
    @param planner the planner
    @param pb the problem data
    @return the phase indices, sorted potential surfaces and sorted surface indices
    """
    indices = []
    surfaces = []
    surfaces_indices = []
    for i, phase in enumerate(pb.phaseData):
        for j, n_surface in enumerate(phase.n_surfaces):
            if n_surface > 1:
                if np.array(alphas[i][j]).min() > ALPHA_THRESHOLD:
                    indices.append((i, j))
                    sorted_surfaces = np.argsort(alphas[i][j])
                    surfaces_indices += [sorted_surfaces]
                    surfaces += [[[phase.S[j][idx]] for idx in sorted_surfaces]]
    return indices, surfaces, surfaces_indices


def is_sparsity_fixed(pb, alphas):
    """
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return true if the sparsity is fixed (ie there is one and only one alpha~0 per phase)
    """
    indices, _, _ = get_undecided_surfaces(pb, alphas)
    return len(indices) == 0


def generate_fixed_sparsity_problems(pb, alphas):
    """
    Check if the combinatorial is not too big, if not return all the problems
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of problems
    """
    indices, surfaces, surfaces_indices = get_undecided_surfaces(pb, alphas)
    all_len = [len(s) for s in surfaces]
    n_pbs = 1
    for l in all_len:
        n_pbs *= l
    if n_pbs > 2000:
        print("Problem probably too big to handle combinatorial", n_pbs)
        return None
    print("Handling combinatorial: ", n_pbs)
    return generate_combinatorials(pb, indices, surfaces, surfaces_indices)


def generate_combinatorials(pb, indices, surfaces, surface_indices):
    """
    Generate all the problems with only one potential surface per undecided phase
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of problems, the phase indices, and the indices of the selected surfaces
    """
    pbs = []
    sorted_combinations = [el for el in itertools.product(*surface_indices)]
    all_indices = [[el for el in range(lens)] for lens in [len(surfs) for surfs in surfaces]]

    combinations = [c for c in itertools.product(*all_indices)]
    for j, combination in enumerate(combinations):
        fixed_pb = copy.deepcopy(pb)
        for i, idx in enumerate(indices):
            fixed_pb.phaseData[idx].S = surfaces[i][combination[i]]
            fixed_pb.phaseData[idx].n_surfaces = 1
        pbs += [[fixed_pb, indices, sorted_combinations[j]]]
    return pbs
