import numpy as np
import itertools
import copy

from sl1m.solver import call_QP_solver, call_LP_solver

from sl1m.problem_data import ProblemData

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
        print ("QP ")
        result = call_QP_solver(P, q, G, h, C, d, QP_SOLVER)
    else:
        result = call_LP_solver(q, G, h, C, d, LP_SOLVER)

    if result.success:
        coms, moving_foot_pos, all_feet_pos = planner.get_result(result.x)
        return ProblemData(True, result.time, coms, moving_foot_pos, all_feet_pos)
    else:
        print("optimize_sparse_L1 failed to solve the QP")
        return ProblemData(False, result.time)


def fix_sparsity_combinatorial_gait(planner, pb, surfaces, LP_SOLVER):
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
    if is_sparsity_fixed_gait(pb, alphas):
        selected_surfaces = planner.selected_surfaces(alphas)
        for i, phase in enumerate(pb.phaseData):
            for j in range(len(phase.n_surfaces)):
                phase.S[j] = [phase.S[j][selected_surfaces[i][j]]]
                phase.n_surfaces[j] = 1
        return True, pb, selected_surfaces, t

    pbs = generate_fixed_sparsity_problems_gait(pb, alphas)
    if pbs is None:
        print("No combinatorial problems was found")
        return False, pb, [], t

    # Handle the combinatorial
    sparsity_fixed = False
    i = 0
    for (fixed_pb, combination) in pbs:
        G, h, C, d = planner.convert_pb_to_LP(fixed_pb, False)
        q = 100. * planner.alphas
        result = call_LP_solver(q, G, h, C, d, LP_SOLVER)
        t += result.time
        if result.success:
            alphas = planner.get_alphas(result.x)
            if is_sparsity_fixed_gait(fixed_pb, alphas):
                sparsity_fixed = True
                pb = fixed_pb
                break
        i += 1

    if not sparsity_fixed:
        print("Sparsity could not be fixed")
        return False, pb, [], t

    selected_surfaces = planner.selected_surfaces(alphas)
    k = 0
    surface_indices = []
    for i, phase in enumerate(fixed_pb.phaseData):
        for j in range(len(phase.n_surfaces)):
            if len(phase.S[j]) > 1:
                phase.S[j] = [phase.S[j][selected_surfaces[i][j]]]
                phase.n_surfaces[j] = 1
                surface_indices.append(selected_surfaces[i][j])
            elif len(surfaces[i][j]) > 1:
                surface_indices.append(combination[k])
                k += 1
            else:
                surface_indices.append(0)

    return sparsity_fixed, fixed_pb, surface_indices, t


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
    print(alphas)
    if is_sparsity_fixed(pb, alphas):
        surface_indices = planner.selected_surfaces(alphas)
        for i, phase in enumerate(pb.phaseData):
            phase.S = [phase.S[surface_indices[i]]]
            phase.n_surfaces = 1
        return True, pb, surface_indices, t

    pbs = generate_fixed_sparsity_problems(pb, alphas)
    if pbs is None:
        print("No combinatorial problems was found")
        return False, pb, [], t

    # Handle the combinatorial
    sparsity_fixed = False
    i = 0
    for (fixed_pb, combination) in pbs:
        G, h, C, d = planner.convert_pb_to_LP(fixed_pb, False)
        q = 100. * planner.alphas
        result = call_LP_solver(q, G, h, C, d, LP_SOLVER)
        t += result.time
        if result.success:
            alphas = planner.get_alphas(result.x)
            if is_sparsity_fixed(fixed_pb, alphas):
                sparsity_fixed = True
                pb = fixed_pb
                break
        i += 1

    if not sparsity_fixed:
        print("Sparsity could not be fixed")
        return False, pb, [], t

    selected_surfaces = planner.selected_surfaces(alphas)
    j = 0
    surface_indices = []
    for i, phase in enumerate(pb.phaseData):
        if phase.n_surfaces > 1:
            phase.S = [phase.S[selected_surfaces[i]]]
            phase.n_surfaces = len(phase.S)
            surface_indices.append(selected_surfaces[i])
        elif len(surfaces[i]) > 1:
            surface_indices.append(combination[j])
            j += 1
        else:
            surface_indices.append(0)

    return sparsity_fixed, fixed_pb, surface_indices, t


def get_undecided_surfaces_gait(pb, alphas):
    """
    Get the surfaces and indices of all the undecided surfaces
    @param planner the planner
    @param pb the problem data
    @return a list of phase indices, and sorted surface indices
    """
    undecided_surfaces = []
    for i, phase in enumerate(pb.phaseData):
        for j, n_surface in enumerate(phase.n_surfaces):
            if n_surface > 1 and np.array(alphas[i][j]).min() > ALPHA_THRESHOLD:
                undecided_surfaces.append([i, j, np.argsort(alphas[i][j])])
    return undecided_surfaces


def get_undecided_surfaces(pb, alphas):
    """
    Get the surfaces and indices of all the undecided surfaces
    @param planner the planner
    @param pb the problem data
    @return a list of phase indices, and sorted surface indices
    """
    undecided_surfaces = []
    for i, phase in enumerate(pb.phaseData):
        if phase.n_surfaces > 1 and np.array(alphas[i]).min() > ALPHA_THRESHOLD:
            undecided_surfaces.append([i, np.argsort(alphas[i])])
    return undecided_surfaces


def is_sparsity_fixed_gait(pb, alphas):
    """
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return true if the sparsity is fixed (ie there is one and only one alpha~0 per phase)
    """
    return len(get_undecided_surfaces_gait(pb, alphas)) == 0


def is_sparsity_fixed(pb, alphas):
    """
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return true if the sparsity is fixed (ie there is one and only one alpha~0 per phase)
    """
    return len(get_undecided_surfaces(pb, alphas)) == 0


def generate_fixed_sparsity_problems_gait(pb, alphas):
    """
    Check if the combinatorial is not too big, if not return all the problems
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of problems
    """
    undecided_surfaces = get_undecided_surfaces_gait(pb, alphas)
    all_len = [len(undecided_surface[-1]) for undecided_surface in undecided_surfaces]
    n_pbs = 1
    for l in all_len:
        n_pbs *= l
    if n_pbs > 2000:
        print("Problem probably too big to handle combinatorial", n_pbs)
        return None
    print("Handling combinatorial: ", n_pbs)
    return generate_combinatorials_gait(pb, undecided_surfaces)


def generate_fixed_sparsity_problems(pb, alphas):
    """
    Check if the combinatorial is not too big, if not return all the problems
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of problems
    """
    undecided_surfaces = get_undecided_surfaces(pb, alphas)
    all_len = [len(undecided_surface[-1]) for undecided_surface in undecided_surfaces]
    n_pbs = 1
    for l in all_len:
        n_pbs *= l
    if n_pbs > 2000:
        print("Problem probably too big to handle combinatorial", n_pbs)
        return None
    print("Handling combinatorial: ", n_pbs)
    return generate_combinatorials(pb, undecided_surfaces)


def generate_combinatorials(pb, undecided_surfaces):
    """
    Generate all the problems with only one potential surface per undecided phase
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of problems, the phase indices, and the indices of the selected surfaces
    """
    pbs = []
    sorted_combinations = list(itertools.product(*[s[-1] for s in undecided_surfaces]))
    for combination in sorted_combinations:
        fixed_pb = copy.deepcopy(pb)
        for i, undecided_surface in enumerate(undecided_surfaces):
            phase = fixed_pb.phaseData[undecided_surface[0]]
            phase.S = [phase.S[combination[i]]]
            phase.n_surfaces = 1
        pbs.append([fixed_pb, combination])
    return pbs


def generate_combinatorials_gait(pb, undecided_surfaces):
    """
    Generate all the problems with only one potential surface per undecided phase
    @param pb the problem data
    @param alphas the list of slack variables found by the planner
    @return the list of fixed problems with the indices of the selected surfaces
    """
    pbs = []
    sorted_combinations = list(itertools.product(*[s[-1] for s in undecided_surfaces]))
    for combination in sorted_combinations:
        fixed_pb = copy.deepcopy(pb)
        for i, undecided_surface in enumerate(undecided_surfaces):
            phase = fixed_pb.phaseData[undecided_surface[0]]
            effector = undecided_surface[1]
            phase.S[effector] = [phase.S[effector][combination[i]]]
            phase.n_surfaces[effector] = 1
        pbs.append([fixed_pb, combination])
    return pbs
