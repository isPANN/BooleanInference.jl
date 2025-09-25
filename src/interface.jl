function convert_cnf_to_bip(cnf::CNF)
    return convert_sat_to_bip(Satisfiability{ProblemReductions.SAT}(cnf))
end

function convert_sat_to_bip(sat::ConstraintSatisfactionProblem)
    problem = GenericTensorNetwork(sat)
    # Get labels of all tensors in the circuit
    he2v = getixsv(problem.code)
    tensors = GenericTensorNetworks.generate_tensors(Tropical(1.0), problem)
    vec_tensors = [vec(t) for t in tensors]
    new_tensors = [replace(t, Tropical(1.0) => zero(Tropical{Float64})) for t in vec_tensors]
    return BooleanInferenceProblem(new_tensors, he2v, length(problem.problem.symbols)), problem.problem.symbols
end

function convert_circuit_to_bip(cir::Circuit)
    return convert_sat_to_bip(CircuitSAT(cir; use_constraints=true))
end

function solve_sat_problem(sat::ConstraintSatisfactionProblem; bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1, 1), measure=NumOfVertices()), reducer=NoReducer())
    p, syms = convert_sat_to_bip(sat)
    return solve_boolean_inference_problem(p; bsconfig, reducer)
end

function solve_factoring(n::Int, m::Int, N::Int; bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1, 1), measure=NumOfVertices()), reducer=NoReducer())
    # global BRANCHNUMBER = 0
    fproblem = Factoring(m, n, N)
    res = reduceto(CircuitSAT, fproblem)
    problem = CircuitSAT(res.circuit.circuit; use_constraints=true)
    ans, vals = solve_sat_problem(problem; bsconfig, reducer)
    a, b = ProblemReductions.read_solution(fproblem, [vals[res.p]..., vals[res.q]...])
    # @show BRANCHNUMBER
    return a, b
end

function solve_sat_with_assignments(sat::ConstraintSatisfactionProblem)
    res, vals = solve_sat_problem(sat)
    return res, Dict(zip(sat.symbols, vals))
end

function solve_boolean_inference_problem(bip::BooleanInferenceProblem; bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KaHyParSelector(10), measure=NumOfVertices()), reducer=NoReducer())
    # Initialize branching status
    bs = initialize_branching_status(bip)
    bs = deduction_reduce(bip, bs, collect(1:length(bip.he2v)))
    ns, res, count_num = branch_and_reduce(bip, bs, bsconfig, reducer)
    answer = get_answer(res, bip.literal_num)
    return ns, answer
end

# """
#     solve_boolean_inference_problem_with_contraction(bip::BooleanInferenceProblem; 
#         bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1,1), measure=NumOfVertices()), 
#         reducer=NoReducer(), 
#         contraction_threshold::Int=10)

# Enhanced solver that uses direct tensor contraction for small subproblems.
# When the number of undecided variables falls below `contraction_threshold`, 
# it switches to direct tensor network contraction instead of further branching.
# This can significantly reduce the number of branches for problems with small subproblems.
# """
# function solve_boolean_inference_problem_with_contraction(bip::BooleanInferenceProblem;
#     bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1, 1), measure=NumOfVertices()),
#     reducer=NoReducer(),
#     contraction_threshold::Int=10)
#     # Initialize branching status
#     bs = initialize_branching_status(bip)
#     bs = deduction_reduce(bip, bs, collect(1:length(bip.he2v)))
#     ns, res, count_num = branch_and_reduce_with_contraction(bip, bs, bsconfig, reducer; contraction_threshold)
#     answer = get_answer(res, bip.literal_num)
#     return ns, answer
# end

# """
#     solve_sat_problem_with_contraction(sat::ConstraintSatisfactionProblem; 
#         bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1,1), measure=NumOfVertices()), 
#         reducer=NoReducer(), 
#         contraction_threshold::Int=10)

# Enhanced SAT solver that uses direct tensor contraction for small subproblems.
# """
# function solve_sat_problem_with_contraction(sat::ConstraintSatisfactionProblem;
#     bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1, 1), measure=NumOfVertices()),
#     reducer=NoReducer(),
#     contraction_threshold::Int=10)
#     p, syms = convert_sat_to_bip(sat)
#     return solve_boolean_inference_problem_with_contraction(p; bsconfig, reducer, contraction_threshold)
# end

# """
#     solve_factoring_with_contraction(n::Int, m::Int, N::Int; 
#         bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1,1), measure=NumOfVertices()), 
#         reducer=NoReducer(), 
#         contraction_threshold::Int=10)

# Enhanced factoring solver that uses direct tensor contraction for small subproblems.
# """
# function solve_factoring_with_contraction(n::Int, m::Int, N::Int;
#     bsconfig::BranchingStrategy=BranchingStrategy(table_solver=TNContractionSolver(), selector=KNeighborSelector(1, 1), measure=NumOfVertices()),
#     reducer=NoReducer(),
#     contraction_threshold::Int=10)
#     fproblem = Factoring(m, n, N)
#     res = reduceto(CircuitSAT, fproblem)
#     problem = CircuitSAT(res.circuit.circuit; use_constraints=true)
#     ans, vals = solve_sat_problem_with_contraction(problem; bsconfig, reducer, contraction_threshold)
#     a, b = ProblemReductions.read_solution(fproblem, [vals[res.p]..., vals[res.q]...])
#     return a, b
# end