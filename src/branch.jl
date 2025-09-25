function OptimalBranchingCore.apply_branch(p::BooleanInferenceProblem, bs::AbstractBranchingStatus, clause::Clause{INT}, vertices::Vector{T}) where {INT<:Integer,T<:Integer}
    bs_new, aedges = decide_literal(bs, p, vertices, clause)
    return deduction_reduce(p, bs_new, aedges)
end

function OptimalBranchingCore.branch_and_reduce(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, config::BranchingStrategy, reducer::AbstractReducer)
    # Determines if search should terminate
    stopped, res, count_num = check_stopped(bs)
    stopped && return res, bs, count_num

    # branch the problem
    # select a subset of variables
    subbip = select_variables(problem, bs, config.measure, config.selector)
    # compute the BranchingTable
    tbl = branching_table(problem, bs, config.table_solver, subbip)
    iszero(tbl.bit_length) && return false, bs, 1

    result = optimal_branching_rule(tbl, subbip.vs, bs, problem, config.measure, config.set_cover_solver)  # compute the optimal branching rule
    for branch in OptimalBranchingCore.get_clauses(result)
        res, bs_new, count_num1 = branch_and_reduce(problem, apply_branch(problem, bs, branch, subbip.vs), config, reducer)
        count_num += count_num1
        if res
            return res, bs_new, count_num
        end
    end
    return false, bs, count_num
end

# """
#     branch_and_reduce_with_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, config::BranchingStrategy, reducer::AbstractReducer; contraction_threshold::Int=10)

# Enhanced branching function that uses direct tensor contraction for small subproblems.
# When the number of undecided variables falls below `contraction_threshold`, 
# it switches to direct tensor network contraction instead of further branching.
# """
# function branch_and_reduce_with_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus,config::BranchingStrategy, reducer::AbstractReducer; contraction_threshold::Int=10)
#     # Determines if search should terminate
#     stopped, res, count_num = check_stopped(bs)
#     stopped && return res, bs, count_num

#     # Check if problem size is below threshold for direct contraction
#     if should_use_direct_contraction(problem, bs, config; threshold=contraction_threshold)
#         @info "Using direct tensor contraction for problem size $(problem.literal_num)"
#         return direct_tensor_contraction(problem, bs)
#     end

#     # branch the problem
#     # select a subset of variables
#     subbip = select_variables(problem, bs, config.measure, config.selector)
#     # compute the BranchingTable
#     tbl = branching_table(problem, bs, config.table_solver, subbip)
#     iszero(tbl.bit_length) && return false,bs,1

#     result = optimal_branching_rule(tbl, subbip.vs, bs,problem, config.measure, config.set_cover_solver)  # compute the optimal branching rule
#     for branch in OptimalBranchingCore.get_clauses(result)
#         res, bs_new ,count_num1= branch_and_reduce_with_contraction(problem, apply_branch(problem,bs, branch, subbip.vs), config, reducer; contraction_threshold)
#         count_num += count_num1
#         if res
#             return res, bs_new, count_num
#         end
#     end
#     return false, bs,count_num
# end

function check_stopped(bs::AbstractBranchingStatus)
    # global BRANCHNUMBER
    # if there is no clause, then the problem is solved.
    if all(bs.undecided_literals .== -1)
        # BRANCHNUMBER += 1
        return true, true, 1
    end
    if any(bs.undecided_literals .== 0)
        # BRANCHNUMBER += 1
        return true, false, 1
    end
    return false, false, 0
end

# function OptimalBranchingCore.optimal_branching_rule(table::BranchingTable, variables::Vector, bs::AbstractBranchingStatus,p::BooleanInferenceProblem, m::AbstractMeasure, solver::AbstractSetCoverSolver)
#     candidates = collect(candidate_clauses(table))
#     size_reductions = [measure(bs, m) - measure((apply_branch(p,bs, candidate, variables)), m) for candidate in candidates]
#     return minimize_γ(table, candidates, size_reductions, solver; γ0=2.0)
# end

function mybranch_and_reduce(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, config::BranchingStrategy, reducer::AbstractReducer)
    stopped, res = check_stopped(bs)
    stopped && return res, bs

    v = findfirst(x -> readbit(bs.decided_mask, x) == 0, 1:problem.literal_num)

    res, bs_new = branch_and_reduce(problem, apply_branch(problem, bs, Clause(0b1, 0b1), [v]), config, reducer)
    if res
        return res, bs_new
    end

    res, bs_new = branch_and_reduce(problem, apply_branch(problem, bs, Clause(0b1, 0b0), [v]), config, reducer)
    if res
        return res, bs_new
    end

    return false, bs
end

function OptimalBranchingCore.optimal_branching_rule(table::BranchingTable, variables::Vector, bs::AbstractBranchingStatus, p::BooleanInferenceProblem, m::AbstractMeasure, solver::AbstractSetCoverSolver)
    candidates = OptimalBranchingCore.bit_clauses(table)
    return OptimalBranchingCore.greedymerge(candidates, p, bs, variables, m)
end

# TODO: NOT DRY!
function OptimalBranchingCore.greedymerge(cls::Vector{Vector{Clause{INT}}}, problem::AbstractProblem, bs::AbstractBranchingStatus, variables::Vector, m::AbstractMeasure) where {INT}
    active_cls = collect(1:length(cls))
    cls = copy(cls)
    merging_pairs = [(i, j) for i in active_cls, j in active_cls if i < j]
    n = length(variables)
    size_reductions = [OptimalBranchingCore.size_reduction(problem, m, bs, candidate[1], variables) for candidate in cls]
    γ = OptimalBranchingCore.complexity_bv(size_reductions)
    while !isempty(merging_pairs)
        i, j = popfirst!(merging_pairs)
        if i in active_cls && j in active_cls
            for ii in 1:length(cls[i]), jj in 1:length(cls[j])
                if OptimalBranchingCore.bdistance(cls[i][ii], cls[j][jj]) == 1
                    cl12 = OptimalBranchingCore.gather2(n, cls[i][ii], cls[j][jj])
                    if cl12.mask == 0
                        continue
                    end
                    l12 = OptimalBranchingCore.size_reduction(problem, m, bs, cl12, variables)
                    if γ^(-size_reductions[i]) + γ^(-size_reductions[j]) >= γ^(-l12) + 1e-8
                        push!(cls, [cl12])
                        k = length(cls)
                        deleteat!(active_cls, findfirst(==(i), active_cls))
                        deleteat!(active_cls, findfirst(==(j), active_cls))
                        for ii in active_cls
                            push!(merging_pairs, (ii, k))
                        end
                        push!(active_cls, k)
                        push!(size_reductions, l12)
                        γ = OptimalBranchingCore.complexity_bv(size_reductions[active_cls])
                        break
                    end
                end
            end
        end
    end
    return [cl[1] for cl in cls[active_cls]]
end

function OptimalBranchingCore.size_reduction(p::AbstractProblem, m::AbstractMeasure, bs::AbstractBranchingStatus, cl::Clause{INT}, variables::Vector) where {INT}
    return measure(bs, m) - measure(apply_branch(p, bs, cl, variables), m)
end

# """
#     should_use_direct_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, config::BranchingStrategy; threshold::Int=10) -> Bool

# Determines whether to use direct enumeration instead of further branching.
# Only returns true for very small subproblems where enumeration is feasible.
# """
# function should_use_direct_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, config::BranchingStrategy; threshold::Int=10)
#     # Count undecided variables
#     undecided_count = problem.literal_num - count_ones(bs.decided_mask)

#     # Only use direct enumeration for very small problems
#     # Threshold should be much smaller than the original problem size
#     effective_threshold = min(threshold, 15)  # Cap at 15 variables max

#     return undecided_count <= threshold
# end

# """
#     direct_tensor_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus) -> (Bool, AbstractBranchingStatus, Int)

# Directly solve the remaining subproblem using tensor network contraction with OMEinsum.
# This uses the SubBIP infrastructure to create a tensor network for all remaining variables
# and contracts it to find valid assignments efficiently.
# """
# function direct_tensor_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus)
#     # Get undecided variables
#     undecided_variables = Vector{Int}()
#     undecided_variables_length = count(==(0), bs.decided_mask)
#     sizehint!(undecided_variables, undecided_variables_length)

#     for i in 1:problem.literal_num
#         if readbit(bs.decided_mask, i) == 0
#             push!(undecided_variables, i)
#         end
#     end

#     if isempty(undecided_variables)
#         # All variables decided, check if solution is valid
#         return all(bs.undecided_literals .>= 0), bs, 1
#     end

#     # Create a SubBIP for all remaining undecided variables
#     # This will construct the tensor network for the remaining subproblem
#     subbip = SubBIP(problem, bs, undecided_variables)

#     # Use tensor network contraction to find valid assignments
#     return solve_with_tensor_contraction(problem, bs, subbip)
# end

# """
#     solve_with_tensor_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, subbip::SubBIP) -> (Bool, AbstractBranchingStatus, Int)

# Use tensor network contraction to solve the subproblem represented by SubBIP.
# This function contracts the tensor network to find all valid assignments and returns the first one found.
# """
# function solve_with_tensor_contraction(problem::BooleanInferenceProblem, bs::AbstractBranchingStatus, subbip::SubBIP)
#     # Contract the entire tensor network to get all feasible assignments
#     # The result will be a tensor where each entry corresponds to a variable assignment
#     # and the value indicates whether that assignment is feasible (Tropical(0.0)) or not
#     contracted_tensor = subbip.sub_tensors

#     # Find all feasible assignments (entries equal to Tropical(0.0))
#     feasible_indices = findall(==(Tropical(0.0)), contracted_tensor)

#     if isempty(feasible_indices)
#         # No feasible assignment found
#         return false, bs, 1
#     end

#     # Take the first feasible assignment
#     first_feasible = feasible_indices[1]

#     # Convert the linear index back to variable assignments
#     # The tensor dimensions correspond to variables in subbip.vs
#     tensor_size = size(contracted_tensor)
#     assignment = ind2sub_custom(tensor_size, first_feasible)

#     # Apply the assignment to create the solution
#     bs_solution = BranchingStatus(bs.config, bs.decided_mask, copy(bs.undecided_literals))
#     for (i, var) in enumerate(subbip.vs)
#         # Convert tensor index (1 or 2) to boolean value (false or true)
#         value = assignment[i] == 2
#         bs_solution = apply_single_variable_assignment(bs_solution, problem, var, value)
#     end

#     return true, bs_solution, 1
# end

# """
# Helper function to convert linear index to subscripts (similar to ind2sub in MATLAB)
# """
# function ind2sub_custom(dims::Tuple, linear_idx::CartesianIndex)
#     return [linear_idx[i] for i in 1:length(dims)]
# end

# """
# Helper function to apply a single variable assignment to branching status
# """
# function apply_single_variable_assignment(bs::AbstractBranchingStatus, problem::BooleanInferenceProblem, var::Int, value::Bool)
#     # This is a simplified version - you may need to adapt based on your exact BranchingStatus implementation
#     clause = value ? Clause(0b1, 0b1) : Clause(0b1, 0b0)
#     bs_new, _ = decide_literal(bs, problem, [var], clause)
#     return bs_new
# end

# """
# Helper function to check if boundary assignment has feasible internal assignments
# """
# function has_feasible_internal_assignment(subbip::SubBIP, boundary_assignment::Vector{Bool})
#     out_vs_num = length(subbip.outside_vs_ind)
#     vs_num = length(subbip.vs)

#     # Create index mapping
#     ind_pos = [i ∈ subbip.outside_vs_ind ? findfirst(==(i), subbip.outside_vs_ind) : findfirst(==(i), setdiff(1:vs_num, subbip.outside_vs_ind)) for i in 1:vs_num]

#     # Convert boundary assignment to tensor indices
#     out_index = [bit ? 2 : 1 for bit in boundary_assignment]

#     # Create tensor access vector
#     vec = [var_idx ∈ subbip.outside_vs_ind ? out_index[ind_pos[var_idx]] : (:) for var_idx in 1:vs_num]

#     # Check for feasible internal assignments
#     feasible_indices = findall(==(Tropical(0.0)), subbip.sub_tensors[vec...])

#     return !isempty(feasible_indices)
# end

# """
# Helper function to apply boundary assignment and find a complete solution
# """
# function apply_boundary_assignment(bs::AbstractBranchingStatus, problem::BooleanInferenceProblem, subbip::SubBIP, boundary_assignment::Vector{Bool}, config_idx::Int)
#     # This is a placeholder - you'll need to implement the full logic
#     # to apply both boundary and internal variable assignments
#     bs_new = BranchingStatus(bs.config, bs.decided_mask, copy(bs.undecided_literals))

#     # Apply boundary assignments
#     for (i, var_idx) in enumerate(subbip.outside_vs_ind)
#         var = subbip.vs[var_idx]
#         value = boundary_assignment[i]
#         bs_new = apply_single_variable_assignment(bs_new, problem, var, value)
#     end

#     # Find and apply internal assignments
#     # (This would need the same logic as in the tensor contraction solver)

#     return bs_new
# end
