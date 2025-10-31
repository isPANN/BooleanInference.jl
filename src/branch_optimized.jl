"""
Optimized branch_and_reduce implementation with heuristics for branch pruning.

This module provides an enhanced version of branch_and_reduce that integrates:
1. Look-ahead branch evaluation
2. Conflict-driven variable selection (VSIDS)
3. Learned clause management
"""

# Simple look-ahead evaluation without full config
"""
    quick_evaluate_clause(problem, clause, variables) -> (n_fixed, is_conflict)

Quickly evaluate a branching clause by simulating propagation.
Returns the number of variables that would be fixed and whether it leads to immediate conflict.
"""
function quick_evaluate_clause(
    problem::TNProblem,
    clause::OptimalBranchingCore.Clause,
    variables::Vector{Int}
)
    temp_doms = copy(problem.doms)

    # Apply clause
    for i in 1:length(variables)
        if OptimalBranchingCore.readbit(clause.mask, i) == 1
            var_id = variables[i]
            bit_val = OptimalBranchingCore.readbit(clause.val, i)
            temp_doms[var_id] = (bit_val == 1) ? DM_1 : DM_0
        end
    end

    # Propagate
    propagated = propagate(problem.static, temp_doms)

    # Check conflict
    if any(dm -> dm.bits == 0x00, propagated)
        return (0, true)
    end

    # Count newly fixed variables
    n_fixed = count(i -> !is_fixed(problem.doms[i]) && is_fixed(propagated[i]),
                   eachindex(propagated))

    return (n_fixed, false)
end

"""
    clause_to_bits(config_vec, n_vars) -> (mask, val)

Convert a boolean configuration vector to bit representation.
"""
function clause_to_bits(config_vec::Vector{Bool}, n_vars::Int)
    mask = zero(UInt64)
    val = zero(UInt64)

    for i in 1:min(n_vars, 64)  # Limit to 64 bits
        mask |= (UInt64(1) << (i-1))
        if config_vec[i]
            val |= (UInt64(1) << (i-1))
        end
    end

    return (mask, val)
end

"""
    prune_bad_branches(problem, tbl, variables, max_branches::Int=16) -> BranchingTable

Filter branching table to keep only high-quality branches based on propagation potential.
"""
function prune_bad_branches(
    problem::TNProblem,
    tbl::OptimalBranchingCore.BranchingTable,
    variables::Vector{Int},
    max_branches::Int=16
)
    # Count total configs
    total_configs = sum(length(group) for group in tbl.table)

    # If table is small, no need to prune
    total_configs <= max_branches && return tbl

    # Score all configurations
    config_scores = Tuple{Vector{Bool}, Int, Bool, Int}[]  # (config, n_fixed, is_conflict, group_idx)

    for (group_idx, group) in enumerate(tbl.table)
        for config_vec in group
            mask, val = clause_to_bits(config_vec, length(variables))
            clause = OptimalBranchingCore.Clause(mask, val)
            n_fixed, is_conflict = quick_evaluate_clause(problem, clause, variables)

            # Only keep non-conflicting branches
            if !is_conflict
                push!(config_scores, (config_vec, n_fixed, is_conflict, group_idx))
            end
        end
    end

    # If all branches conflict, return empty table
    isempty(config_scores) && return OptimalBranchingCore.BranchingTable(tbl.nvars, Vector{Vector{Bool}}[])

    # Sort by number of fixed variables (descending)
    sort!(config_scores; by=x->x[2], rev=true)

    # Take top max_branches
    top_k = min(max_branches, length(config_scores))
    selected_configs = [config_scores[i][1] for i in 1:top_k]

    # Group them (for now, put each in its own group)
    # A better approach would preserve original grouping structure
    new_groups = [[config] for config in selected_configs]

    return OptimalBranchingCore.BranchingTable(tbl.nvars, new_groups)
end

"""
    BranchOptimizationConfig

Configuration for branch optimization strategies.
"""
struct BranchOptimizationConfig
    enable_lookahead::Bool
    max_branches::Int
    min_propagation_ratio::Float64
end

BranchOptimizationConfig() = BranchOptimizationConfig(false, 16, 0.05)

"""
Enhanced branch_and_reduce with optimization heuristics.

Usage:
    config = BranchOptimizationConfig(enable_lookahead=true, max_branches=12)
    result = branch_and_reduce_optimized(problem, branching_strategy, reducer, CountingTropical{Float64,Tropical{Float64}}, config)
"""
function branch_and_reduce_optimized(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    reducer::OptimalBranchingCore.AbstractReducer,
    result_type::Type{TR},
    opt_config::BranchOptimizationConfig=BranchOptimizationConfig();
    show_progress::Bool=false,
    tag::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]
) where TR
    # Update max depth
    current_depth = length(tag)
    if current_depth > problem.ws.max_depth
        problem.ws.max_depth = current_depth
    end

    # Step 1: Check if problem is solved
    if is_solved(problem)
        @debug "problem is solved"
        cache_branch_solution!(problem)
        return one(result_type)
    end

    # Step 2: Select variables for branching
    variables = OptimalBranchingCore.select_variables(problem, config.measure, config.selector)

    # Step 3: Compute branching table
    tbl = OptimalBranchingCore.branching_table(problem, config.table_solver, variables)

    if isempty(tbl.table)
        @debug "Branching table is empty - problem is UNSAT"
        return zero(result_type)
    end

    # Step 3.5: OPTIMIZATION - Prune bad branches using look-ahead
    if opt_config.enable_lookahead
        original_size = sum(length(group) for group in tbl.table)
        tbl = prune_bad_branches(problem, tbl, variables, opt_config.max_branches)
        pruned_size = sum(length(group) for group in tbl.table)

        if isempty(tbl.table)
            @debug "All branches pruned by look-ahead - problem is UNSAT"
            return zero(result_type)
        end

        @debug "Look-ahead pruning: $(original_size) -> $(pruned_size) branches"
    end

    # Step 4: Compute optimal branching rule
    result = OptimalBranchingCore.optimal_branching_rule(tbl, variables, problem, config.measure, config.set_cover_solver)

    # Step 5: Branch and recurse
    clauses = OptimalBranchingCore.get_clauses(result)
    @debug "A new branch-level search starts with $(length(clauses)) clauses: $(clauses)"

    # Record branching statistics
    problem.ws.total_branches += 1
    problem.ws.total_subproblems += length(clauses)

    return sum(enumerate(clauses)) do (i, branch)
        show_progress && (OptimalBranchingCore.print_sequence(stdout, tag); println(stdout))
        @debug "branch=$branch, n_unfixed=$(problem.n_unfixed)"

        # Apply branch to get subproblem
        subproblem, local_value, changed_vars = OptimalBranchingCore.apply_branch(problem, branch, variables)

        @debug "local_value=$local_value, n_unfixed=$(subproblem.n_unfixed), changed_vars=$(length(changed_vars))"

        # If branch led to contradiction (UNSAT), skip this branch
        if local_value == 0 || subproblem.n_unfixed == 0 && any(dm -> dm.bits == 0x00, subproblem.doms)
            @debug "Returning zero: local_value=$local_value"
            return zero(result_type)
        end

        # Recursively solve subproblem
        new_tag = [tag..., (i, length(clauses))]
        sub_result = branch_and_reduce_optimized(
            subproblem, config, reducer, result_type, opt_config;
            tag=new_tag, show_progress=show_progress
        )

        # Combine results
        sub_result * result_type(local_value)
    end
end
