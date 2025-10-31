"""
Look-ahead heuristics for intelligent branch pruning.
Evaluate potential branches before full exploration to eliminate poor choices.
"""

"""
    LookaheadConfig

Configuration for look-ahead branch evaluation:
- `enabled`: Whether to use look-ahead
- `max_lookahead_depth`: How deep to simulate propagation
- `min_fixed_ratio`: Minimum ratio of newly fixed variables to consider a branch "good"
- `max_branches`: Maximum number of branches to generate (prune rest)
"""
struct LookaheadConfig
    enabled::Bool
    max_lookahead_depth::Int
    min_fixed_ratio::Float64
    max_branches::Int
end

LookaheadConfig() = LookaheadConfig(true, 1, 0.1, 16)

"""
    evaluate_branch_propagation(problem, clause, variables) -> (n_fixed, has_conflict)

Simulate applying a branch and propagating to see how much it reduces the problem.
Returns:
- n_fixed: Number of newly fixed variables after propagation
- has_conflict: Whether this branch immediately leads to conflict
"""
function evaluate_branch_propagation(
    problem::TNProblem,
    clause::OptimalBranchingCore.Clause,
    variables::Vector{Int}
)
    # Apply the clause to create a temporary domain
    temp_doms = copy(problem.doms)
    n_directly_fixed = 0

    for i in 1:length(variables)
        if OptimalBranchingCore.readbit(clause.mask, i) == 1
            var_id = variables[i]
            bit_val = OptimalBranchingCore.readbit(clause.val, i)
            new_val = (bit_val == 1) ? DM_1 : DM_0

            if !is_fixed(problem.doms[var_id])
                temp_doms[var_id] = new_val
                n_directly_fixed += 1
            end
        end
    end

    # Propagate to see secondary effects
    propagated = propagate(problem.static, temp_doms)

    # Check for immediate conflict
    if any(dm -> dm.bits == 0x00, propagated)
        return (0, true)
    end

    # Count how many additional variables got fixed by propagation
    n_propagated_fixed = count(eachindex(propagated)) do i
        !is_fixed(problem.doms[i]) && is_fixed(propagated[i])
    end

    return (n_propagated_fixed, false)
end

"""
    score_clause(problem, clause, variables, config) -> Float64

Compute a quality score for a branching clause.
Higher score = better branch (more propagation, less search space).
"""
function score_clause(
    problem::TNProblem,
    clause::OptimalBranchingCore.Clause,
    variables::Vector{Int},
    config::LookaheadConfig
)
    config.enabled || return 0.0

    n_fixed, has_conflict = evaluate_branch_propagation(problem, clause, variables)

    # Conflicting branches get negative score (will be pruned)
    has_conflict && return -Inf

    # Score based on how much the branch reduces the problem
    # More fixed variables = better branch
    reduction_ratio = n_fixed / max(1, problem.n_unfixed)

    return reduction_ratio
end

"""
    filter_branching_table(problem, tbl, variables, config) -> BranchingTable

Filter the branching table to keep only high-quality branches.
This reduces the branching factor by eliminating branches that don't propagate well.
"""
function filter_branching_table(
    problem::TNProblem,
    tbl::OptimalBranchingCore.BranchingTable,
    variables::Vector{Int},
    config::LookaheadConfig
)
    config.enabled || return tbl

    # If table is small enough, no need to filter
    total_clauses = sum(length(group) for group in tbl.table)
    total_clauses <= config.max_branches && return tbl

    # Score all clauses
    clause_scores = Tuple{Vector{Bool}, Float64}[]

    for group in tbl.table
        for config_vec in group
            # Convert config vector to a clause for scoring
            clause = vector_to_clause(config_vec, variables)
            score = score_clause(problem, clause, variables, config)

            # Skip branches that immediately conflict
            score > -Inf && push!(clause_scores, (config_vec, score))
        end
    end

    # If all branches conflict, return empty table
    isempty(clause_scores) && return OptimalBranchingCore.BranchingTable(tbl.nvars, Vector{Vector{Int}}())

    # Sort by score (descending) and take top k
    sort!(clause_scores; by=x->x[2], rev=true)
    top_k = min(config.max_branches, length(clause_scores))

    # Keep only good branches that meet the minimum threshold
    filtered_configs = Vector{Vector{Bool}}[]
    for i in 1:top_k
        config_vec, score = clause_scores[i]
        if score >= config.min_fixed_ratio
            push!(filtered_configs, [config_vec])
        end
    end

    # If nothing meets threshold, keep at least the best one
    if isempty(filtered_configs) && !isempty(clause_scores)
        push!(filtered_configs, [clause_scores[1][1]])
    end

    return OptimalBranchingCore.BranchingTable(tbl.nvars, filtered_configs)
end

"""
    vector_to_clause(config_vec, variables) -> Clause

Convert a configuration vector to a Clause object for evaluation.
Note: variables parameter is kept for interface consistency but not used in current implementation.
"""
function vector_to_clause(config_vec::Vector{Bool}, variables::Vector{Int})
    _ = variables  # Part of interface, may be used in future
    n = length(config_vec)
    mask = zero(UInt64)
    val = zero(UInt64)

    for i in 1:n
        mask |= (UInt64(1) << (i-1))
        if config_vec[i]
            val |= (UInt64(1) << (i-1))
        end
    end

    return OptimalBranchingCore.Clause(mask, val)
end
