"""
Learned clause management for conflict-driven search.
Stores and reuses conflicts discovered during search to prune future branches.
"""

"""
    LearnedClause

Represents a learned constraint: a set of variable assignments that leads to UNSAT.
- vars: Variable IDs involved in the conflict
- vals: Values (0 or 1) that lead to conflict
- activity: How recently this clause was used (for clause database management)
"""
mutable struct LearnedClause
    vars::Vector{Int}
    vals::Vector{Bool}
    activity::Float64
end

"""
    LearnedClauseDB

Database of learned clauses with periodic cleanup.
"""
mutable struct LearnedClauseDB
    clauses::Vector{LearnedClause}
    max_size::Int
    decay_factor::Float64
end

LearnedClauseDB(max_size::Int=1000) = LearnedClauseDB(LearnedClause[], max_size, 0.95)

"""
    matches_conflict(clause, doms) -> Bool

Check if current domains match a learned conflict clause.
"""
function matches_conflict(clause::LearnedClause, doms::Vector{DomainMask})
    for (i, var_id) in enumerate(clause.vars)
        dm = doms[var_id]

        # If variable is fixed to the conflicting value, we might hit the same conflict
        if is_fixed(dm)
            expected_val = clause.vals[i]
            actual_val = has1(dm)
            if expected_val != actual_val
                return false  # Different assignment, won't match this conflict
            end
        end
    end
    return true
end

"""
    check_learned_conflicts(db, doms) -> Bool

Check if current state matches any learned conflict.
Returns true if we should prune this branch.
"""
function check_learned_conflicts(db::LearnedClauseDB, doms::Vector{DomainMask})
    for clause in db.clauses
        if matches_conflict(clause, doms)
            # Bump activity of this clause since it was useful
            clause.activity += 1.0
            return true
        end
    end
    return false
end

"""
    learn_conflict!(db, problem, branch_path)

Learn from a conflict: extract the minimal set of assignments that led to UNSAT.
Note: problem parameter is reserved for future conflict analysis improvements.
"""
function learn_conflict!(
    db::LearnedClauseDB,
    problem::TNProblem,
    branch_path::Vector{Tuple{Int,Bool}}  # Sequence of (var_id, value) decisions
)
    _ = problem  # Reserved for future conflict analysis

    # Only learn from non-trivial conflicts (at least 2 variables)
    length(branch_path) < 2 && return

    # For now, learn the full path as a conflict clause
    # TODO: Implement conflict analysis to minimize the clause
    vars = [var_id for (var_id, _) in branch_path]
    vals = [val for (_, val) in branch_path]

    learned = LearnedClause(vars, vals, 1.0)
    push!(db.clauses, learned)

    # Cleanup if database is too large
    if length(db.clauses) > db.max_size
        cleanup_clause_db!(db)
    end
end

"""
    cleanup_clause_db!(db)

Remove low-activity learned clauses to keep database size manageable.
"""
function cleanup_clause_db!(db::LearnedClauseDB)
    # Decay all activities
    for clause in db.clauses
        clause.activity *= db.decay_factor
    end

    # Keep top 50% by activity
    sort!(db.clauses; by=c->c.activity, rev=true)
    target_size = div(db.max_size, 2)
    resize!(db.clauses, target_size)
end

"""
    analyze_conflict(problem, doms) -> Vector{Int}

Analyze why the current state is UNSAT and extract the minimal conflict set.
Returns the variable IDs involved in the conflict.

This is a simplified version - a full implementation would do:
1. Trace back through propagation to find the decision variables
2. Use resolution to minimize the conflict clause
3. Apply first-UIP scheme from modern SAT solvers
"""
function analyze_conflict(problem::TNProblem, doms::Vector{DomainMask})
    _ = problem  # Will be used for implication graph tracing in future

    # Find all fixed variables that might contribute to the conflict
    fixed_vars = Int[]

    for var_id in eachindex(doms)
        if is_fixed(doms[var_id])
            push!(fixed_vars, var_id)
        end
    end

    # For now, return all fixed variables
    # A better implementation would trace the implication graph
    return fixed_vars
end
