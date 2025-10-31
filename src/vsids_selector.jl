"""
Variable State Independent Decaying Sum (VSIDS) inspired selector.
Prioritizes variables involved in recent conflicts and constraint propagations.
"""

mutable struct VSIDSSelector <: AbstractSelector
    k::Int  # number of hops for k-neighboring
    max_tensors::Int
    scores::Dict{Int,Float64}  # Variable activity scores
    decay::Float64  # Score decay factor (e.g., 0.95)
    bump_amount::Float64  # How much to increase score on conflict
end

function VSIDSSelector(k::Int=2; max_tensors::Int=2, decay::Float64=0.95, bump_amount::Float64=1.0)
    return VSIDSSelector(k, max_tensors, Dict{Int,Float64}(), decay, bump_amount)
end

"""
    bump_variable_score!(selector, var_id)

Increase the activity score of a variable (called when involved in conflict/propagation).
"""
function bump_variable_score!(selector::VSIDSSelector, var_id::Int)
    current = get(selector.scores, var_id, 0.0)
    selector.scores[var_id] = current + selector.bump_amount
end

"""
    decay_scores!(selector)

Decay all variable scores (called periodically, e.g., after each branching decision).
"""
function decay_scores!(selector::VSIDSSelector)
    for (var_id, score) in selector.scores
        selector.scores[var_id] = score * selector.decay
    end
end

"""
    bump_conflicting_variables!(selector, problem, changed_vars)

Bump scores of variables involved in recent propagation/conflict.
This should be called after constraint propagation with the list of changed variables.
"""
function bump_conflicting_variables!(selector::VSIDSSelector, changed_vars::Vector{Int})
    for var_id in changed_vars
        bump_variable_score!(selector, var_id)
    end
end

function OptimalBranchingCore.select_variables(
    problem::TNProblem,
    measure::AbstractMeasure,
    selector::VSIDSSelector
)
    # Suppress unused argument warning - measure is part of the interface
    _ = measure

    unfixed_vars = get_unfixed_vars(problem.doms)
    isempty(unfixed_vars) && return Int[]

    # Find the variable with highest activity score
    best_var_idx = argmax(unfixed_vars) do var_id
        get(selector.scores, var_id, 0.0)
    end

    # If all scores are zero, fall back to least occurrence
    if all(get(selector.scores, v, 0.0) == 0.0 for v in unfixed_vars)
        best_var_idx = argmin(length(problem.static.v2t[u]) for u in unfixed_vars)
    end

    # Compute k-neighboring region around the high-activity variable
    region = k_neighboring(
        problem.static,
        problem.doms,
        unfixed_vars[best_var_idx];
        max_tensors = selector.max_tensors,
        k = selector.k
    )
    cache_region!(problem, region)

    # Decay scores periodically
    decay_scores!(selector)

    return vcat(region.boundary_vars, region.inner_vars)
end
