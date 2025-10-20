# Register candidate ownership for variables
function register_candidate_regions!(state::GammaQueueState, owner_var::Int, candidate::BranchingCandidate)
    candidate.region === nothing && return
    owner = owner_var
    var_regions = state.var_regions
    @inbounds for var_id in candidate.variables
        owners = var_regions[var_id]
        owner ∉ owners && push!(owners, owner)
    end
end

# Unregister candidate ownership when candidate is invalidated
function unregister_candidate_regions!(state::GammaQueueState, owner_var::Int, candidate::BranchingCandidate)
    candidate.region === nothing && return
    owner = owner_var
    var_regions = state.var_regions
    @inbounds for var_id in candidate.variables
        owners = var_regions[var_id]
        idx = findfirst(==(owner), owners)
        idx !== nothing && deleteat!(owners, idx)
    end
end

function invalidate_candidate!(state::GammaQueueState, var_id::Int)
    candidate = pop!(state.candidates, var_id, nothing)
    if !isnothing(candidate)
        unregister_candidate_regions!(state, var_id, candidate)
        # Clear pending flag if this candidate had a region
        delete!(state.pending_candidates, var_id)
    end
    # Lazy deletion: mark as deleted instead of removing from heap
    push!(state.deleted_vars, var_id)
end

function mark_candidates_pending!(state::GammaQueueState, vars::Vector{Int})
    var_regions = state.var_regions
    pending = state.pending_candidates
    n = length(var_regions)
    @inbounds for var_id in vars
        (var_id < 1 || var_id > n) && continue
        owners = var_regions[var_id]
        for owner_var in owners
            push!(pending, owner_var)
        end
    end
end

function gamma_queue_state(problem::TNProblem)::GammaQueueState
    ws = problem.ws
    state = ws.branch_queue
    var_num = length(problem.static.vars)
    if !(state isa GammaQueueState) || length(state.var_regions) != var_num
        # initialize the branch queue
        new_state = GammaQueueState(var_num)
        ws.branch_queue = new_state
        return new_state
    end
    return state::GammaQueueState
end

@inline function dom_versions(problem::TNProblem, vars::Vector{Int})
    versions = Vector{UInt32}(undef, length(vars))
    # Cache field access to reduce getproperty overhead
    var_versions = problem.ws.var_versions
    @inbounds for (i, var_id) in enumerate(vars)
        versions[i] = var_versions[var_id]
    end
    return versions
end

function build_branching_candidate(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    variables::Vector{Int},
    region::Union{Nothing,Region},
)
    tbl = OptimalBranchingCore.branching_table(problem, config.table_solver, variables)
    result = OptimalBranchingCore.optimal_branching_rule(
        tbl,
        variables,
        problem,
        config.measure,
        config.set_cover_solver,
    )
    focus_var = region === nothing ? (isempty(variables) ? 0 : variables[1]) : region.id
    versions = dom_versions(problem, variables)
    return BranchingCandidate(focus_var, region, variables, result, versions)
end

@inline function candidate_valid(candidate::BranchingCandidate, problem::TNProblem)::Bool
    vars = candidate.variables
    versions = candidate.versions
    length(vars) != length(versions) && return false
    
    # Cache field access to reduce getproperty overhead
    var_versions = problem.ws.var_versions
    @inbounds for (i, var_id) in enumerate(vars)
        var_versions[var_id] != versions[i] && return false
    end
    return true
end

function build_region_for_variable(problem::TNProblem, selector::LeastOccurrenceSelector, var_id::Int)
    is_fixed(problem.doms[var_id]) && return nothing
    return k_neighboring(problem.static, problem.doms, var_id; 
                        max_tensors=selector.max_tensors, k=selector.k)
end

function compute_branching_candidate_for_var(problem::TNProblem, config, selector, var_id::Int)
    region = build_region_for_variable(problem, selector, var_id)
    isnothing(region) && return nothing
    
    cache_region!(problem, region)
    variables = vcat(region.boundary_vars, region.inner_vars)
    isempty(variables) && return nothing
    
    candidate = build_branching_candidate(problem, config, variables, region)
    @debug "candidate_gamma" focus=region.id γ=candidate.result.γ
    return candidate
end

function default_branching_candidate(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
)
    variables = OptimalBranchingCore.select_variables(
        problem,
        config.measure,
        config.selector,
    )
    region = get_cached_region(problem)
    if region !== nothing
        set_last_region!(problem, region.id)
    end
    return build_branching_candidate(problem, config, variables, region)
end

function ensure_candidate!(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    selector::LeastOccurrenceSelector,
    state::GammaQueueState,
    var_id::Int,
)
    # Check if we can reuse existing candidate
    candidate = get(state.candidates, var_id, nothing)
    
    if !isnothing(candidate)
        # Cache region access to avoid repeated getproperty
        cand_region = candidate.region
        owner_pending = var_id ∈ state.pending_candidates
        
        if !owner_pending && candidate_valid(candidate, problem)
            return candidate
        end
        
        # Invalidate old candidate
        invalidate_candidate!(state, var_id)
    end
    
    # Compute new candidate
    candidate = compute_branching_candidate_for_var(problem, config, selector, var_id)
    isnothing(candidate) && return nothing
    
    # Register new candidate
    state.candidates[var_id] = candidate
    # Push to heap with negative gamma for max-heap behavior
    push!(state.queue, (-candidate.result.γ, var_id))
    delete!(state.deleted_vars, var_id)  # Clear deletion flag if present
    register_candidate_regions!(state, var_id, candidate)
    
    # Clear pending flag for the newly computed candidate's region
    delete!(state.pending_candidates, var_id)
    
    return candidate
end

function prune_fixed!(state::GammaQueueState, problem::TNProblem)
    for var_id in collect(keys(state.candidates))
        is_fixed(problem.doms[var_id]) && invalidate_candidate!(state, var_id)
    end
    return nothing
end

function prune_queue_fixed!(state::GammaQueueState, problem::TNProblem)
    # Mark fixed variables for lazy deletion
    # Cache field access to reduce getproperty overhead
    doms = problem.doms
    deleted_vars = state.deleted_vars
    for var_id in keys(state.candidates)
        is_fixed(doms[var_id]) && push!(deleted_vars, var_id)
    end
    return nothing
end

function initialize_queue_lazily!(problem, state::GammaQueueState)
    # Lazily initialize queue with all unfixed variables using a simple heuristic
    # The actual gamma will be computed only when needed
    queue = state.queue
    deleted_vars = state.deleted_vars
    for var_id in get_unfixed_vars(problem.doms)
        var_id ∈ deleted_vars && continue
        haskey(state.candidates, var_id) && continue
        # Use a simple heuristic: degree (higher degree = higher priority)
        # Store as (-heuristic, var_id) for max-heap behavior
        deg = Float64(problem.static.vars[var_id].deg)
        push!(queue, (-deg, var_id))
    end
    return nothing
end

function select_branching_candidate(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
)
    selector = config.selector
    selector isa LeastOccurrenceSelector || return default_branching_candidate(problem, config)

    state = gamma_queue_state(problem)
    prune_fixed!(state, problem)
    prune_queue_fixed!(state, problem)
    
    # Lazily initialize queue if empty
    isempty(state.queue) && initialize_queue_lazily!(problem, state)

    while !isempty(state.queue)
        # Pop from heap (lazy deletion - skip invalid entries)
        neg_gamma, focus_var = pop!(state.queue)
        
        # Skip if marked as deleted
        if focus_var ∈ state.deleted_vars
            delete!(state.deleted_vars, focus_var)
            continue
        end
        
        # Check if we already have a valid candidate
        candidate = get(state.candidates, focus_var, nothing)
        if candidate !== nothing
            owner_pending = focus_var ∈ state.pending_candidates
            
            if !owner_pending && candidate_valid(candidate, problem)
                # Valid cached candidate - use it directly
                cand_region = candidate.region
                if !isnothing(cand_region)
                    set_last_region!(problem, cand_region.id)
                end
                # No need to delete pending - already not pending
                return candidate
            end
        end

        # Candidate is missing, invalid, or pending - compute/refresh it now (lazily)
        # ensure_candidate! will delete from pending_candidates upon success
        candidate = ensure_candidate!(problem, config, selector, state, focus_var)
        
        if candidate === nothing
            # Failed to compute candidate (e.g., var became fixed) - skip and try next
            continue
        end

        # Successfully computed/refreshed candidate - use it
        cand_region = candidate.region
        if !isnothing(cand_region)
            set_last_region!(problem, cand_region.id)
        end
        # No need to delete pending - ensure_candidate! already did it
        return candidate
    end

    return default_branching_candidate(problem, config)
end

function OptimalBranchingCore.branch_and_reduce(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    reducer::OptimalBranchingCore.AbstractReducer,
    result_type::Type{TR};
    show_progress::Bool=false,
    tag::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]
) where TR
    if is_solved(problem)
        @debug "problem is solved"
        cache_branch_solution!(problem)
        return one(result_type)
    end
        
    @assert reducer isa NoReducer
    
    candidate = select_branching_candidate(problem, config)
    variables = candidate.variables
    result = candidate.result
    @debug "selected_gamma" focus=(candidate.region === nothing ? nothing : candidate.region.id) γ=result.γ
    
    clauses = OptimalBranchingCore.get_clauses(result)
    @debug "A new branch-level search starts with $(length(clauses)) clauses: $(clauses)"
    problem.ws.total_branches += 1
    problem.ws.total_subproblems += length(clauses)

    # Save workspace state before branching for rollback
    saved_var_versions = copy(problem.ws.var_versions)
    
    result = sum(enumerate(clauses)) do (i, branch)
        show_progress && (OptimalBranchingCore.print_sequence(stdout, tag); println(stdout))
        @debug "branch=$branch, n_unfixed=$(problem.n_unfixed)"
        
        # Save state before this branch
        branch_saved_versions = copy(problem.ws.var_versions)
        
        subproblem, local_value = OptimalBranchingCore.apply_branch(
            problem, 
            branch, 
            variables
        )
        
        @debug "local_value=$local_value, n_unfixed=$(subproblem.n_unfixed)"
        
        if local_value == 0 || subproblem.n_unfixed == 0 && any(dm -> dm.bits == 0x00, subproblem.doms)
            @debug "Returning zero: local_value=$local_value, n_unfixed=$(subproblem.n_unfixed), has_contradiction=$(any(dm -> dm.bits == 0x00, subproblem.doms))"
            # Restore state before returning
            copyto!(problem.ws.var_versions, branch_saved_versions)
            return zero(result_type)
        end
        
        new_tag = show_progress ? [tag..., (i, length(clauses))] : tag
        sub_result = OptimalBranchingCore.branch_and_reduce(
            subproblem, 
            config, 
            reducer, 
            result_type;
            tag=new_tag,
            show_progress=show_progress
        )
        
        # Restore state after this branch completes (for next sibling branch)
        copyto!(problem.ws.var_versions, branch_saved_versions)
        
        sub_result * result_type(local_value)
    end
    
    # Restore original state after all branches complete
    copyto!(problem.ws.var_versions, saved_var_versions)
    
    return result
end

function OptimalBranchingCore.optimal_branching_rule(
    tbl::OptimalBranchingCore.BranchingTable,
    variables::Vector{T},
    problem::TNProblem,
    measure::OptimalBranchingCore.AbstractMeasure,
    solver::OptimalBranchingCore.AbstractSetCoverSolver
) where T
    candidates = OptimalBranchingCore.bit_clauses(tbl)
    return OptimalBranchingCore.greedymerge(candidates, problem, variables, measure)
end

function OptimalBranchingCore.apply_branch(problem::TNProblem, clause, variables::Vector{Int})
    # Apply clause to fix variables
    new_doms = copy(problem.doms)
    n_fixed = 0
    
    for i in 1:length(variables)
        OptimalBranchingCore.readbit(clause.mask, i) != 1 && continue
        
        var_id = variables[i]
        is_fixed(problem.doms[var_id]) && continue
        
        bit_val = OptimalBranchingCore.readbit(clause.val, i)
        new_doms[var_id] = bit_val == 1 ? DM_1 : DM_0
        n_fixed += 1
    end
    
    @debug "apply_branch: Fixed $n_fixed variables"
    
    # Propagate and check validity
    propagated_doms = propagate(problem.static, new_doms)
    any(dm -> dm.bits == 0x00, propagated_doms) && return (create_failed_problem(problem), 0)
    
    new_n_unfixed = count_unfixed(propagated_doms)
    @debug "apply_branch: n_unfixed: $(problem.n_unfixed) -> $new_n_unfixed"
    
    # Validate progress
    (new_n_unfixed > problem.n_unfixed || (new_n_unfixed == problem.n_unfixed && n_fixed == 0)) && 
        return (create_failed_problem(problem), 0)
    
    # Single-pass: update version numbers and mark candidates pending (zero allocation)
    state = problem.ws.branch_queue
    if state isa GammaQueueState
        # Fast path with pending marking
        var_regions = state.var_regions
        pending = state.pending_candidates
        var_versions = problem.ws.var_versions
        old_doms = problem.doms
        n = length(var_regions)
        @inbounds for idx in 1:length(propagated_doms)
            propagated_doms[idx].bits == old_doms[idx].bits && continue
            # Variable changed - update version
            var_versions[idx] += 1
            # Mark affected candidates as pending
            (idx < 1 || idx > n) && continue
            owners = var_regions[idx]
            for owner_var in owners
                push!(pending, owner_var)
            end
        end
    else
        # No queue state - just update versions
        var_versions = problem.ws.var_versions
        old_doms = problem.doms
        @inbounds for idx in 1:length(propagated_doms)
            propagated_doms[idx].bits != old_doms[idx].bits && (var_versions[idx] += 1)
        end
    end
    
    # Create new problem and update caches
    new_problem = TNProblem(problem.static, propagated_doms, new_n_unfixed, problem.ws)
    clear_region_cache!(problem)
    
    return (new_problem, 1)
end

@inline create_failed_problem(problem::TNProblem) = 
    TNProblem(problem.static, fill(DomainMask(0x00), length(problem.doms)), 0, problem.ws)

function OptimalBranchingCore.reduce_problem(
    ::Type{T},
    problem::TNProblem,
    ::OptimalBranchingCore.NoReducer
) where T
    return (problem, one(T))
end
