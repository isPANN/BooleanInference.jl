struct BranchingCandidate
    focus_var::Int
    region::Union{Nothing,Region}
    variables::Vector{Int}
    table::OptimalBranchingCore.BranchingTable
    result::OptimalBranchingCore.OptimalBranchingResult
    signature::Vector{UInt8}
end

mutable struct RegionMembership
    regions::Set{Int}
end
RegionMembership() = RegionMembership(Set{Int}())

mutable struct GammaQueueState
    queue::PriorityQueue{Int,Float64}
    candidates::Dict{Int,BranchingCandidate}
    var_regions::Vector{RegionMembership}
    stale_regions::Set{Int}
    pending_regions::Set{Int}
end
function GammaQueueState(var_num::Int)
    memberships = [RegionMembership() for _ in 1:var_num]
    GammaQueueState(
        PriorityQueue{Int,Float64}(),
        Dict{Int,BranchingCandidate}(),
        memberships,
        Set{Int}(),
        Set{Int}(),
    )
end

@inline function register_candidate_regions!(
    state::GammaQueueState,
    candidate::BranchingCandidate,
)
    region = candidate.region
    region === nothing && return nothing
    region_id = region.id
    for var_id in candidate.variables
        push!(state.var_regions[var_id].regions, region_id)
    end
    return nothing
end

@inline function unregister_candidate_regions!(
    state::GammaQueueState,
    candidate::BranchingCandidate,
)
    region = candidate.region
    region === nothing && return nothing
    region_id = region.id
    for var_id in candidate.variables
        delete!(state.var_regions[var_id].regions, region_id)
    end
    return nothing
end

@inline function invalidate_candidate!(state::GammaQueueState, var_id::Int)
    candidate = pop!(state.candidates, var_id, nothing)
    if candidate !== nothing
        unregister_candidate_regions!(state, candidate)
    end
    haskey(state.queue, var_id) && delete!(state.queue, var_id)
    delete!(state.stale_regions, var_id)
    delete!(state.pending_regions, var_id)
    return nothing
end

@inline function mark_regions_stale!(state::GammaQueueState, vars::Vector{Int})
    isempty(vars) && return nothing
    for var_id in vars
        if var_id < 1 || var_id > length(state.var_regions)
            continue
        end
        memberships = state.var_regions[var_id].regions
        isempty(memberships) && continue
        for region_id in memberships
            push!(state.stale_regions, region_id)
            push!(state.pending_regions, region_id)
        end
    end
    return nothing
end

@inline function refresh_pending_regions!(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    selector::LeastOccurrenceSelector,
    state::GammaQueueState,
)
    isempty(state.pending_regions) && return nothing
    # Avoid mutating while iterating by copying ids
    pending_ids = collect(state.pending_regions)
    for region_id in pending_ids
        ensure_candidate!(problem, config, selector, state, region_id)
    end
    return nothing
end

@inline function gamma_queue_state(problem::TNProblem)
    ws = problem.ws
    state = ws.branch_queue
    var_num = length(problem.static.vars)
    if !(state isa GammaQueueState) || length(state.var_regions) != var_num
        new_state = GammaQueueState(var_num)
        ws.branch_queue = new_state
        return new_state
    end
    return state
end

@inline function dom_signature(problem::TNProblem, vars::Vector{Int})
    sig = Vector{UInt8}(undef, length(vars))
    @inbounds for (i, var_id) in enumerate(vars)
        sig[i] = problem.doms[var_id].bits
    end
    return sig
end

@inline function build_branching_candidate(
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
    signature = dom_signature(problem, variables)
    return BranchingCandidate(focus_var, region, variables, tbl, result, signature)
end

@inline function candidate_valid(candidate::BranchingCandidate, problem::TNProblem)::Bool
    vars = candidate.variables
    sig = candidate.signature
    length(vars) == length(sig) || return false
    @inbounds for (i, var_id) in enumerate(vars)
        dm = problem.doms[var_id]
        if dm.bits == 0x00 || dm.bits != sig[i]
            return false
        end
    end
    return true
end

function build_region_for_variable(
    problem::TNProblem,
    selector::LeastOccurrenceSelector,
    var_id::Int,
)
    is_fixed(problem.doms[var_id]) && return nothing
    return k_neighboring(
        problem.static,
        problem.doms,
        var_id;
        max_tensors = selector.max_tensors,
        k = selector.k,
    )
end

function compute_branching_candidate_for_var(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    selector::LeastOccurrenceSelector,
    var_id::Int,
)
    region = build_region_for_variable(problem, selector, var_id)
    region === nothing && return nothing
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
    candidate = get(state.candidates, var_id, nothing)
    needs_refresh = candidate === nothing || in(var_id, state.stale_regions)
    if !needs_refresh && candidate !== nothing
        if candidate_valid(candidate, problem)
            return candidate
        else
            needs_refresh = true
        end
    end

    if needs_refresh && candidate !== nothing
        invalidate_candidate!(state, var_id)
    end

    candidate = compute_branching_candidate_for_var(problem, config, selector, var_id)
    if candidate === nothing
        delete!(state.stale_regions, var_id)
        delete!(state.pending_regions, var_id)
        return nothing
    end

    state.candidates[var_id] = candidate
    state.queue[var_id] = candidate.result.γ
    register_candidate_regions!(state, candidate)
    delete!(state.stale_regions, var_id)
    delete!(state.pending_regions, var_id)
    return candidate
end

function prune_fixed!(state::GammaQueueState, problem::TNProblem)
    for var_id in collect(keys(state.candidates))
        is_fixed(problem.doms[var_id]) || continue
        invalidate_candidate!(state, var_id)
    end
    return nothing
end

function populate_queue!(
    problem::TNProblem,
    config::OptimalBranchingCore.BranchingStrategy,
    selector::LeastOccurrenceSelector,
    state::GammaQueueState,
)
    for var_id in get_unfixed_vars(problem.doms)
        if !haskey(state.candidates, var_id) || !haskey(state.queue, var_id)
            ensure_candidate!(problem, config, selector, state, var_id)
        end
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
    refresh_pending_regions!(problem, config, selector, state)
    isempty(state.queue) && populate_queue!(problem, config, selector, state)

    while !isempty(state.queue)
        focus_var, _ = peek(state.queue)
        candidate = ensure_candidate!(problem, config, selector, state, focus_var)
        candidate === nothing && continue

        region = candidate.region
        if region !== nothing
            set_last_region!(problem, region.id)
        end
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
    tbl = candidate.table
    result = candidate.result
    @debug "selected_gamma" focus=(candidate.region === nothing ? nothing : candidate.region.id) γ=result.γ
    
    clauses = OptimalBranchingCore.get_clauses(result)
    @debug "A new branch-level search starts with $(length(clauses)) clauses: $(clauses)"
    
    return sum(enumerate(clauses)) do (i, branch)
        show_progress && (OptimalBranchingCore.print_sequence(stdout, tag); println(stdout))
        @debug "branch=$branch, n_unfixed=$(problem.n_unfixed)"
        
        subproblem, local_value = OptimalBranchingCore.apply_branch(
            problem, 
            branch, 
            variables
        )
        
        @debug "local_value=$local_value, n_unfixed=$(subproblem.n_unfixed)"
        
        if local_value == 0 || subproblem.n_unfixed == 0 && any(dm -> dm.bits == 0x00, subproblem.doms)
            @debug "Returning zero: local_value=$local_value, n_unfixed=$(subproblem.n_unfixed), has_contradiction=$(any(dm -> dm.bits == 0x00, subproblem.doms))"
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
        
        sub_result * result_type(local_value)
    end
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

function OptimalBranchingCore.apply_branch(
    problem::TNProblem, 
    clause::OptimalBranchingCore.Clause{INT}, 
    variables::Vector{Int}
) where {INT<:Integer}
    new_doms = copy(problem.doms)
    n_fixed = 0
    for i in 1:length(variables)
        if OptimalBranchingCore.readbit(clause.mask, i) == 1
            var_id = variables[i]
            bit_val = OptimalBranchingCore.readbit(clause.val, i)
            new_val = (bit_val == 1) ? DM_1 : DM_0
            
            if !is_fixed(problem.doms[var_id])
                new_doms[var_id] = new_val
                n_fixed += 1
            end
        end
    end
    
    @debug "apply_branch: Fixed $n_fixed variables"

    propagated_doms = propagate(problem.static, new_doms)
    
    if any(dm -> dm.bits == 0x00, propagated_doms)
        @debug "apply_branch: Contradiction detected during propagation"
        return (TNProblem(problem.static, fill(DomainMask(0x00), length(propagated_doms)), 0, problem.ws), 0)
    end
    
    new_n_unfixed = count_unfixed(propagated_doms)

    changed_vars = Int[]
    @inbounds for (idx, dm) in enumerate(propagated_doms)
        if dm.bits != problem.doms[idx].bits
            push!(changed_vars, idx)
        end
    end
    
    @debug "apply_branch: n_unfixed: $(problem.n_unfixed) -> $new_n_unfixed"
    
    if new_n_unfixed > problem.n_unfixed
        @debug "apply_branch: ERROR - unfixed count increased!"
        return (TNProblem(problem.static, fill(DomainMask(0x00), length(propagated_doms)), 0, problem.ws), 0)
    elseif new_n_unfixed == problem.n_unfixed && n_fixed == 0
        @debug "apply_branch: No progress made (n_unfixed same and n_fixed=0)"
        return (TNProblem(problem.static, fill(DomainMask(0x00), length(propagated_doms)), 0, problem.ws), 0)
    end

    new_problem = TNProblem(
        problem.static,
        propagated_doms,
        new_n_unfixed,
        problem.ws
    )
    
    clear_region_cache!(problem)
    
    state = problem.ws.branch_queue
    if state isa GammaQueueState
        mark_regions_stale!(state, changed_vars)
    end
    
    return (new_problem, 1)
end

function OptimalBranchingCore.reduce_problem(
    ::Type{T},
    problem::TNProblem,
    ::OptimalBranchingCore.NoReducer
) where T
    return (problem, one(T))
end
