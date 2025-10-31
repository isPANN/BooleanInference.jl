using OptimalBranchingCore: BranchingTable

struct TNContractionSolver <: AbstractTableSolver end
struct TNPropagationSolver <: AbstractTableSolver end

function _build_axismap(output_vars::Vector{Int})
    isempty(output_vars) && return Int[]
    max_var = maximum(output_vars)
    axismap = zeros(Int, max_var)
    @inbounds for (idx, var_id) in enumerate(output_vars)
        axismap[var_id] = idx
    end
    return axismap
end

function separate_fixed_free_boundary(region::Region, doms::Vector{DomainMask})
    fixed_vars = Int[]; fixed_vals = Int[]; free_vars = Int[]; free_indices = Int[]
    
    for (i, var_id) in enumerate(region.boundary_vars)
        if is_fixed(doms[var_id])
            push!(fixed_vars, var_id)
            push!(fixed_vals, has1(doms[var_id]) ? 1 : 0)
        else
            push!(free_vars, var_id)
            push!(free_indices, i)
        end
    end
    return fixed_vars, fixed_vals, free_vars, free_indices
end

function construct_boundary_config(region::Region, doms::Vector{DomainMask}, free_pos_of_boundary::Vector{Int}, free_config::Int)
    n_boundary = length(region.boundary_vars)
    boundary_config = Vector{Bool}(undef, n_boundary)
    
    @inbounds for i in 1:n_boundary
        var_id = region.boundary_vars[i]
        if is_fixed(doms[var_id])
            boundary_config[i] = has1(doms[var_id])
        else
            j = free_pos_of_boundary[i]
            bit = (free_config >> (j-1)) & 0x1
            boundary_config[i] = (bit == 1)
        end
    end
    
    return boundary_config
end

function construct_inner_config(region::Region, doms::Vector{DomainMask}, inner_posmap::Vector{Int}, free_inner_configs::Vector{Vector{Bool}})
    n_inner = length(region.inner_vars)
    n_config = length(free_inner_configs)
    inner_configs = Vector{Vector{Bool}}(undef, n_config)

    @inbounds for config_idx in 1:n_config
        inner_config = Vector{Bool}(undef, n_inner)
        for i in 1:n_inner
            var_id = region.inner_vars[i]
            if is_fixed(doms[var_id])
                inner_config[i] = has1(doms[var_id])
            else
                j = inner_posmap[var_id]
                inner_config[i] = free_inner_configs[config_idx][j]
            end
        end
        inner_configs[config_idx] = inner_config
    end
    
    return inner_configs
end

function extract_inner_configs(contracted::AbstractArray{T, N}, n_inner::Int) where {T, N}
    @assert n_inner == N
    configs = Vector{Vector{Bool}}()
    one_tropical = one(Tropical{Float64})
    
    @inbounds for lin in eachindex(contracted)
        if contracted[lin] == one_tropical
            config = Vector{Bool}(undef, n_inner)
            linear_idx = LinearIndices(contracted)[lin] - 1
            for i in 1:n_inner
                config[i] = ((linear_idx >> (i-1)) & 0x1) == 1
            end
            push!(configs, config)
        end
    end
    return configs
end

function extract_boundary_configs(contracted::AbstractArray{T, N}, n_boundary::Int) where {T, N}
    @assert n_boundary == N
    configs = Vector{Vector{Bool}}()
    one_tropical = one(Tropical{Float64})
    
    @inbounds for lin in eachindex(contracted)
        if contracted[lin] == one_tropical
            config = Vector{Bool}(undef, n_boundary)
            linear_idx = LinearIndices(contracted)[lin] - 1
            for i in 1:n_boundary
                config[i] = ((linear_idx >> (i-1)) & 0x1) == 1
            end
            push!(configs, config)
        end
    end
    return configs
end

function combine_configs(boundary_config::Vector{Bool}, inner_configs::Vector{Vector{Bool}})
    n_configs = length(inner_configs)
    n_boundary = length(boundary_config)
    n_inner = isempty(inner_configs) ? 0 : length(inner_configs[1])
    n_total = n_boundary + n_inner
    
    full_configs = Vector{Vector{Bool}}(undef, n_configs)
    
    @inbounds for i in 1:n_configs
        full_config = Vector{Bool}(undef, n_total)
        copyto!(full_config, 1, boundary_config, 1, n_boundary)
        copyto!(full_config, n_boundary + 1, inner_configs[i], 1, n_inner)
        full_configs[i] = full_config
    end
    
    return full_configs
end

function get_region_contraction(problem::TNProblem, region::Region)
    cached = get_cached_region_contraction(problem; region_id=region.id)
    if isnothing(cached)
        @debug "contract_region_cache_miss" focus=region.id
        # Cache the contraction with all region variables unfixed
        # This way the cache can be reused across different domain states
        region_vars_set = Set(vcat(region.boundary_vars, region.inner_vars))
        unfixed_doms = copy(problem.doms)
        for var_id in region_vars_set
            unfixed_doms[var_id] = DomainMask(0x03)  # Unfix region variables
        end

        contracted, output_vars = contract_region(problem.static, region, unfixed_doms)
        cache_region_contraction!(problem, contracted, output_vars; region_id=region.id)

        # Now slice the cached tensor according to current doms
        return slice_cached_contraction(contracted, output_vars, region, problem.doms)
    else
        @debug "contract_region_cache_hit" focus=region.id
        cached_tensor, cached_vars = cached
        # Slice the cached tensor according to current doms
        return slice_cached_contraction(cached_tensor, cached_vars, region, problem.doms)
    end
end

"""
    slice_cached_contraction(tensor, output_vars, region, doms)

Slice a cached tensor contraction according to current domain constraints.
The cached tensor was computed with all region variables unfixed.
"""
function slice_cached_contraction(
    tensor::AbstractArray{Tropical{Float64}},
    output_vars::Vector{Int},
    region::Region,
    doms::Vector{DomainMask}
)
    # Build list of variables to fix (those in region that are now fixed in doms)
    fixed_vars = Tuple{Int,Bool}[]
    for (i, var_id) in enumerate(output_vars)
        if is_fixed(doms[var_id])
            push!(fixed_vars, (var_id, has1(doms[var_id])))
        end
    end

    if isempty(fixed_vars)
        # No variables to fix, return tensor as-is
        return tensor, output_vars
    end

    # Build axis map
    axismap = zeros(Int, maximum(output_vars))
    for (i, var_id) in enumerate(output_vars)
        axismap[var_id] = i
    end

    # Slice the tensor
    sliced = slice_region_contraction(tensor, fixed_vars, axismap)

    # Update output_vars to only include unfixed variables
    new_output_vars = Int[]
    for var_id in output_vars
        if !is_fixed(doms[var_id])
            push!(new_output_vars, var_id)
        end
    end

    # If all variables are now fixed, the result should be a scalar
    if isempty(new_output_vars)
        @assert ndims(sliced) == 0 || length(sliced) == 1
        scalar_val = ndims(sliced) == 0 ? sliced : sliced[1]
        return [scalar_val], Int[]
    end

    return collect(sliced), new_output_vars
end

function slice_region_contraction(
    tensor::AbstractArray{Tropical{Float64}},
    assignments::Vector{Tuple{Int,Bool}},
    axismap::Vector{Int},
)
    isempty(assignments) && return tensor
    nd = ndims(tensor)
    nd == 0 && return tensor

    indices = Any[Colon() for _ in 1:nd]
    
    @inbounds for (var_id, value) in assignments
        axis = axismap[var_id]
        @assert axis > 0 "Boundary variable $var_id not found in contraction axes"
        indices[axis] = value ? 2 : 1
    end

    return @view tensor[indices...]
end

function find_inner_configs(
    problem::TNProblem,
    region::Region,
    boundary_assignments::Vector{Tuple{Int,Bool}}
)
    # Create a temporary domain mask with boundary variables fixed
    temp_doms = copy(problem.doms)
    @inbounds for (var_id, value) in boundary_assignments
        temp_doms[var_id] = value ? DM_1 : DM_0
    end
    
    # Contract the region with boundary variables fixed
    # This will slice out the boundary variables and return only inner variables
    n_tensors = length(region.tensors)
    tn = problem.static
    
    sliced_tensors = Vector{Vector{Tropical{Float64}}}(undef, n_tensors)
    tensor_indices = Vector{Vector{Int}}(undef, n_tensors)
    
    @inbounds for (i, tensor_id) in enumerate(region.tensors)
        original_tensor = tn.tensors[tensor_id]
        var_axes = original_tensor.var_axes
        tensor_data = original_tensor.tensor
        
        # Slice according to fixed variables (including newly fixed boundary vars)
        sliced = slicing(tensor_data, temp_doms, var_axes)
        sliced_tensors[i] = sliced
        
        # Collect remaining unfixed variables (only inner vars should remain)
        remaining_vars = Int[]
        for var_id in var_axes
            if !is_fixed(temp_doms[var_id])
                push!(remaining_vars, var_id)
            end
        end
        tensor_indices[i] = remaining_vars
    end
    
    # Collect unfixed inner variables
    inner_output_vars = Int[]
    for var_id in region.inner_vars
        if !is_fixed(temp_doms[var_id])
            push!(inner_output_vars, var_id)
        end
    end
    
    if isempty(inner_output_vars)
        # All inner variables are fixed
        contracted = contract_tensors(sliced_tensors, tensor_indices, Int[])
        # Should be satisfiable (we already checked the boundary config)
        return Vector{Bool}[[]]
    end
    
    # Contract to get inner variable tensor
    contracted = contract_tensors(sliced_tensors, tensor_indices, inner_output_vars)
    
    # Extract satisfying inner configurations
    free_inner_configs = extract_inner_configs(contracted, length(inner_output_vars))
    
    if isempty(free_inner_configs)
        # This shouldn't happen if the boundary config is satisfiable
        @warn "No inner configs found for satisfiable boundary config"
        return Vector{Bool}[[]]
    end
    
    # Construct full inner configurations (including fixed inner vars)
    inner_posmap = _build_axismap(inner_output_vars)
    inner_configs = construct_inner_config(region, temp_doms, inner_posmap, free_inner_configs)
    
    return inner_configs
end

function handle_no_boundary_case(
    problem::TNProblem,
    region::Region,
    contracted::AbstractArray{Tropical{Float64}},
)
    # When there are no boundary variables, all variables are inner variables
    # After contraction, inner variables are contracted out, so the result is a scalar
    # If the scalar is 1 (tropical), the problem is satisfiable
    @assert ndims(contracted) == 0 || (length(contracted) == 1) "Expected scalar for no boundary case"
    
    n_inner = length(region.inner_vars)
    one_tropical = one(Tropical{Float64})
    scalar_value = ndims(contracted) == 0 ? contracted : contracted[1]
    
    if scalar_value == one_tropical
        # Problem is satisfiable. Need to find actual inner variable assignments.
        # Since boundary vars are all fixed, we can contract with those fixings to find inner vars.
        # But actually, if n_inner == 0, there are no variables to assign!
        if n_inner == 0
            return BranchingTable(0, [[]])  # Empty config for empty problem
        end
        
        # Otherwise, we need to find the inner configs by further contracting
        # For now, return a trivial satisfying assignment (all zeros)
        # This is a simplification - ideally we should extract the actual configs
        dummy_config = zeros(Bool, n_inner)
        return BranchingTable(n_inner, [[dummy_config]])
    else
        # Problem is unsatisfiable - return completely empty table
        return BranchingTable(n_inner, Vector{Vector{Int}}())
    end
end

function OptimalBranchingCore.branching_table(
    problem::TNProblem, 
    solver::TNContractionSolver, 
    variables::Vector{T}
) where T
    region = get_cached_region(problem)
    isnothing(region) && error("No cached region found. Make sure `select_variables` is called first.")
    
    n_boundary = length(region.boundary_vars)
    n_inner = length(region.inner_vars)
    n_total = n_boundary + n_inner
    @assert n_total == length(variables)
    
    contracted_tensor, output_vars = get_region_contraction(problem, region)
    axismap = _build_axismap(output_vars)

    # After the contraction optimization, output_vars should only contain boundary variables
    # Inner variables have been contracted out
    @assert all(var_id in region.boundary_vars for var_id in output_vars) "Output vars should only contain boundary variables"

    n_boundary == 0 && return handle_no_boundary_case(problem, region, contracted_tensor)
    
    # Extract all satisfiable boundary configurations directly from the contracted tensor
    # This is much more efficient than enumerating all 2^n_boundary configurations
    _, _, free_boundary_vars, free_boundary_indices = separate_fixed_free_boundary(region, problem.doms)
    n_free_boundary = length(free_boundary_vars)
    
    # Extract satisfiable free boundary variable configurations
    free_boundary_configs = extract_boundary_configs(contracted_tensor, n_free_boundary)
    
    # If no satisfiable configurations, return empty table
    # Use empty vector to indicate UNSAT (no valid configurations at all)
    isempty(free_boundary_configs) && return BranchingTable(n_total, Vector{Vector{Int}}())
    
    # Build position mapping for constructing full boundary configs
    free_pos_of_boundary = zeros(Int, n_boundary)
    @inbounds for (j, boundary_idx) in enumerate(free_boundary_indices)
        free_pos_of_boundary[boundary_idx] = j
    end
    
    valid_config_groups = Vector{Vector{Vector{Bool}}}()
    
    # For each satisfiable boundary configuration
    @inbounds for free_boundary_config in free_boundary_configs
        # Construct the full boundary configuration (including fixed vars)
        boundary_config = Vector{Bool}(undef, n_boundary)
        for i in 1:n_boundary
            var_id = region.boundary_vars[i]
            if is_fixed(problem.doms[var_id])
                boundary_config[i] = has1(problem.doms[var_id])
            else
                j = free_pos_of_boundary[i]
                boundary_config[i] = free_boundary_config[j]
            end
        end
        
        # Build assignments for finding inner configs
        assignments = Tuple{Int,Bool}[]
        for (j, var_id) in enumerate(free_boundary_vars)
            push!(assignments, (var_id, free_boundary_config[j]))
        end
        
        # Find satisfying inner variable assignments for this boundary config
        inner_configs = find_inner_configs(problem, region, assignments)
        
        # Combine boundary and inner configs
        full_configs = combine_configs(boundary_config, inner_configs)
        
        push!(valid_config_groups, full_configs)
    end
    
    # If no valid config groups, return empty table (UNSAT)
    isempty(valid_config_groups) && return BranchingTable(n_total, Vector{Vector{Int}}())
    return BranchingTable(n_total, valid_config_groups)
end

# ============================================================================
# Propagation-based branching table (more efficient alternative)
# ============================================================================

"""
    solve_remaining_inner_vars_simple(static, doms, inner_vars) -> Vector{Vector{Bool}}

Simple backtracking search for remaining unfixed inner variables.
Returns all satisfying assignments for the inner variables.
"""
function solve_remaining_inner_vars_simple(static::TNStatic, doms::Vector{DomainMask}, inner_vars::Vector{Int})
    # Find unfixed inner variables
    unfixed_inner = Int[]
    for var_id in inner_vars
        if !is_fixed(doms[var_id])
            push!(unfixed_inner, var_id)
        end
    end
    
    if isempty(unfixed_inner)
        # All inner vars are fixed, extract the assignment
        config = [has1(doms[v]) for v in inner_vars]
        return [config]
    end
    
    # Use simple enumeration for now (can be optimized with DPLL later)
    n_unfixed = length(unfixed_inner)
    satisfying_configs = Vector{Vector{Bool}}()
    
    # Try all 2^n_unfixed assignments
    for assignment in 0:(2^n_unfixed - 1)
        test_doms = copy(doms)
        
        # Apply assignment
        for (i, var_id) in enumerate(unfixed_inner)
            bit = (assignment >> (i-1)) & 0x1
            test_doms[var_id] = bit == 1 ? DM_1 : DM_0
        end
        
        # Propagate
        propagated = propagate(static, test_doms)
        
        # Check for conflicts
        if !any(dm -> dm.bits == 0x00, propagated)
            # No conflict - this is a satisfying assignment
            config = [has1(propagated[v]) for v in inner_vars]
            push!(satisfying_configs, config)
        end
    end
    
    return satisfying_configs
end

"""
    OptimalBranchingCore.branching_table(problem, solver::TNPropagationSolver, variables)

Compute branching table using constraint propagation instead of tensor network contraction.
This is generally more efficient as it avoids tensor construction and einsum operations.
"""
function OptimalBranchingCore.branching_table(
    problem::TNProblem, 
    solver::TNPropagationSolver, 
    variables::Vector{T}
) where T
    region = get_cached_region(problem)
    isnothing(region) && error("No cached region found. Make sure `select_variables` is called first.")
    
    n_boundary = length(region.boundary_vars)
    n_inner = length(region.inner_vars)
    n_total = n_boundary + n_inner
    @assert n_total == length(variables)
    
    # Handle special case: no variables to branch on
    if n_total == 0
        return BranchingTable(0, [[]])
    end
    
    # Separate fixed and free boundary variables
    fixed_boundary_vars, fixed_boundary_vals, free_boundary_vars, _ = 
        separate_fixed_free_boundary(region, problem.doms)
    
    n_free_boundary = length(free_boundary_vars)
    
    # Handle case: all boundary variables are fixed
    if n_free_boundary == 0
        # Just need to solve for inner variables with all boundary vars fixed
        inner_configs = solve_remaining_inner_vars_simple(
            problem.static, problem.doms, region.inner_vars
        )
        
        if isempty(inner_configs)
            return BranchingTable(n_total, Vector{Vector{Int}}())
        end
        
        # Construct full configs (boundary + inner)
        boundary_config = [has1(problem.doms[v]) for v in region.boundary_vars]
        full_configs = [vcat(boundary_config, inner_cfg) for inner_cfg in inner_configs]
        return BranchingTable(n_total, [full_configs])
    end
    
    # Main algorithm: enumerate free boundary configurations
    valid_config_groups = Vector{Vector{Vector{Bool}}}()
    
    for boundary_assignment in 0:(2^n_free_boundary - 1)
        # Create temporary domain with this boundary assignment
        temp_doms = copy(problem.doms)
        
        for (i, var_id) in enumerate(free_boundary_vars)
            bit = (boundary_assignment >> (i-1)) & 0x1
            temp_doms[var_id] = bit == 1 ? DM_1 : DM_0
        end
        
        # Propagate constraints
        propagated = propagate(problem.static, temp_doms)
        
        # Check for conflicts
        if any(dm -> dm.bits == 0x00, propagated)
            continue  # This boundary config is unsatisfiable
        end
        
        # Solve for remaining inner variables
        inner_configs = solve_remaining_inner_vars_simple(
            problem.static, propagated, region.inner_vars
        )
        
        if isempty(inner_configs)
            continue  # No satisfying assignment for inner vars
        end
        
        # Construct full boundary configuration
        boundary_config = [has1(propagated[v]) for v in region.boundary_vars]
        
        # Combine boundary and inner configs
        full_configs = [vcat(boundary_config, inner_cfg) for inner_cfg in inner_configs]
        push!(valid_config_groups, full_configs)
    end
    
    isempty(valid_config_groups) && return BranchingTable(n_total, Vector{Vector{Int}}())
    return BranchingTable(n_total, valid_config_groups)
end
