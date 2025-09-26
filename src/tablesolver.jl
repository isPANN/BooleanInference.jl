"""
    TNContractionSolver <: AbstractTableSolver

A table solver that uses tensor network contraction to compute branching tables.
This solver enumerates all possible configurations of boundary variables and uses
tensor network contraction to determine feasible internal variable assignments.
"""
struct TNContractionSolver <: AbstractTableSolver end

"""
    branching_table(bip::BooleanInferenceProblem, bs::AbstractBranchingStatus, solver::TNContractionSolver, subbip::SubBIP)

Compute a branching table for the given subproblem using tensor network contraction.

# Arguments
- `bip`: The boolean inference problem
- `bs`: Current branching status
- `solver`: The TNContractionSolver instance
- `subbip`: Subproblem containing selected variables and their associated tensor data

# Returns
- `BranchingTable`: Contains all possible branching configurations for the subproblem

# Algorithm
1. Enumerate all possible configurations of boundary variables (variables connected to external clauses)
2. For each boundary configuration, use tensor network contraction to find feasible internal variable assignments
3. Group configurations by boundary variable assignments
4. Return a branching table containing all valid branching decisions

The function uses the following key concepts:
- **Boundary variables**: Variables that appear in clauses outside the current subproblem
- **Internal variables**: Variables that only appear in clauses within the current subproblem  
- **Tensor contraction**: Efficiently computes the feasibility of variable assignments using tropical arithmetic
"""
function OptimalBranchingCore.branching_table(
    bip::BooleanInferenceProblem,
    bs::AbstractBranchingStatus,
    solver::TNContractionSolver,
    subbip::SubBIP
)
    out_vs_num = length(subbip.outside_vs_ind)
    vs_num = length(subbip.vs)

    # pick packed integer type for bit_length = vs_num
    INT = _best_uint(Val(vs_num))  # returns a Type like UInt64 or LongLongUInt{N}

    # Early return (empty): table::Vector{Vector{INT}} with one empty group
    if out_vs_num == 0
        empty_group = Vector{INT}()                        # Vector{INT}()
        empty_table = Vector{Vector{INT}}([empty_group])   # Vector{Vector{INT}}
        return BranchingTable{INT}(0, empty_table)
    end

    outside = subbip.outside_vs_ind

    # pos_outside[v] = position of v in `outside` (0 if not boundary)
    pos_outside = zeros(Int, vs_num)
    @inbounds for (i, v) in pairs(outside)
        pos_outside[v] = i
    end

    # internal list by one linear scan
    internal_len = vs_num - out_vs_num
    internal_vs_ind = Vector{Int}(undef, internal_len)
    k = 0
    @inbounds for v in 1:vs_num
        if pos_outside[v] == 0
            k += 1
            internal_vs_ind[k] = v
        end
    end

    # pos_internal[v] = position of v in internal list (0 if not internal)
    pos_internal = zeros(Int, vs_num)
    @inbounds for (i, v) in pairs(internal_vs_ind)
        pos_internal[v] = i
    end

    # ---- preallocs ----
    possible_configurations = Vector{Vector{INT}}()  # groups; each group is Vector{INT}
    local_hint = (out_vs_num <= 60) ? min(1 << out_vs_num, 1000) : 1000
    sizehint!(possible_configurations, local_hint)

    # out_index[j] ∈ {1,2} for boundary variables, maps bit -> tensor index
    out_index = Vector{Int}(undef, out_vs_num)

    # vec holds slicing indices for subbip.sub_tensors[vec...]
    # initialize all to : ; only boundary dims will be overwritten to 1/2 per config
    vec = Vector{Union{Colon,Int}}(undef, vs_num)
    fill!(vec, Colon())

    # ---- enumerate all boundary configurations ----
    @inbounds for config_idx in 0:(2^out_vs_num-1)
        # integer -> tensor indices {1,2}
        @inbounds for j in 0:out_vs_num-1
            out_index[j+1] = (config_idx & (1 << j)) != 0 ? 2 : 1
        end

        # update boundary slice indices
        @inbounds for v in outside
            vec[v] = out_index[pos_outside[v]]
        end

        # feasible internal assignments => Tropical(0.0)
        # Use view to avoid copying the entire sliced array
        tensor_view = @view subbip.sub_tensors[vec...]
        in_indies = findall(==(Tropical(0.0)), tensor_view)
        isempty(in_indies) && continue

        # pack each feasible full assignment to an INT
        pcs = Vector{INT}(undef, length(in_indies))
        @inbounds for (k2, in_index) in pairs(in_indies)
            x = zero(INT)

            # boundary bits
            @inbounds for v in outside
                if out_index[pos_outside[v]] == 2
                    x |= (one(INT) << (v - 1))
                end
            end
            # internal bits (in_index is over internal dims)
            @inbounds for v in internal_vs_ind
                if in_index[pos_internal[v]] == 2
                    x |= (one(INT) << (v - 1))
                end
            end

            pcs[k2] = x
        end

        push!(possible_configurations, pcs)
    end

    if isempty(possible_configurations)
        empty_group = Vector{INT}()
        empty_table = Vector{Vector{INT}}([empty_group])
        return BranchingTable{INT}(0, empty_table)
    else
        return BranchingTable{INT}(vs_num, possible_configurations)
    end
end

# choose smallest unsigned for bit_length
@generated function _best_uint(nbits::Val{N}) where {N}
    if N <= 8
        :(UInt8)
    elseif N <= 16
        :(UInt16)
    elseif N <= 32
        :(UInt32)
    elseif N <= 64
        :(UInt64)
    elseif N <= 128
        :(UInt128)
    else
        :(LongLongUInt{N})
    end
end
