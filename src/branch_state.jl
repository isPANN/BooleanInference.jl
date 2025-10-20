# Branching state structures (loaded after problems.jl, before branch.jl)

struct BranchingCandidate
    focus_var::Int
    region::Union{Nothing,Region}
    variables::Vector{Int}
    result::OptimalBranchingCore.OptimalBranchingResult
    versions::Vector{UInt32}  # version snapshot of the variables
end

mutable struct GammaQueueState
    # BinaryMinHeap stores (-gamma, var_id) tuples for better performance
    # Negative gamma to implement max-heap (higher gamma = higher priority)
    queue::BinaryMinHeap{Tuple{Float64,Int}}
    candidates::Dict{Int,BranchingCandidate}  # var_id -> candidate
    var_regions::Vector{Vector{Int}}  # var_id -> list of candidate ids depending on this var
    pending_candidates::Set{Int}  # candidate var_ids waiting for refresh
    deleted_vars::Set{Int}  # vars marked for lazy deletion
end

function GammaQueueState(var_num::Int)
    GammaQueueState(
        BinaryMinHeap{Tuple{Float64,Int}}(),
        Dict{Int,BranchingCandidate}(),
        [Int[] for _ in 1:var_num],
        Set{Int}(),
        Set{Int}(),
    )
end
