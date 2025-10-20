# Simple contraction result storage
struct RegionContraction
    tensor::AbstractArray{Tropical{Float64}}
    vars::Vector{Int}
end

# Hash key based on region structure (independent of variable assignments)
struct RegionKey
    tensors::Vector{Int}
    inner_vars::Vector{Int}
    boundary_vars::Vector{Int}
end

RegionKey(region::Region) = RegionKey(region.tensors, region.inner_vars, region.boundary_vars)

function Base.hash(k::RegionKey, h::UInt)
    h = hash(k.tensors, h)
    h = hash(k.inner_vars, h)
    h = hash(k.boundary_vars, h)
    return h
end

function Base.:(==)(a::RegionKey, b::RegionKey)
    return a.tensors == b.tensors &&
           a.inner_vars == b.inner_vars &&
           a.boundary_vars == b.boundary_vars
end

# Global caches (assuming only one problem is solved at a time)
const REGION_CONTRACTION_CACHE = Dict{RegionKey, RegionContraction}()
const LAST_REGION = Ref{Union{Nothing, Region}}(nothing)

# ============ Core API Functions ============

function cache_region!(problem::TNProblem, region::Region)
    LAST_REGION[] = region
    return nothing
end

function get_cached_region(problem::TNProblem)
    return LAST_REGION[]
end

function get_cached_contraction(problem::TNProblem, region::Region)
    region_key = RegionKey(region)
    return get(REGION_CONTRACTION_CACHE, region_key, nothing)
end

function cache_contraction!(problem::TNProblem, region::Region, tensor::AbstractArray{Tropical{Float64}}, vars::Vector{Int})
    region_key = RegionKey(region)
    contraction = RegionContraction(tensor, vars)
    REGION_CONTRACTION_CACHE[region_key] = contraction
    return contraction
end

function clear_region_cache!(problem::TNProblem)
    # Only clear temporary region reference, keep persistent contractions
    LAST_REGION[] = nothing
    return nothing
end

function clear_all_region_caches!()
    empty!(REGION_CONTRACTION_CACHE)
    LAST_REGION[] = nothing
    return nothing
end

function set_last_region!(problem::TNProblem, region_id::Int)
    problem.ws.last_region_id = region_id
    return nothing
end
