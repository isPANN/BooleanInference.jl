module BooleanInference

using TropicalNumbers
using SparseArrays
using OptimalBranchingCore
using OptimalBranchingCore: AbstractProblem,select_variables,apply_branch,reduce_problem,_vec2int,optimal_branching_rule,candidate_clauses
using OptimalBranchingCore.BitBasis
using GenericTensorNetworks
using GenericTensorNetworks.OMEinsum
import ProblemReductions
import ProblemReductions: CircuitSAT, Circuit, Factoring, reduceto, SatisfiabilityHard, Satisfiability
using KaHyPar

# using GenericTensorNetworks: ∧, ∨, ¬

# status
export BranchingStatus, initialize_branching_status
# stride
export get_tensor_number,slice_tensor, vec2tensor,indices_to_mask,lluint2vec
# types
export BooleanInferenceProblem,BooleanResultBranchCount,NumOfVertices,NumOfClauses,NumOfDegrees

# interface
export convert_cnf_to_bip,convert_circuit_to_bip,convert_sat_to_bip,solve_boolean_inference_problem,solve_factoring,solve_sat_with_assignments
export solve_boolean_inference_problem_with_contraction,solve_sat_problem_with_contraction,solve_factoring_with_contraction,branch_and_reduce_with_contraction

# reducer
export NoReducer,decide_literal

# selector
export KNeighborSelector,neighboring,k_neighboring,subhg,Smallest2NeighborSelector,KaHyParSelector

# tablesolver
export TNContractionSolver

# legacy extras moved under legacy/

include("status.jl")
include("stride.jl")
include("types.jl")
include("interface.jl")
include("reducer.jl")
include("selector.jl")
include("tablesolver.jl")
include("branch.jl")

## legacy modules parked under src/legacy; not loaded by default

end
