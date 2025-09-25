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
function OptimalBranchingCore.branching_table(bip::BooleanInferenceProblem, bs::AbstractBranchingStatus, solver::TNContractionSolver, subbip::SubBIP)
	out_vs_num = length(subbip.outside_vs_ind)
	vs_num = length(subbip.vs)
	
	# Early return for edge cases
	out_vs_num == 0 && return BranchingTable(0, [Int[]])
	
	# Precompute internal variable indices
	internal_vs_ind = setdiff(1:vs_num, subbip.outside_vs_ind)
	
	# Create optimized index mappings using BitSet for faster membership testing
	outside_vs_set = BitSet(subbip.outside_vs_ind)
	
	# Precompute position mappings for both boundary and internal variables
	outside_pos = Dict{Int, Int}()
	internal_pos = Dict{Int, Int}()
	for (i, idx) in enumerate(subbip.outside_vs_ind)
		outside_pos[idx] = i
	end
	for (i, idx) in enumerate(internal_vs_ind)
		internal_pos[idx] = i
	end
	
	possible_configurations = Vector{Vector{Bool}}[]
	sizehint!(possible_configurations, min(2^out_vs_num, 1000))  # Cap the hint at reasonable size
	
	# Pre-allocate reusable vectors
	out_index = Vector{Int}(undef, out_vs_num)
	vec = Vector{Any}(undef, vs_num)
	
	# Fill the constant parts of vec (internal variables always get :)
	for var_idx in 1:vs_num
		if !(var_idx in outside_vs_set)
			vec[var_idx] = (:)
		end
	end
	
	# Enumerate all possible configurations of boundary variables (2^out_vs_num possibilities)
	for config_idx in 0:(2^out_vs_num-1)
		# Convert integer config_idx to tensor indices directly
		@inbounds for j in 0:out_vs_num-1
			out_index[j+1] = (config_idx & (1 << j) != 0) ? 2 : 1
		end
		
		# Update only the boundary variable positions in vec
		@inbounds for var_idx in subbip.outside_vs_ind
			vec[var_idx] = out_index[outside_pos[var_idx]]
		end
		
		# Slice by the fixed boundary 
		# Find all positions where the contracted tensor equals Tropical(0.0) (feasible)
		in_indies = findall(==(Tropical(0.0)), subbip.sub_tensors[vec...])
		
		# Skip this boundary configuration if no feasible internal assignments exist
		isempty(in_indies) && continue
		
		# Pre-allocate the configuration vector for this boundary assignment
		pcs = Vector{Vector{Bool}}(undef, length(in_indies))
		
		# Generate all possible configurations for this boundary assignment
		@inbounds for (k, in_index) in pairs(in_indies)
			assignment = Vector{Bool}(undef, vs_num)
			
			# Process boundary variables
			for var_idx in subbip.outside_vs_ind
				assignment[var_idx] = out_index[outside_pos[var_idx]] == 2
			end
			
			# Process internal variables
			for var_idx in internal_vs_ind
				assignment[var_idx] = in_index[internal_pos[var_idx]] == 2
			end
			
			pcs[k] = assignment
		end
		push!(possible_configurations, pcs)
	end
	
	# Return empty branching table if no feasible configurations found, otherwise return the computed table
	return isempty(possible_configurations) ? BranchingTable(0, [Int[]]) : BranchingTable(vs_num, possible_configurations)
end