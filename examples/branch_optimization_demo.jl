"""
Demo: Branch Optimization Strategies

This example demonstrates how to use various branch optimization heuristics
to reduce the number of branches explored during solving.
"""

using BooleanInference
using GenericTensorNetworks
using ProblemReductions

# Load optimization modules (if they're included in main module)
# include("../src/branch_optimized.jl")

"""
Example 1: Basic usage with look-ahead pruning
"""
function demo_lookahead()
    println("\n=== Demo 1: Look-ahead Branch Pruning ===")

    # Create a simple SAT problem
    cnf = CNF([
        [1, 2, 3],
        [-1, 2],
        [-2, 3],
        [1, -3]
    ])

    problem = tnproblem(cnf)

    # Solve WITHOUT optimization
    println("\n1. Baseline (no optimization):")
    config_baseline = BranchingStrategy(
        selector=LeastOccurrenceSelector(2),
        measure=NumUnfixedVars(),
        table_solver=TNContractionSolver(),
        set_cover_solver=GreedyMerge()
    )

    ws_baseline = Workspace()
    problem_baseline = TNProblem(problem.static, copy(problem.doms), problem.n_unfixed, ws_baseline)

    result_baseline = solve(BooleanInference.CNF, problem_baseline, config_baseline)
    println("   Branches: $(ws_baseline.total_branches)")
    println("   Subproblems: $(ws_baseline.total_subproblems)")
    println("   Avg branching factor: $(ws_baseline.total_subproblems / max(1, ws_baseline.total_branches))")

    # Solve WITH look-ahead optimization
    println("\n2. With Look-ahead (max_branches=8):")
    ws_optimized = Workspace()
    problem_optimized = TNProblem(problem.static, copy(problem.doms), problem.n_unfixed, ws_optimized)

    opt_config = BranchOptimizationConfig(
        enable_lookahead=true,
        max_branches=8,
        min_propagation_ratio=0.05
    )

    result_optimized = branch_and_reduce_optimized(
        problem_optimized,
        config_baseline,
        NoReducer(),
        CountingTropical{Float64,Tropical{Float64}},
        opt_config
    )

    println("   Branches: $(ws_optimized.total_branches)")
    println("   Subproblems: $(ws_optimized.total_subproblems)")
    println("   Avg branching factor: $(ws_optimized.total_subproblems / max(1, ws_optimized.total_branches))")

    reduction = 100 * (1 - ws_optimized.total_branches / max(1, ws_baseline.total_branches))
    println("\n   Branch reduction: $(round(reduction, digits=1))%")

    return result_baseline, result_optimized
end

"""
Example 2: Comparison of different max_branches values
"""
function demo_parameter_sweep()
    println("\n=== Demo 2: Parameter Sweep (max_branches) ===")

    # Create a moderately complex problem
    cnf = CNF([
        [1, 2, 3], [-1, 2, 4], [-2, 3, 5],
        [1, -3, -5], [-1, -4, 5], [2, 4, -5],
        [-2, -3, 4], [1, 3, -4]
    ])

    problem = tnproblem(cnf)
    config = BranchingStrategy(
        selector=LeastOccurrenceSelector(2),
        measure=NumUnfixedVars(),
        table_solver=TNContractionSolver(),
        set_cover_solver=GreedyMerge()
    )

    println("\nTesting different max_branches values:")
    println("max_branches | total_branches | total_subproblems | avg_factor | reduction")
    println("-" ^ 75)

    baseline_branches = nothing

    for max_branches in [32, 16, 8, 4, 2]
        ws = Workspace()
        test_problem = TNProblem(problem.static, copy(problem.doms), problem.n_unfixed, ws)

        opt_config = BranchOptimizationConfig(
            enable_lookahead=true,
            max_branches=max_branches,
            min_propagation_ratio=0.05
        )

        try
            result = branch_and_reduce_optimized(
                test_problem, config, NoReducer(),
                CountingTropical{Float64,Tropical{Float64}},
                opt_config
            )

            avg_factor = ws.total_subproblems / max(1, ws.total_branches)

            if isnothing(baseline_branches)
                baseline_branches = ws.total_branches
                reduction_str = "baseline"
            else
                reduction = 100 * (1 - ws.total_branches / baseline_branches)
                reduction_str = "$(round(reduction, digits=1))%"
            end

            @printf("%-12d | %-14d | %-17d | %-10.2f | %s\n",
                    max_branches, ws.total_branches, ws.total_subproblems,
                    avg_factor, reduction_str)
        catch e
            println("Error with max_branches=$max_branches: $e")
        end
    end
end

"""
Example 3: Analyzing branch quality
"""
function demo_branch_quality_analysis()
    println("\n=== Demo 3: Branch Quality Analysis ===")

    # Simple problem to analyze in detail
    cnf = CNF([
        [1, 2], [-1, 3], [2, -3], [-2, 3]
    ])

    problem = tnproblem(cnf)

    println("\nAnalyzing branch propagation effects:")

    # Select variables
    config = BranchingStrategy(
        selector=LeastOccurrenceSelector(1),
        measure=NumUnfixedVars(),
        table_solver=TNContractionSolver(),
        set_cover_solver=GreedyMerge()
    )

    variables = OptimalBranchingCore.select_variables(problem, config.measure, config.selector)
    println("Selected variables: $variables")

    # Build branching table
    tbl = OptimalBranchingCore.branching_table(problem, config.table_solver, variables)
    println("Total configurations: $(sum(length(group) for group in tbl.table))")

    # Analyze each configuration
    println("\nConfiguration analysis:")
    println("Config | Fixed vars | Has conflict | Quality")
    println("-" ^ 50)

    for (group_idx, group) in enumerate(tbl.table)
        for (config_idx, config_vec) in enumerate(group)
            mask, val = clause_to_bits(config_vec, length(variables))
            clause = OptimalBranchingCore.Clause(mask, val)

            n_fixed, is_conflict = quick_evaluate_clause(problem, clause, variables)

            quality = is_conflict ? "BAD" : (n_fixed > 0 ? "GOOD" : "OK")
            conflict_str = is_conflict ? "YES" : "NO"

            println("$group_idx.$config_idx   | $n_fixed         | $conflict_str         | $quality")
        end
    end
end

"""
Example 4: Comparison with different selectors
"""
function demo_selector_comparison()
    println("\n=== Demo 4: Variable Selector Comparison ===")

    # Medium complexity problem
    cnf = CNF([
        [1, 2, 3], [-1, 2, 4], [-2, 3, 5], [1, -3, -5],
        [-1, -4, 5], [2, 4, -5], [-2, -3, 4], [1, 3, -4],
        [1, 4, 5], [-2, -4, -5]
    ])

    problem = tnproblem(cnf)

    println("\nComparing different variable selection strategies:")
    println("Selector | k | Branches | Subproblems | Avg Factor")
    println("-" ^ 60)

    opt_config = BranchOptimizationConfig(enable_lookahead=true, max_branches=12)

    for k in [1, 2, 3]
        selector = LeastOccurrenceSelector(k, 2)
        config = BranchingStrategy(
            selector=selector,
            measure=NumUnfixedVars(),
            table_solver=TNContractionSolver(),
            set_cover_solver=GreedyMerge()
        )

        ws = Workspace()
        test_problem = TNProblem(problem.static, copy(problem.doms), problem.n_unfixed, ws)

        try
            result = branch_and_reduce_optimized(
                test_problem, config, NoReducer(),
                CountingTropical{Float64,Tropical{Float64}},
                opt_config
            )

            avg_factor = ws.total_subproblems / max(1, ws.total_branches)
            @printf("%-8s | %d | %-8d | %-11d | %.2f\n",
                    "LeastOcc", k, ws.total_branches, ws.total_subproblems, avg_factor)
        catch e
            println("Error with k=$k: $e")
        end
    end
end

"""
Run all demos
"""
function run_all_demos()
    println("=" ^ 80)
    println("Branch Optimization Strategy Demonstrations")
    println("=" ^ 80)

    try
        demo_lookahead()
    catch e
        println("Demo 1 failed: $e")
    end

    try
        demo_parameter_sweep()
    catch e
        println("Demo 2 failed: $e")
    end

    try
        demo_branch_quality_analysis()
    catch e
        println("Demo 3 failed: $e")
    end

    try
        demo_selector_comparison()
    catch e
        println("Demo 4 failed: $e")
    end

    println("\n" * "=" ^ 80)
    println("Demos complete!")
    println("=" ^ 80)
end

# Uncomment to run when loaded
# run_all_demos()
