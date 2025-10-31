using BooleanInference
using BooleanInference: TNProblem, TNContractionSolver, LeastOccurrenceSelector, NumUnfixedVars
using BooleanInference: BranchingStrategy, NoReducer
using OptimalBranchingCore: branch_and_reduce
using ProblemReductions: Factoring, reduceto, CircuitSAT
using GenericTensorNetworks
using TropicalNumbers: Tropical
using Printf

"""
    benchmark_factoring(n, m, N; num_runs=3)

Benchmark factoring problem and return statistics.
"""
function benchmark_factoring(n, m, N; num_runs=3)
    times = Float64[]
    stats_list = []

    for i in 1:num_runs
        # Setup problem
        fproblem = Factoring(n, m, N)
        circuit_sat = reduceto(CircuitSAT, fproblem)
        problem = CircuitSAT(circuit_sat.circuit.circuit; use_constraints=true)

        tn = GenericTensorNetwork(problem)
        tn_static = BooleanInference.setup_from_tensor_network(tn)
        tn_problem = TNProblem(tn_static)

        br_strategy = BranchingStrategy(
            table_solver = TNContractionSolver(),
            selector = LeastOccurrenceSelector(2, 5),
            measure = NumUnfixedVars()
        )

        # Run and time
        GC.gc()  # Force garbage collection before timing
        t_start = time()
        res = branch_and_reduce_optimized(tn_problem, br_strategy, NoReducer(), Tropical{Float64}; show_progress=false)
        t_elapsed = time() - t_start

        push!(times, t_elapsed)
        push!(stats_list, BooleanInference.get_branching_stats(tn_problem))
    end

    return (
        mean_time = sum(times) / length(times),
        median_time = sort(times)[div(length(times)+1, 2)],
        min_time = minimum(times),
        max_time = maximum(times),
        stats = stats_list[1]  # Return first run's stats
    )
end

"""
    run_benchmark_suite()

Run a suite of benchmark problems and print results.
"""
function run_benchmark_suite()
    println("=" ^ 80)
    println("Branch Optimization Benchmark Suite")
    println("=" ^ 80)

    # Test cases: (n, m, N)
    test_cases = [
        (4, 4, 143),
        (6, 6, 1147),
        (8, 8, 40759),
        (10, 10, 614891),
        (12, 12, 10371761),
        (14, 14, 183974111),
    ]

    println("\nRunning benchmarks (3 runs each, reporting median)...\n")

    results = []
    for (n, m, N) in test_cases
        print("Testing Factoring($n, $m, $N)... ")
        flush(stdout)

        try
            result = benchmark_factoring(n, m, N; num_runs=3)
            push!(results, (n=n, m=m, N=N, result=result))
            @printf("%.3f s\n", result.median_time)
        catch e
            println("FAILED: $e")
        end
    end

    # Print detailed results
    println("\n" * "=" ^ 80)
    println("Detailed Results")
    println("=" ^ 80)

    for r in results
        println("\nFactoring($(r.n), $(r.m), $(r.N)):")
        @printf("  Median time: %.3f s\n", r.result.median_time)
        @printf("  Min time: %.3f s\n", r.result.min_time)
        @printf("  Max time: %.3f s\n", r.result.max_time)
        println("  Total branches: $(r.result.stats.total_branches)")
        println("  Total subproblems: $(r.result.stats.total_subproblems)")
        println("  Max depth: $(r.result.stats.max_depth)")
        @printf("  Avg branching factor: %.2f\n", r.result.stats.avg_branching_factor)
    end

    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark_suite()
end
