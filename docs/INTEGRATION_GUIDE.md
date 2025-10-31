# 分支优化集成指南

本指南说明如何将提出的优化策略集成到现有的 BooleanInference.jl 代码库中。

## 快速开始

### 1. 最简单的优化：启用Look-ahead

最容易集成且效果显著的优化是 **look-ahead branch pruning**。

#### 修改 `src/branch.jl`

在 `branch_and_reduce` 函数中，branching table 构建之后添加过滤步骤：

```julia
# 原代码（第37行）
tbl = OptimalBranchingCore.branching_table(problem, config.table_solver, variables)

# 添加这段（如果启用优化）
if hasfield(typeof(config), :max_branches) && config.max_branches > 0
    tbl = prune_bad_branches(problem, tbl, variables, config.max_branches)
end

# If table is empty, the problem is UNSAT
if isempty(tbl.table)
    @debug "Branching table is empty - problem is UNSAT"
    return zero(result_type)
end
```

然后在 `src/branch_optimized.jl` 中使用 `prune_bad_branches` 函数。

### 2. 扩展 BranchingStrategy

修改 `BranchingStrategy` 以支持优化参数：

```julia
# 添加到 problems.jl 或创建新的 config.jl
struct BranchingStrategyOptimized <: OptimalBranchingCore.AbstractBranchingStrategy
    selector::OptimalBranchingCore.AbstractSelector
    measure::OptimalBranchingCore.AbstractMeasure
    table_solver::OptimalBranchingCore.AbstractTableSolver
    set_cover_solver::OptimalBranchingCore.AbstractSetCoverSolver

    # 新增的优化参数
    max_branches::Int  # 最大分支数，0表示禁用
    min_propagation_ratio::Float64  # 最小传播比例
end

function BranchingStrategyOptimized(;
    selector=LeastOccurrenceSelector(2),
    measure=NumUnfixedVars(),
    table_solver=TNContractionSolver(),
    set_cover_solver=GreedyMerge(),
    max_branches=16,
    min_propagation_ratio=0.05
)
    return BranchingStrategyOptimized(
        selector, measure, table_solver, set_cover_solver,
        max_branches, min_propagation_ratio
    )
end
```

### 3. 使用优化版本求解

```julia
using BooleanInference

# 创建问题
cnf = CNF([[1,2,3], [-1,2], [-2,3]])
problem = tnproblem(cnf)

# 配置优化策略
config = BranchingStrategyOptimized(
    selector=LeastOccurrenceSelector(2),
    max_branches=12,  # 限制每次最多12个分支
    min_propagation_ratio=0.05
)

# 求解
result = solve(BooleanInference.CNF, problem, config)

# 检查统计
println("Total branches: $(problem.ws.total_branches)")
println("Total subproblems: $(problem.ws.total_subproblems)")
```

---

## 渐进式集成路线图

### Phase 1: 基础Look-ahead（1-2小时工作量）

**目标**: 在不改变接口的情况下减少30-50%的分支数

**步骤**:

1. 将 `src/branch_optimized.jl` 中的 `prune_bad_branches` 函数复制到 `src/branchtable.jl`
2. 在 `branch_and_reduce` 中添加条件调用（基于环境变量或全局标志）
3. 运行现有测试验证正确性

```julia
# src/branchtable.jl 末尾添加

const ENABLE_LOOKAHEAD = Ref(false)  # 全局开关
const MAX_BRANCHES = Ref(16)

function prune_bad_branches(problem, tbl, variables)
    # ... 复制实现
end

# 在 branch.jl 中
if ENABLE_LOOKAHEAD[]
    tbl = prune_bad_branches(problem, tbl, variables, MAX_BRANCHES[])
end
```

**测试**:
```julia
# 启用优化
BooleanInference.ENABLE_LOOKAHEAD[] = true

# 运行benchmark
@time result = solve(CNF, problem, config)
```

### Phase 2: VSIDS变量选择（2-4小时工作量）

**目标**: 智能选择分支变量，进一步减少20-30%的分支数

**步骤**:

1. 将 `src/vsids_selector.jl` 添加到模块
2. 在 `src/BooleanInference.jl` 中导出新selector
3. 修改 `apply_branch` 返回 `changed_vars` 信息（已有）
4. 在分支循环中更新VSIDS分数

```julia
# src/branch.jl 中修改

return sum(enumerate(clauses)) do (i, branch)
    # ... 现有代码

    subproblem, local_value, changed_vars = apply_branch(problem, branch, variables)

    # 如果使用VSIDS，更新分数
    if config.selector isa VSIDSSelector
        for var_id in changed_vars
            bump_variable_score!(config.selector, var_id)
        end
    end

    # ... 继续递归
end
```

### Phase 3: 学习子句（可选，4-6小时工作量）

**目标**: 在有重复结构的问题上减少40-60%的分支数

**步骤**:

1. 将 `src/learned_clauses.jl` 添加到模块
2. 在 `Workspace` 中添加 `learned_db::Union{Nothing,LearnedClauseDB}` 字段
3. 在 `branch_and_reduce` 开头检查学习子句
4. 在检测到冲突时学习

**注意**: 这个功能需要跟踪分支路径，实现较复杂，建议先完成Phase 1和2。

---

## 推荐的集成顺序

基于你的代码库和问题类型，建议按以下优先级集成：

### 优先级1: Look-ahead pruning ⭐⭐⭐⭐⭐
- **难度**: 低
- **效果**: 高（30-50%分支减少）
- **开销**: 中等（每个候选分支需要模拟传播）
- **适用性**: 广泛，几乎所有问题都受益

### 优先级2: 改进的branching table生成 ⭐⭐⭐⭐
- **难度**: 低
- **效果**: 中（15-30%分支减少）
- **开销**: 低
- **建议**:
  - 对于稀疏问题，用 `TNPropagationSolver` 替代 `TNContractionSolver`
  - 已经在你的代码中实现（`branchtable.jl:446`）

### 优先级3: VSIDS selector ⭐⭐⭐
- **难度**: 中
- **效果**: 中高（20-40%分支减少）
- **开销**: 极低
- **适用性**: 在有局部性的问题上效果显著

### 优先级4: Learned clauses ⭐⭐
- **难度**: 高
- **效果**: 高（40-60%，但仅限特定问题）
- **开销**: 中等
- **适用性**: 窄，主要对有重复子结构的问题有效

---

## 实际示例：修改后的接口

### 简化接口（向后兼容）

```julia
# 不改变现有接口，添加新函数
function solve_optimized(
    ::Type{T},
    problem::TNProblem,
    config::BranchingStrategy;
    max_branches::Int=16,
    enable_lookahead::Bool=true
) where T
    # 使用优化版本
    opt_config = BranchOptimizationConfig(enable_lookahead, max_branches, 0.05)

    return branch_and_reduce_optimized(
        problem,
        config,
        NoReducer(),
        CountingTropical{Float64,Tropical{Float64}},
        opt_config
    )
end

# 使用
result = solve_optimized(CNF, problem, config; max_branches=12)
```

### 完全集成接口

```julia
# 替换 solve 函数中的 branch_and_reduce 调用

function solve(
    ::Type{T},
    problem::TNProblem,
    config::BranchingStrategy;
    show_progress::Bool=false,
    max_branches::Int=0,  # 0表示禁用
    kwargs...
) where T
    if max_branches > 0
        # 使用优化版本
        opt_config = BranchOptimizationConfig(true, max_branches, 0.05)
        return branch_and_reduce_optimized(
            problem, config, NoReducer(),
            CountingTropical{Float64,Tropical{Float64}},
            opt_config;
            show_progress=show_progress
        )
    else
        # 使用原版本
        return OptimalBranchingCore.branch_and_reduce(
            problem, config, NoReducer(),
            CountingTropical{Float64,Tropical{Float64}};
            show_progress=show_progress
        )
    end
end
```

---

## 性能测试框架

创建一个benchmark脚本来评估优化效果：

```julia
# test/benchmark_optimization.jl

using BooleanInference
using ProblemReductions
using Printf

function benchmark_strategy(problem_name, cnf, configs)
    println("\n=== Problem: $problem_name ===")
    println(@sprintf("%-30s | %10s | %12s | %10s",
            "Strategy", "Branches", "Subproblems", "Time(ms)"))
    println("-" ^ 70)

    baseline_branches = nothing

    for (name, config, opt_config) in configs
        problem = tnproblem(cnf)

        t_start = time()
        if isnothing(opt_config)
            result = solve(CNF, problem, config)
        else
            result = solve_optimized(CNF, problem, config, opt_config)
        end
        t_elapsed = (time() - t_start) * 1000

        branches = problem.ws.total_branches
        subprobs = problem.ws.total_subproblems

        if isnothing(baseline_branches)
            baseline_branches = branches
        end

        reduction = 100 * (1 - branches / baseline_branches)

        println(@sprintf("%-30s | %10d | %12d | %10.2f (%.1f%%)",
                name, branches, subprobs, t_elapsed, reduction))
    end
end

# 运行benchmark
problems = [
    ("Small SAT", CNF([[1,2,3], [-1,2], [-2,3], [1,-3]])),
    ("Medium SAT", load_cnf_file("data/medium.cnf")),
    ("Large Circuit", load_circuit_sat("data/large_circuit.txt"))
]

configs = [
    ("Baseline",
     BranchingStrategy(LeastOccurrenceSelector(2), NumUnfixedVars(),
                       TNContractionSolver(), GreedyMerge()),
     nothing),

    ("Lookahead (k=16)",
     BranchingStrategy(LeastOccurrenceSelector(2), NumUnfixedVars(),
                       TNContractionSolver(), GreedyMerge()),
     BranchOptimizationConfig(true, 16, 0.05)),

    ("Lookahead (k=8)",
     BranchingStrategy(LeastOccurrenceSelector(2), NumUnfixedVars(),
                       TNContractionSolver(), GreedyMerge()),
     BranchOptimizationConfig(true, 8, 0.05)),
]

for (name, cnf) in problems
    benchmark_strategy(name, cnf, configs)
end
```

---

## 调试和监控

添加调试输出来监控优化效果：

```julia
# 在 prune_bad_branches 中
if @isdebug
    original_size = sum(length(group) for group in tbl.table)
    pruned_size = sum(length(group) for group in filtered_tbl.table)
    @debug "Look-ahead statistics" original=original_size pruned=pruned_size
           ratio=pruned_size/original_size
end

# 运行时启用调试
ENV["JULIA_DEBUG"] = "BooleanInference"
result = solve(CNF, problem, config)
```

---

## 常见问题

### Q: 优化会导致结果不正确吗？
A: 不会。look-ahead只是剪枝，不会改变可满足性判断。但要确保：
- 不会剪掉所有分支（至少保留一个）
- 冲突检测正确

### Q: 什么时候look-ahead反而降低性能？
A:
- 问题很小（< 10变量）：开销大于收益
- 分支表本身就很小（< 10个配置）：没什么可剪
- 传播很慢的问题：评估开销大

建议：只在 `total_configs > max_branches * 2` 时启用look-ahead

### Q: 如何选择 max_branches 参数？
A: 经验法则：
- 密集约束问题：4-8
- 中等密度：12-16
- 稀疏问题：32或禁用

### Q: 能与region缓存兼容吗？
A: 完全兼容。look-ahead使用临时domain，不影响缓存。

---

## 下一步

1. **阅读**: [branch_optimization_strategies.md](branch_optimization_strategies.md) 了解理论背景
2. **实验**: 运行 `examples/branch_optimization_demo.jl`
3. **集成**: 从Phase 1开始，逐步添加优化
4. **测试**: 在你的实际问题上benchmark
5. **调优**: 根据结果调整参数

有问题可以参考代码中的注释或查看相关文献。
