# 分支优化策略指南

本文档总结了减少不必要分支数的启发式优化策略。

## 核心思想

减少分支数的关键在于：
1. **选择好的分支变量** - 选择能最大化约束传播的变量
2. **过滤劣质分支** - 提前评估分支质量，剪枝无用分支
3. **学习失败模式** - 记住导致冲突的路径，避免重复探索
4. **动态调整策略** - 根据搜索历史调整变量选择

## 优化策略总览

### 1. VSIDS变量选择（Variable State Independent Decaying Sum）

**原理**: 优先选择最近参与冲突/传播的变量，因为这些变量更可能快速减少搜索空间。

**实现**: `src/vsids_selector.jl`

```julia
# 使用示例
selector = VSIDSSelector(2; max_tensors=2, decay=0.95)
config = BranchingStrategy(
    selector=selector,
    measure=NumUnfixedVars(),
    table_solver=TNContractionSolver(),
    set_cover_solver=GreedyMerge()
)
```

**关键机制**:
- 每次变量参与传播时，增加其活跃度分数
- 定期衰减所有分数，保持对最近活动的偏好
- 选择高分变量作为分支点

**预期效果**: 减少20-40%的分支数，特别是在有结构的SAT问题上。

---

### 2. Look-ahead分支评估

**原理**: 在真正分支前，模拟每个候选分支的传播效果，选择能固定最多变量的分支。

**实现**: `src/lookahead.jl`

```julia
# 配置look-ahead
lookahead_config = LookaheadConfig(
    enabled=true,
    max_lookahead_depth=1,
    min_fixed_ratio=0.1,  # 至少固定10%的变量
    max_branches=16        # 最多保留16个最优分支
)

# 在branching_table后应用过滤
tbl = branching_table(problem, solver, variables)
filtered_tbl = filter_branching_table(problem, tbl, variables, lookahead_config)
```

**关键机制**:
- 对每个候选分支，模拟应用+传播
- 计算传播后新固定的变量数
- 保留top-k个传播效果最好的分支
- 丢弃立即导致冲突的分支

**预期效果**:
- 在密集约束问题上减少30-50%的分支数
- 轻微增加每个分支的计算开销（但总体更快）

---

### 3. 学习子句（Learned Clauses）

**原理**: 当一个分支路径导致UNSAT时，提取导致冲突的变量赋值模式，避免将来再次尝试类似路径。

**实现**: `src/learned_clauses.jl`

```julia
# 创建学习子句数据库
learned_db = LearnedClauseDB(max_size=1000)

# 在检测到冲突时学习
if local_value == 0
    learn_conflict!(learned_db, problem, current_branch_path)
    return zero(result_type)
end

# 在分支前检查已知冲突
if check_learned_conflicts(learned_db, problem.doms)
    return zero(result_type)  # 提前剪枝
end
```

**关键机制**:
- 记录导致UNSAT的变量赋值组合
- 在新分支前检查是否匹配已知冲突模式
- 定期清理低活跃度的学习子句

**预期效果**:
- 在有重复子结构的问题上减少40-60%的分支数
- 在随机SAT问题上效果较弱（约10-20%改进）

---

### 4. 改进的变量选择启发式

除了VSIDS，还可以尝试以下启发式：

#### 4.1 最多约束优先（Most Constrained First）

```julia
# 选择domain最小（最受约束）的变量
function select_most_constrained(problem)
    unfixed_vars = get_unfixed_vars(problem.doms)
    return argmin(unfixed_vars) do var_id
        count_feasible_values(problem, var_id)
    end
end
```

#### 4.2 最大传播潜力（Maximum Propagation Potential）

```julia
# 选择预期能触发最多传播的变量
function select_max_propagation(problem)
    unfixed_vars = get_unfixed_vars(problem.doms)
    return argmax(unfixed_vars) do var_id
        estimate_propagation_impact(problem, var_id)
    end
end
```

#### 4.3 最小度数优先（Minimum Width）

```julia
# 选择参与最少约束的变量（你目前的策略）
# 优点：region较小，branching table构建快
# 缺点：可能选到不重要的变量
function select_min_degree(problem)
    unfixed_vars = get_unfixed_vars(problem.doms)
    return argmin(u -> length(problem.static.v2t[u]), unfixed_vars)
end
```

---

## 综合策略：多层次剪枝

最佳实践是组合多种策略：

```julia
"""
综合优化的branch_and_reduce实现
"""
function optimized_branch_and_reduce(
    problem::TNProblem,
    config::BranchingStrategy,
    reducer::AbstractReducer,
    result_type::Type{TR};
    learned_db::Union{Nothing,LearnedClauseDB}=nothing,
    lookahead_config::LookaheadConfig=LookaheadConfig(),
    branch_path::Vector{Tuple{Int,Bool}}=Tuple{Int,Bool}[]
) where TR

    # Level 1: 检查学习子句（最快）
    if !isnothing(learned_db) && check_learned_conflicts(learned_db, problem.doms)
        return zero(result_type)
    end

    # Level 2: 检查是否已解决
    if is_solved(problem)
        cache_branch_solution!(problem)
        return one(result_type)
    end

    # Level 3: 智能变量选择（VSIDS或其他启发式）
    variables = select_variables(problem, config.measure, config.selector)

    # Level 4: 构建分支表
    tbl = branching_table(problem, config.table_solver, variables)
    isempty(tbl.table) && return zero(result_type)

    # Level 5: Look-ahead过滤劣质分支
    if lookahead_config.enabled
        tbl = filter_branching_table(problem, tbl, variables, lookahead_config)
        isempty(tbl.table) && return zero(result_type)
    end

    # Level 6: 计算最优分支规则
    result = optimal_branching_rule(tbl, variables, problem, config.measure, config.set_cover_solver)
    clauses = get_clauses(result)

    # Level 7: 递归求解各分支
    return sum(enumerate(clauses)) do (i, branch)
        subproblem, local_value, changed_vars = apply_branch(problem, branch, variables)

        # 记录冲突并学习
        if local_value == 0
            !isnothing(learned_db) && learn_conflict!(learned_db, problem, branch_path)
            return zero(result_type)
        end

        # 更新VSIDS分数
        if config.selector isa VSIDSSelector
            bump_conflicting_variables!(config.selector, changed_vars)
        end

        # 递归
        new_branch_path = [branch_path..., extract_assignments(branch, variables)...]
        optimized_branch_and_reduce(
            subproblem, config, reducer, result_type;
            learned_db=learned_db,
            lookahead_config=lookahead_config,
            branch_path=new_branch_path
        )
    end
end
```

---

## 性能调优建议

### 针对不同问题类型的配置

#### 1. 密集约束问题（如电路SAT）
```julia
# 启用look-ahead，限制分支数
config = LookaheadConfig(
    enabled=true,
    max_branches=8,
    min_fixed_ratio=0.2
)
selector = VSIDSSelector(2; max_tensors=3)
```

#### 2. 稀疏随机SAT
```julia
# 轻量级look-ahead，依赖学习子句
config = LookaheadConfig(
    enabled=true,
    max_branches=32,
    min_fixed_ratio=0.05
)
learned_db = LearnedClauseDB(max_size=2000)
```

#### 3. 结构化问题（如因数分解）
```julia
# 重度look-ahead，小region
config = LookaheadConfig(
    enabled=true,
    max_branches=4,
    min_fixed_ratio=0.3
)
selector = VSIDSSelector(1; max_tensors=2)  # 更小的k
```

---

## 实验建议

### 基准测试框架

```julia
function benchmark_branching_strategies()
    problems = load_test_problems()
    strategies = [
        ("Baseline", LeastOccurrenceSelector(2)),
        ("VSIDS", VSIDSSelector(2)),
        ("VSIDS+Lookahead", VSIDSSelector(2), LookaheadConfig(enabled=true)),
        ("VSIDS+Lookahead+Learn", VSIDSSelector(2), LookaheadConfig(enabled=true), LearnedClauseDB())
    ]

    for (name, selector, lookahead, learned) in strategies
        println("Testing: $name")
        for problem in problems
            result = solve_with_config(problem, selector, lookahead, learned)
            println("  $(problem.name): branches=$(result.total_branches), time=$(result.time)")
        end
    end
end
```

### 关键指标

1. **总分支数** (`problem.ws.total_branches`)
2. **总子问题数** (`problem.ws.total_subproblems`)
3. **平均分支因子** (`total_subproblems / total_branches`)
4. **求解时间** - 包含启发式开销
5. **缓存命中率** - Region contraction缓存效率

---

## 理论分析

### 分支因子的影响

假设平均分支因子为 `b`，搜索深度为 `d`：
- **时间复杂度**: O(b^d)
- **空间复杂度**: O(d)

减少分支因子的效果：
- `b=10` -> `b=5`: 在深度20时，节点数从 10^20 降到 10^14（百万倍加速）
- `b=5` -> `b=2`: 在深度20时，节点数从 10^14 降到 10^6（亿倍加速）

### 各策略的理论分支因子减少

| 策略 | 预期分支因子减少 | 开销 |
|------|------------------|------|
| VSIDS | 20-40% | 极低 |
| Look-ahead (k=16) | 30-50% | 中等 |
| Learned clauses | 40-60% (结构化问题) | 低 |
| 组合策略 | 60-80% | 中等 |

---

## 下一步优化方向

1. **自适应策略选择**: 根据问题特征自动选择最佳策略组合
2. **并行分支评估**: 多线程评估look-ahead
3. **高级冲突分析**: 实现1-UIP学习方案
4. **重启策略**: 周期性重启以跳出局部困难区域
5. **Phase saving**: 记住变量的历史赋值倾向
