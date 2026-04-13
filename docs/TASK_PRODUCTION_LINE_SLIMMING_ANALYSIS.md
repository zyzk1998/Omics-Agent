# Task 产线瘦身分析报告（Phase 2）

目标：保留 **`_generate_analysis_summary`**，把前置 LLM 压到约 2～3 次。范围：`orchestrator.py`、`planner.py`、`agentic.py`、`base_agent.py`、`data_diagnostician.py`。

| 标签 | 环节 | 说明 |
|------|------|------|
| A | `QueryRewriter.rewrite` | Task 段 Step 2，`achat` |
| B | `_classify_intent` | **编排器在 `generate_plan` 之前**已调一次；`generate_plan` 传 `domain_name` 可跳过内部第二次 |
| C | `_analyze_user_intent` | 编排器先调一次；`generate_plan` 传 `target_steps` 可跳过内部第二次 |
| D | 数据诊断 LLM | `DataDiagnostician.analyze` 无 LLM；`BaseAgent._perform_data_diagnosis` 内 `achat` 出 Markdown |

**关键点**：非 `target_domain` 时，**只改 `generate_plan` 不能省 B/C**，必须先改编排器前置意图块。

---

## A：旁路 QueryRewriter

- **触发**：`__init__` 有 LLM 则 `self.query_rewriter = QueryRewriter(...)`；`stream_process` 澄清合并分支不调用；否则 `await self.query_rewriter.rewrite(query, history)`。
- **风险**：`query_rewriter is None` 或跳过 rewrite → 仅用原始 `query`，关键词/Planner 匹配可能变差；**无异常**。
- **做法**：`GIBH_LITE_TASK_MODE` 为真时 `refined_query = query`，等价于不初始化 rewriter；澄清逻辑不变。

```python
elif _lite or not self.query_rewriter:
    refined_query = query
else:
    refined_query = await self.query_rewriter.rewrite(query, history)
```

---

## B：旁路 `_classify_intent`（astream）

- **现状**：`planner.generate_plan`：`domain_name` 非空则跳过 `_classify_intent`（`planner.py`）。但 `orchestrator` 在进 `generate_plan` 前已 `await planner._classify_intent(...)`。
- **风险**：错/空 `domain_name` → 不支持域或 `get_workflow` 失败。SemanticRouter 仅 `task`，**不能**直接当域名。
- **做法**：在编排器**调用 `_classify_intent` 之前**做 `domain_resolved`：`target_domain` 映射、`file_metadata.domain`、资产嗅探、或扩展 Router 输出；若 `is_supported(domain_resolved)` 则赋值 `domain_name` 并跳过 astream，否则回退原调用。`generate_plan` 继续传 `domain_name=`。

---

## C：降级 `_analyze_user_intent`

- **产出**：`{"target_steps": [step_id...], "skip_steps": [step_id...]}`，空 `target_steps` 表示全量 DAG；ID 须落在 `workflow.steps_dag`。
- **现状**：`planner._analyze_user_intent_core` 一次 `achat`；已有 `_fallback_intent_analysis`（关键词→步骤）；编排器在 steps 仍空时也会调回退。
- **风险**：纯规则难以覆盖「全分析但跳过某步」等 `skip_steps`；无回退则常落全量 DAG，与用户「只做 PCA」可能不一致。
- **做法**：LITE 下先 `_fallback_intent_analysis(q, list(workflow.steps_dag.keys()))`，得到 `{"target_steps": [...], "skip_steps": []}`；可选匹配不到再 `achat` 一次（混合）。

---

## D：规则拦截诊断 LLM

- **现状**：`DataDiagnostician.analyze` 纯统计；LLM 在 `BaseAgent._perform_data_diagnosis` 拼 prompt 后 **`achat`**。
- **风险**：规则直接返回短 Markdown 时，版式/参数表与现网不一致；规则过严误报。
- **做法**：在 **`diagnostician.analyze` 成功之后、`achat` 之前** 若 LITE 且 `_rule_based_diagnosis_markdown(stats, omics_type)` 非 `None` 则写缓存并 `return`，否则走原 LLM。

```python
if _lite and (md := _rule_based_diagnosis_markdown(stats, omics_type)):
    return md  # 例：missing_rate>0.5 / n_samples<3 等硬规则
```

---

## 特性开关 `GIBH_LITE_TASK_MODE`

- `1`/`true`：启用 A～D 旁路/降级；未设/`0`：与现网一致。实现方式对齐 `GIBH_ENABLE_SEMANTIC_ROUTER`（环境变量 + 默认关）。所有分支 **`if not _lite:` 保留原路径**。
- 测试：LITE 开/关各一条 Task 冒烟（可 mock LLM 调用次数）。

---

## 实施顺序

1. 开关 + **A**  
2. **B**（编排器域名决议）  
3. **C**（词典 + 可选混合 LLM）  
4. **D**（需产品确认无 LLM 报告形态）

---

*行号以仓库搜索为准。*
