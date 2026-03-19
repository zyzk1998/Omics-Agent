# 方向 A 深度重构设计：工作流内 LLM 流式穿透

## 一、已定位的 LLM 消耗点（第一步结论）

| 位置 | 文件 | 函数 | 说明 |
|------|------|------|------|
| **专家解读报告** | `gibh_agent/agents/base_agent.py` | `_generate_analysis_summary` | 约 865–2010 行；约 1902 行 `completion = await llm_client_to_use.achat(messages, temperature=0.3, max_tokens=2500)`，生成最终分析报告正文 |
| **调用方** | `gibh_agent/core/orchestrator.py` | `stream_process`（直接执行路径） | 约 657–718 行：在 `WorkflowExecutor.execute_workflow()` 返回后，调用 `summary = await target_agent._generate_analysis_summary(...)`，再将 `summary` 填入 diagnosis 事件 |

**结论**：专家解读报告由 **orchestrator 直接调用 agent._generate_analysis_summary** 生成，不经过 WorkflowExecutor。因此「Generator 穿透」只需：**Agent 层 yield thought → Orchestrator 层消费并 _emit_sse(thought)**。WorkflowExecutor 无需改动。

---

## 二、第二步核心：achat → astream + 边流边存 + yield thought

### 2.1 设计要点

- 在 **base_agent** 新增异步生成器 **`_stream_analysis_summary`**，与 `_generate_analysis_summary` 同参、同前置逻辑（构建 `messages`、`_clean_report_content`、retry 所需上下文）。
- 将「单次 achat + 解析」替换为：
  - 使用 **`llm_client.astream(messages, ...)`** 作为流式输入；
  - 使用现有 **`stream_from_llm_chunks(chunk_iter, model_name)`** 做双轨解析（`reasoning_content` / `<think>` → thought，普通 content → message）；
  - **边流边存**：对 `("message", data)` 将 `data["content"]` 追加到本地 `content_buffer`；
  - **实时向上抛出思考**：对 `("thought", data)` 直接 **yield ("thought", data)**，供上层转发 SSE；
  - 流结束后用 **`_clean_report_content(content_buffer)`** 得到与重构前等价的最终报告文本；若长度过短则走现有 retry（achat）逻辑，最终 **yield ("summary", final_report_str)**。

### 2.2 防回退

- **最终报告内容**：`content_buffer` 累积全部 `message` 片段后，再经同一套 `_clean_report_content`（含 <<<SUGGESTIONS>>> 剥离与 context 写入），与当前 achat 路径一致，**不因流式丢失字符**。
- **status / process_log / step_result**：不在此函数内改动；orchestrator 侧仅在「生成报告」这一段改为消费 `_stream_analysis_summary` 并多转发 thought，其余 status/step_result/diagnosis 逻辑不变。

### 2.3 核心代码片段（base_agent.py）

以下为 **第二步** 的核心逻辑：在 base_agent 中新增 `_stream_analysis_summary`，在「已有 messages、_clean_report_content、llm_client_to_use」的前提下，用 astream + stream_from_llm_chunks 替代 achat，并 yield thought / summary。

```python
# 在 base_agent.py 中新增（建议放在 _generate_analysis_summary 之后）

async def _stream_analysis_summary(
    self,
    steps_results: List[Dict[str, Any]],
    omics_type: str = "Metabolomics",
    workflow_name: str = "Analysis Pipeline",
    summary_context: Optional[Dict[str, Any]] = None,
    output_dir: Optional[str] = None,
) -> AsyncIterator[Tuple[str, Any]]:
    """
    流式生成分析摘要：边流边向调用方 yield ("thought", data)，
    流结束后 yield ("summary", final_report_str)。
    与 _generate_analysis_summary 同参、同前置逻辑，仅将 achat 改为 astream + 双轨解析。
    """
    from ..core.stream_utils import stream_from_llm_chunks

    # ----- 与 _generate_analysis_summary 完全一致的前置逻辑 -----
    # 1) 读取 execution results、构建 results_summary、key_findings、messages、_clean_report_content、
    #    llm_client_to_use、retry 所需上下文（key_findings_json、successful_steps、steps_results 等）
    #    此处省略：与现有 _generate_analysis_summary 中 try 块内、在 "开始LLM调用" 之前的部分相同。
    #    建议抽成 _prepare_analysis_summary_request() 供两处共用，返回 (messages, llm_client_to_use, _clean_report_content, retry_context)。
    messages, llm_client_to_use, _clean_report_content, retry_context = await self._prepare_analysis_summary_request(
        steps_results, omics_type, workflow_name, summary_context, output_dir
    )
    if not messages or not llm_client_to_use:
        yield ("summary", None)
        return

    content_buffer = ""

    # ----- 将 achat 替换为 astream + stream_from_llm_chunks -----
    chunk_iter = llm_client_to_use.astream(messages, temperature=0.3, max_tokens=2500)
    model_name = getattr(llm_client_to_use, "model", None)

    async for event_type, data in stream_from_llm_chunks(chunk_iter, model_name=model_name):
        if event_type == "thought":
            # 立即向上抛出，orchestrator 将 _emit_sse(state_snapshot, "thought", data)
            yield ("thought", data)
        elif event_type == "message":
            content_buffer += (data.get("content") or "")
        elif event_type == "suggestions" and isinstance(data, list):
            if not hasattr(self, "context"):
                self.context = {}
            self.context["report_suggestions"] = data

    # ----- 流结束：与 achat 路径等价的最终报告 -----
    final_summary = _clean_report_content(content_buffer) if content_buffer else ""

    # 过短时走现有 retry（achat）逻辑，保证与重构前行为一致
    if not final_summary or len((final_summary or "").strip()) < 100:
        retry_result = await self._retry_analysis_summary_achat(retry_context, llm_client_to_use, _clean_report_content)
        final_summary = retry_result or final_summary

    yield ("summary", final_summary)
```

- **`_prepare_analysis_summary_request`**：从当前 `_generate_analysis_summary` 中抽出「从 steps_results 到构建 messages、_clean_report_content、retry 上下文」的全部逻辑，返回上述四元组，供 `_generate_analysis_summary`（achat 路径）与 `_stream_analysis_summary`（astream 路径）共用，避免重复与行为分叉。
- **`_retry_analysis_summary_achat`**：将当前 `_generate_analysis_summary` 内「重试 achat + extract_think_and_content + 返回 cleaned 字符串」封装为独立方法，供流式路径在「content_buffer 过短」时调用，保证最终报告与重构前一致。

---

## 三、第三步核心：Orchestrator 侧 Generator 穿透（Event Bubbling）

### 3.1 设计要点

- **不改变**现有 status / step_result / diagnosis 的触发顺序与 payload 结构。
- 仅在「生成专家解读报告」这一段：若 agent 支持 `_stream_analysis_summary`，则 **async for 消费其事件**；对 `("thought", data)` 立即 **yield self._emit_sse(state_snapshot, "thought", data)**；对 `("summary", data)` 赋给局部变量 `summary`，后续原有 diagnosis/result/done 逻辑**完全沿用**该 `summary`。

### 3.2 核心代码片段（orchestrator.py）

在 `stream_process` 中，将当前「调用 _generate_analysis_summary 并得到 summary」的整块替换为下面逻辑；**前后文的 status、step_result、diagnosis、result、done 保持不变**。

```python
# 在 orchestrator.py 中，替换当前「summary = await target_agent._generate_analysis_summary(...)」所在块

# 🔥 流式穿透：优先使用 _stream_analysis_summary，将 thought 实时 _emit_sse 给前端
stream_method = getattr(target_agent, "_stream_analysis_summary", None)
if stream_method and callable(stream_method):
    summary = None
    try:
        async for evt_type, evt_data in stream_method(
            steps_results=steps_results,
            omics_type=domain_name,
            workflow_name=workflow_config.get("workflow_name", "工作流"),
            summary_context=summary_context,
            output_dir=output_dir,
        ):
            if evt_type == "thought":
                yield self._emit_sse(state_snapshot, "thought", evt_data)
                await asyncio.sleep(0.01)
            elif evt_type == "summary":
                summary = evt_data
                break
    except Exception as stream_err:
        logger.warning("⚠️ [Orchestrator] _stream_analysis_summary 异常，回退 achat: %s", stream_err)
        summary = await target_agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type=domain_name,
            workflow_name=workflow_config.get("workflow_name", "工作流"),
            summary_context=summary_context,
            output_dir=output_dir,
        )
    if summary is None:
        summary = await target_agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type=domain_name,
            workflow_name=workflow_config.get("workflow_name", "工作流"),
            summary_context=summary_context,
            output_dir=output_dir,
        )
else:
    summary = await target_agent._generate_analysis_summary(
        steps_results=steps_results,
        omics_type=domain_name,
        workflow_name=workflow_config.get("workflow_name", "工作流"),
        summary_context=summary_context,
        output_dir=output_dir,
    )

# 此后保持不变：elapsed_time、last_llm_error、后备 summary、evaluation、step_result、diagnosis、result、done
```

- **防回退**：若 `_stream_analysis_summary` 抛异常或未 yield `"summary"`，则回退到一次 `await _generate_analysis_summary(...)`，保证原有业务逻辑与报告内容一致。
- **不阻断状态机**：thought 是**额外**的 yield，不替代、不跳过任何既有 status / step_result / diagnosis；最终仍用同一个 `summary` 写入 diagnosis 与 state_snapshot。

---

## 四、与现有状态机的兼容性（防回退检查）

| 项目 | 说明 |
|------|------|
| **status / process_log** | 仍由 orchestrator 在「正在生成专家解读报告...」等处 yield；流式报告阶段仅**新增**对 thought 的 yield，不删不改既有 status。 |
| **step_result** | 仍在 execute_workflow 返回后、报告生成前/后按原逻辑 yield，不受 _stream_analysis_summary 影响。 |
| **diagnosis / result / done** | 仍使用同一变量 `summary` 构造 diagnosis_response 和 result；summary 来源从「单次 await _generate_analysis_summary」改为「async for _stream_analysis_summary 得到 summary 或回退到 _generate_analysis_summary」。 |
| **报告内容一致性** | 最终报告 = _clean_report_content(content_buffer)，与当前 achat 路径下对 completion.choices[0].message.content 的清洗逻辑一致；retry 条件与实现与现有一致，不因流式而丢失字符或改变语义。 |

---

## 五、建议实施顺序

1. **base_agent**：实现 `_prepare_analysis_summary_request`、`_retry_analysis_summary_achat`，并实现 `_stream_analysis_summary`（如上第二节片段）；保留 `_generate_analysis_summary` 改为内部调用 `_prepare_analysis_summary_request` + achat + 现有解析与 return。
2. **orchestrator**：在直接执行路径中，将「生成专家解读报告」块替换为第三节的「流式优先 + 回退 achat」逻辑。
3. **回归**：跑现有分析任务（含 STED-EC 与专家解读报告），确认 thought 事件在控制台/前端出现，且 diagnosis 内容与重构前一致。

---

请您重点审查：  
- 第二节中 **content_buffer 只累积 message、thought 仅向上 yield、最终 summary 来自 _clean_report_content(content_buffer) 或 retry** 是否满足「边流边存、报告与重构前一致」；  
- 第三节中 **仅多 yield thought、不改变 summary 的用途与后续 diagnosis/result/done** 是否满足「不阻断原有工作流状态机」。  
确认无误后，再进行全量代码替换与提交。

---

## 六、Generator 穿透链一览（供审查 yield 是否阻断状态机）

```
[LLM 流] astream(messages)
    → stream_from_llm_chunks(chunk_iter)   [stream_utils，已有]
        → yield ("thought", {"content": reasoning_chunk})   # reasoning_content 或 <think> 内
        → yield ("message", {"content": text_chunk})        # 正文
        → yield ("suggestions", list)                       # 可选

[Agent] _stream_analysis_summary
    async for event_type, data in stream_from_llm_chunks(...):
        if event_type == "thought":
            yield ("thought", data)          # 直接向上冒泡，不存
        elif event_type == "message":
            content_buffer += data.get("content", "")   # 边流边存
    # 流结束
    final_summary = _clean_report_content(content_buffer)
    yield ("summary", final_summary)         # 仅 yield 一次，作为该步骤结果

[Orchestrator] stream_process（直接执行路径）
    async for evt_type, evt_data in target_agent._stream_analysis_summary(...):
        if evt_type == "thought":
            yield self._emit_sse(state_snapshot, "thought", evt_data)   # 实时推前端，不阻断
        elif evt_type == "summary":
            summary = evt_data
            break
    # 后续：summary 照常用于 diagnosis / result / done，与现有一致
```

- **thought**：从 stream_utils → Agent yield → Orchestrator 消费并 _emit_sse → 前端；不参与 content_buffer，不改变 step_result/diagnosis 内容。
- **message**：仅在 Agent 内累积到 content_buffer，最后一次性变为 summary；Orchestrator 不收到 "message" 事件，只收到一次 "summary"。
- **原有 status/step_result/diagnosis/result/done**：在 orchestrator 中照常按原顺序 yield，summary 来源从「单次 await」改为「async for 收尾的 summary」或回退的 await，**不增加、不跳过任何既有事件**。
