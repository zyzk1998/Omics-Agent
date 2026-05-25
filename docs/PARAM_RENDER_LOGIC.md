# 参数推荐 UI 渲染全链路审查摘录

> 用途：供架构师 Code Review。本文档仅描述当前仓库**已实现**的数据流，不改动业务代码。  
> **类名勘误**：对话中出现的 `diagnosis-structured-rec-block` 在仓库中**不存在**；实际挂载「💡 参数推荐」子块的容器类名为 `**diagnosis-params-rec-block`**（见 `services/nginx/html/index.html`）。

---

## A. 数据源头：后端下发的 JSON 形态

### A1. SSE `diagnosis` 事件（上传文件后经编排器诊断）

编排器在诊断成功后 `_emit_sse(..., "diagnosis", payload)`，payload 顶层同时携带 `**diagnosis_report`（Markdown 字符串）** 与 `**recommendation`（结构化对象）**：

```2787:2803:gibh_agent/core/orchestrator.py
                        payload = {
                            "message": diagnosis_message,
                            "n_samples": n_samples,
                            "n_features": n_features,
                            "file_type": file_metadata.get('file_type'),
                            "status": "data_ready",
                            "recommendation": recommendation_data,  # 🔥 TASK 3: 添加参数推荐
                            "diagnosis_report": diagnosis_message,  # 🔥 TASK 2: 添加完整诊断报告
                        }
                        if integrity_status:
                            payload["integrity_status"] = integrity_status
                        # 🔥 DRY: Include suggestions from diagnosis so frontend can render chips
                        if agent_instance and getattr(agent_instance, "context", None):
                            sug = agent_instance.context.get("diagnosis_suggestions")
                            if sug:
                                payload["suggestions"] = sug
                        yield self._emit_sse(state_snapshot,"diagnosis", payload)
```

- `**data.diagnosis_report**`：完整诊断 Markdown（含 LLM 输出的「### 💡 参数推荐」表格正文等）。
- `**data.recommendation**`：来自 Agent `context["parameter_recommendation"]`，由 `_extract_parameter_recommendations` 解析 Markdown 表格得到，**规范形状**为 `{ "summary": string, "params": { "<param_name>": { "value", "reason", ... } } }`：

```2554:2564:gibh_agent/agents/base_agent.py
        Returns:
            参数推荐字典，格式：
            {
                "summary": "推荐摘要",
                "params": {
                    "param_name": {
                        "value": "推荐值",
                        "reason": "推荐理由"
                    }
                }
            }
```

**结论**：前端同时可能收到**两条通道**——叙述性 Markdown（`diagnosis_report`）与**已解析的** `recommendation.params`；二者在 UI 层由 `buildDiagnosisReportBodyHtml` 组合（见 B/C），**不是** `data.recommendation.params` 嵌在 `diagnosis_report` JSON 里（后者在典型 Path A 下是**字符串**）。

### A2. SSE `workflow` 事件（规划完成）

编排器在 `workflow_event_data` 顶层**可选**再挂一份 `recommendation` / `diagnosis_report`（来自 `self.agent.context`），与 `workflow_config` 并列（见 `gibh_agent/core/orchestrator.py` 中 `workflow_event_data["recommendation"]` / `["diagnosis_report"]` 附近逻辑）。

### A3. Planner 直接生成的 `recommendation`（无 Agent 诊断表解析时）

当 `recommended_params` 为**扁平**键值时，Planner 会生成与前端表格兼容的 `recommendation`：

```485:494:gibh_agent/core/planner.py
        # 前端“参数推荐”卡片：从 recommended_params 生成 recommendation（前端展示与后端执行一致）
        rec = workflow_plan.get("recommended_params")
        if isinstance(rec, dict) and rec and _is_flat_params(rec):
            workflow_config["recommendation"] = {
                "summary": "基于数据特征生成的参数推荐",
                "params": {
                    k: {"value": v, "reason": "基于数据特征推荐"}
                    for k, v in rec.items()
                },
            }
```

---

## B. 真实调用栈（至蓝底报告区 + `diagnosis-params-rec-block`）

右栏「数据诊断报告」外框的**蓝底**来自 `renderReportSection` 中 `report-content-box` 的**内联样式**（`background:#f0f7ff`），其**子内容**由 `buildDiagnosisReportBodyHtml` 生成；其中「💡 参数推荐」子块类名为 `**diagnosis-params-rec-block`**。

### B1. 主路径：在线 SSE

1. `fetch` 流式 Reader 循环解析 `event:` / `data:` → `JSON.parse` 得到 `data`（见 `index.html` 中约 `11742–11777` 行段）。
2. `**handleServerEvent(eventType, data)`**（约 `9780` 行起）。
3. `case 'diagnosis':` 当 `**data.diagnosis_report != null`** 时：
  `**renderReportSection('diagnosis', data.diagnosis_report, data)**`（约 `9878–9881` 行）。  
   （备用：若只有 `data.report_data.diagnosis` 则走 `report_data` 分支。）
4. `**renderReportSection**`（约 `10553` 行起）中 `section === 'diagnosis'`：
  - `**extractDiagnosisSectionInput(content, fullData)**` → 得到 `{ body, mergedFull }`（`mergedFull` 保留顶层 `recommendation` 等供子块使用）。  
  - `**buildDiagnosisReportBodyHtml(_dx.body, _dx.mergedFull)**` → 返回完整 HTML 字符串。
5. 拼接标题与蓝底容器后：
  `**target.innerHTML = bodyHtml**`（`target` 为 `.report-diagnosis-content`，约 `10594–10596` 行）。

**并行路径（聊天气泡里的工作流卡片，非右栏）**：`case 'workflow':` → `**renderWorkflowCard(data)`**（约 `9909–9915` 行）→ 内部若存在诊断/推荐则再次调用 `**buildDiagnosisReportBodyHtml(diagnosisReport, workflowData)`** 拼入 `formHtml`（约 `4572–4577` 行），最后写入卡片 DOM（同一文件 `_renderWorkflowCardImpl`）。

`**renderWorkflowForm(data)**`（约 `13653` 行起）：在工作流表单顶部同样用 `**buildDiagnosisReportBodyHtml(diagnosisReport, data)**`（约 `13730–13731` 行），标题文案为「数据诊断与参数推荐报告」，与右栏略有不同。

### B2. 降级路径：沉浸式 Demo / `handleServerEvent` 未挂载前

`**immersiveDemoFallbackHandleServerEvent**` → `case 'diagnosis':` → `**extractDiagnosisSectionInput**` → `**buildDiagnosisReportBodyHtml**` → 写入 `.report-diagnosis-content`（约 `1929–2016` 行段）。

---

## C. 嫌疑代码段：拼接「💡 参数推荐」的 HTML（原文照录）

仓库内**仅一处**字面量「💡 参数推荐」用于该子标题，位于 `**buildDiagnosisReportBodyHtml`**：

```1801:1818:services/nginx/html/index.html
            function buildDiagnosisReportBodyHtml(content, fullData) {
                fullData = fullData || {};
                var parts = [];
                var norm = normalizeReportMarkdownInput(content);
                if (/^\s*<div class="smart-content-wrapper\b/.test(norm)) {
                    parts.push('<div class="diagnosis-markdown-block markdown-body">' + norm + '</div>');
                } else if (norm) {
                    parts.push('<div class="diagnosis-markdown-block markdown-body">' + safeMarkedParse(norm) + '</div>');
                }
                var prm = fullData.recommendation && fullData.recommendation.params;
                if (prm && typeof prm === 'object' && !Array.isArray(prm) && Object.keys(prm).length) {
                    parts.push('<div class="diagnosis-params-rec-block mt-3"><h6 class="text-muted mb-2">💡 参数推荐</h6>' +
                        renderRecommendationParamsDict(prm) + '</div>');
                }
                if (!parts.length) {
                    parts.push('<p class="text-muted small">暂无诊断内容</p>');
                }
                return parts.join('');
            }
```

附：**结构化 `params` 字典 → 表格 HTML** 由 `**renderRecommendationParamsDict`** 生成（与「💡」标题拼接在同一分支内调用）。

---

## D. Markdown 引擎介入点：`marked.parse` / `safeMarkedParse`

- `**safeMarkedParse`** 定义（内部在满足条件时调用 `**marked.parse(s)`**）：

```1722:1741:services/nginx/html/index.html
            function safeMarkedParse(text) {
                if (text === null || text === undefined) return '';
                var s = normalizeReportMarkdownInput(text);
                // Sidecar / 挂载调试用前缀行不透出 UI
                s = s.replace(/^\s*(?:\[Local[^\]\n]*\]|local_workspace_mounted\s*=\s*\S+)\s*$/gim, '');
                // 剔除大模型生成的黄色警示图标
                s = s.replace(/⚠️/g, '');
                // #region agent log
                // 已禁用调试请求: fetch('http://localhost:7242/ingest/...',...).catch(()=>{});
                // #endregion
                if (/^\s*<div class="smart-content-wrapper\b/.test(s)) {
                    return s;
                }
                if (typeof marked !== 'undefined' && marked.parse) {
                    return marked.parse(s);
                }
                // 降级：简单的 HTML 转义和换行处理
                return s.replace(/&/g, '&').replace(/</g, '<').replace(/>/g, '>').replace(/\n/g, '<br>');
            }
```

- **喂给 Markdown 的那一段 payload**：在 `**buildDiagnosisReportBodyHtml`** 中，先 `**normalizeReportMarkdownInput(content)`** 得到 `norm`，再走 `**safeMarkedParse(norm)`** 进入 `**marked.parse**`（见上文 C 段 `else if (norm)` 分支）。  
- **不经过 Markdown 的部分**：同一函数内，若存在 `**fullData.recommendation.params`**，则由 `**renderRecommendationParamsDict`** 直接产出表格 HTML，**不调用** `marked.parse`。
- **右栏兜底**：若不存在 `buildDiagnosisReportBodyHtml`，`renderReportSection` 对 diagnosis 分支可退化为 `**safeMarkedParse(...)`**（约 `10571–10575` 行）。

---

## E. 根因与修复（实机 `[object Object],[object Object]`）

**根因（前端）**：`index.html` 中自定义 `**renderer.table`** 按旧版 **marked** 约定把 `header` / `body` 当作**已拼好的 HTML 字符串**；当前 CDN 未锁版本的 **marked** 在解析 Markdown 表格时，向 `renderer.table` 传入的是 **table token**（`header` / `rows` 为单元格 token 数组）。原实现使用 `String(tok.header)`，在数组上退化为 `**[object Object],[object Object]`**，经 `safeMarkedParse` → `marked.parse` 进入诊断区后，**LLM 生成的「### 💡 参数推荐」Markdown 表**整表损坏。

**修复**：

- 将 `**renderer.table`** 改为与 **marked 默认 `Renderer.table(token)`** 一致：用 `**this.tablecell` / `this.tablerow**` 展开 token，再包上原有 `table-responsive` + `md-table` 外壳；并保留 **两参数均为字符串** 时的旧版兼容分支。
- `**renderRecommendationParamsDict`** 中「推荐理由」列由 `**_escRecommendationCell(reason)`** 改为 `**_fmtRecommendationCell(reason)`**，避免 `reason` 为对象时再次退化为 `[object Object]`。

（行号以当前 `services/nginx/html/index.html` 为准，若后续增删行请以 `renderer.table` / `renderRecommendationParamsDict` 锚点搜索。）

---

## 审计备注（与 `[object Object]` 现象的关系）

- 诊断区「💡 参数推荐」**Markdown 表**走 `**buildDiagnosisReportBodyHtml` → `safeMarkedParse` → `marked.parse`**，与 结构化 `recommendation.params` 表（`renderRecommendationParamsDict`）为两条独立通道；上表问题主要由 §E 的 `renderer.table` 与 marked 版本契约不一致引起，下表需检查 `**recommendation.params[*].value` / `reason` 类型**是否非预期嵌套对象。

