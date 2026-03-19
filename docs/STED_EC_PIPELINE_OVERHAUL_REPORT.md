# STED-EC 技能全链路重构 — 三项任务排查与方案汇报

**状态**：待架构师批准后再动代码。

---

## 任务一：后端全量分析与内存防御 (Backend Full Analysis)

### 1.1 检索结论：是否存在数据拆分/降采样？

- **检索范围**：`gibh_agent/tools/sted_ec_tools.py`、`test_sted_ec_pipeline.py`、`gibh_agent/core/executor.py` 中 STED-EC 相关分支。
- **结论**：
  - **工具层无 downsample**：`sted_ec_data_validation`、`sted_ec_time_series_formatting`、`sted_ec_moscot_trajectory`、`sted_ec_plot_trajectory` 均直接对传入的 `h5ad_path` 做 `sc.read_h5ad(path)`，没有任何 `n_cells` 截断、`adata[:n]` 或随机采样。
  - **测试脚本的“切片”**：`test_sted_ec_pipeline.py` 注释中的「自动切片大文件」仅指**选用 mini_sted_ec.h5ad（约 3000 细胞）做快速测试**，工具本身仍读取「传入的路径」对应的文件；生产/前端传入全量 `sted-ec.h5ad` 时，工具会全量读入并计算。
- **代码切入点**：无需在工具内“移除”downsample，只需**保持全量读取**，并在以下位置增加**内存防护与中间释放**。

### 1.2 全量计算下的 OOM 风险与防护方案

| 风险点 | 位置 | 防护手段 |
|--------|------|----------|
| 单次读入大 h5ad | 各工具内 `sc.read_h5ad` | 保持 `adata.X` 为 **scipy.sparse**（Scanpy 默认即 sparse）；写入后若不再需要可 `del adata; gc.collect()`。 |
| 时间对循环内大对象 | `sted_ec_moscot_trajectory` 中 `merged = anndata.concat([adata1, adata2])`、`TemporalProblem(merged).solve(...)`、`cps_adata` | 每处理完一对时间点 (t1, t2)：写入 tmap 后立即 `del merged, sub, tp, cps_adata, transport`（若存在），再 `import gc; gc.collect()`，再进入下一对。 |
| PCA/邻域/UMAP | `sted_ec_plot_trajectory` 中 `sc.tl.pca`、`sc.pp.neighbors`、`sc.tl.umap` | 全量时 `obsm['X_pca']` 会增大内存；保持 `adata.X` 为 sparse；计算完成后若后续不再需要可对 adata 做轻量引用或只保留作图所需列。 |
| 多时间对累积 | moscot 循环多次 `adata1 = adata[mask_t1].copy()` | 避免在循环外保留多份 adata 副本；循环内用 `adata[mask_t1]` 得到子集，用完后不保留引用。 |

**建议的具体改动（仅思路，不直接改代码）**：

1. **sted_ec_moscot_trajectory**  
   - 在 `for i in range(len(times) - 1):` 循环末尾（写入当前 tmap 的 `cps_adata.write(out_path)` 之后）：  
     - `del merged, sub, tp, cps_adata`（及 `transport` 若仍存在）。  
     - 调用 `gc.collect()`。  
   - 不改变循环逻辑与全量 adata 读入，仅增加显式释放。

2. **sted_ec_time_series_formatting / sted_ec_data_validation**  
   - 在 `adata.write(...)` 或 return 前，若该函数内不再使用 adata，可 `del adata; gc.collect()`（可选，视内存监控情况）。

3. **sted_ec_plot_trajectory**  
   - 若在计算完 UMAP 并完成所有作图后，不再需要完整 adata，可在 return 前 `del adata; gc.collect()`（可选）。

4. **文档/注释**  
   - 在 `sted_ec_tools.py` 顶部或各函数 docstring 中注明：本流程为**全量分析**，无降采样；大文件场景下依赖 sparse 矩阵与循环内显式 `gc.collect()` 降低 OOM 风险。

---

## 任务二：中间态数据流装配 (Intermediate Data Streaming)

### 2.1 执行链路与 step_result 去向

- **执行路径**：`orchestrator` 直接执行 → `WorkflowExecutor.execute_workflow` → 对每个 step 调用 `execute_step` → 从 `tool_registry` 调对应 Tool；Tool 返回的 dict 被包装进 `step_detail`（含 `step_result`），再经 `step_result` 事件通过 SSE 推给前端。
- **现有 step_result 结构**：`executor.py` 中 `step_detail` 含 `step_id`、`name`、`status`、`summary`、`step_result: { step_name, status, logs, data }`；若 Tool 返回 `report_data`（如 `images`、`download_links`），会在 executor 中展平到 `step_result.data` 或 `report_data`，并由 orchestrator 写入 `state_snapshot` 与 SSE。
- **结论**：只要 Tool 返回结构中含有**轻量级**的 `summary`、`report_data.images`、`report_data.download_links` 等，现有链路就会把它们推到前端；**不能**在返回值里塞入几百 MB 的 JSON（如完整 adata 或巨大 DataFrame）。

### 2.2 各步骤当前返回值与建议补充（仅摘要/链接，无大 JSON）

| 步骤 | 当前返回 | 建议补充（均放入 return dict，供 executor 原样转发） |
|------|----------|------------------------------------------------------|
| **sted_ec_data_validation** | status, h5ad_path, time_key, cell_type_key | **n_obs**, **n_vars**；**summary**：一段短文本，如「数据校验通过，共 {n_obs} 细胞、{n_vars} 基因；时间列: {time_key}，细胞类型列: {cell_type_key}」。 |
| **sted_ec_time_series_formatting** | status, h5ad_path | **n_obs**（过滤后）、**summary**：如「时间序列标准化完成，保留 {n_obs} 细胞，已写出 {formatted_path}」。不返回大表。 |
| **sted_ec_moscot_trajectory** | status, message, h5ad_path, tmaps_dir, time_key | **summary**：如「已完成 {N} 个时间对 OT 求解，transport 矩阵已写入 tmaps」；可选 **report_data.download_links**：如 `[{ "title": "tmaps 目录", "path": tmaps_dir }]`（若前端支持目录链接则加，否则仅 summary）。不返回 transport 矩阵或大对象。 |
| **sted_ec_plot_trajectory** | 已有 report_data.images、download_links、summary | 保持现状；可把 summary 改为稍详细的 Markdown（如「已生成时空 UMAP、细胞类型 UMAP、细胞类型演化图；ZIP 已打包」）。 |

- **约定**：所有「表格摘要」仅限**少量统计量**（如 n_obs、n_vars、时间点数、图表数量）；禁止在 return 中附带 `adata.obs`、大 DataFrame、或完整矩阵。若需“中间结果表”，只提供**后端生成的 Markdown 摘要字符串**或**下载链接**。

### 2.3 关键代码切入点

- **sted_ec_tools.py**  
  - `sted_ec_data_validation`：在 return 前增加 `n_obs=adata.n_obs`、`n_vars=adata.n_vars`，以及 `summary=...` 字符串。  
  - `sted_ec_time_series_formatting`：在 write 后增加 `n_obs=adata.n_obs`、`summary=...`。  
  - `sted_ec_moscot_trajectory`：在循环结束后增加 `summary=...`（含时间对数量）、可选 `report_data.download_links`。  
  - `sted_ec_plot_trajectory`：已有 report_data；可加强 `summary` 文本。  
- **executor.py**  
  - 当前已把 Tool 返回的 `result` 写入 `step_result` 和 `step_detail`，并会处理 `report_data.images` 等；只需保证 Tool 侧 return 的 key 与现有逻辑一致（如 `report_data.images`、`report_data.download_links`、`summary`），无需改 executor 的装配逻辑（除非需要把 `summary` 显式挂到 step 的某字段，当前可用 `message`/`logs` 或 step_result 内已有字段承载）。

---

## 任务三：前端预览卡片“强交互”重写 (Frontend Preview Hack)

### 3.1 当前 template_mode 与 disabled 逻辑

- **位置**：`services/nginx/html/index.html` 中 `_renderWorkflowCardImpl` 内，约 6676–6689 行：`isTemplateMode` 由 `workflowData.template_mode` / `workflow_config?.template_mode` / `workflow_data?.template_mode` 或占位符 `hasPlaceholder` 决定。
- **禁用效果**（约 6809、6814、6871、6879、6887 行）：  
  - `isTemplateMode === true` 时：  
    - 步骤 **checkbox** 加上 `disabled`。  
    - 步骤 **accordion 按钮** 加上 `style="pointer-events: none; opacity: 0.6;"`。  
    - 参数 **input**（checkbox/number/text）加上 `disabled` 或 `readonly style="background-color: #e9ecef;"`。  
  - 路径类参数在模板模式下本身为只读占位（如「上传数据以激活」），这是合理的，不改。
- **结论**：预览模式下，**非路径**的 checkbox 与 input 被统一禁用，导致用户无法在「载入官方演示数据」前勾选步骤或填写生信参数。

### 3.2 方案：仅对 STED-EC 放开预览交互，不破坏其他技能

- **思路**：在渲染工作流卡片时，**仅当当前工作流为 STED-EC 时**，在预览模式（template_mode === true）下仍**不**对 checkbox 与 input 施加 disabled/readonly；其余技能保持现状。
- **识别 STED-EC**：  
  - 方式 A：`workflowData.workflow_name` 或 `workflowData.workflow_data?.workflow_name` 含 `"STED"` / `"sted"` / `"轨迹"`。  
  - 方式 B：任意一步的 `tool_id` / `step.step_id` 含 `sted_ec`（与现有「载入官方演示数据」判断一致，见 8451 行）。  
  - **建议**：与「载入官方演示数据」一致，用 **方式 B**（存在 `sted_ec` 的 step_id/tool_id）判定为 STED-EC，避免依赖 workflow_name 命名。
- **具体逻辑**：  
  - 在生成 `formHtml` 的循环前，计算一次：  
    - `const isStedEcWorkflow = workflowSteps.some(s => (s.tool_id || s.step_id || '').indexOf('sted_ec') !== -1);`  
  - 定义：`const allowPreviewInteraction = isStedEcWorkflow;`  
  - 将所有「在 template 下禁用」的条件从 `isTemplateMode` 改为 `isTemplateMode && !allowPreviewInteraction`。  
  - 即：**仅当「模板模式且非 STED-EC」时才 disabled/readonly**；STED-EC 在预览模式下可勾选步骤、可填参数。
- **代码切入点**：  
  - 约 6785 行附近、`workflowSteps.forEach` 之前：增加 `isStedEcWorkflow` 与 `allowPreviewInteraction`。  
  - 6809、6814、6871、6879、6887 行：把 `isTemplateMode` 替换为 `(isTemplateMode && !allowPreviewInteraction)`。

### 3.3 参数注入：生信通用默认值（仅 STED-EC）

- **来源**：根据 `sted_ec_tools.py` 中参数与常见生信取值，拟定一组**前端只读配置**（不请求后端），在渲染 STED-EC 卡片时用于填充未在 `step.params` 中出现的字段。  
- **建议默认值**（可与产品再对齐）：  

  - **sted_ec_data_validation**：time_key=`"day"`，cell_type_key=`"cell_type"`（已有）。  
  - **sted_ec_time_series_formatting**：time_key=`"day"`，cell_type_key=`""`（空则后端嗅探）。  
  - **sted_ec_moscot_trajectory**：n_pcs=`40`，epsilon=`1e-3`，tau_a=`0.99`，tau_b=`0.999`，scale_cost=`"mean"`，max_iterations=`20`，time_key=`"day"`，cell_type_key=`""`。  
  - **sted_ec_plot_trajectory**：plot_type=`"umap"`，time_key=`"day"`，cell_type_key=`""`。

- **注入方式**：  
  - 在 `_renderWorkflowCardImpl` 内、遍历 `stepParams` 生成 input 时，若当前 workflow 为 STED-EC 且 `stepParams[paramName]` 为空/占位，则从一份静态对象（如 `window.STED_EC_DEFAULT_PARAMS` 或内联对象）按 `stepId + paramName` 取默认值作为 `defaultValue` 显示在 input 中。  
  - 这样预览时用户能看到并可修改这些默认值；**载入官方演示数据** 时仍通过 `getWorkflowDataFromForm(formId)` 收集当前表单值（含用户改过的），无需改后端。

### 3.4 状态流转：载入官方演示数据时参数如何带到后端

- **当前逻辑**（约 8475–8491 行）：  
  - 点击「载入官方演示数据」→ 取 `formId`，对 time_key/cell_type_key 做前端赋值（day/空）→ `getWorkflowDataFromForm(formId)` 得到 `data.workflow_data`（含 steps 与每步 params）→ 将路径类参数统一替换为 demo 路径 → `sendMessage(..., workflowPayload, ..., demoFiles)`。  
- **结论**：前端已把**当前卡片上的所有参数**（包括用户修改的）通过 `workflow_data.steps[].params` 发给后端；后端 executor 在 `execute_step` 时使用的是 step 的 `params`，即**前端传来的参数已优先**。  
- **建议**：  
  - 保持该流程不变。  
  - 仅在「强交互」改造后，确保 STED-EC 预览卡片上输入的默认值（见 3.3）在初次渲染时写入 input，这样 `getWorkflowDataFromForm` 会自然带上这些值；若用户修改，则带修改后的值。  
  - 后端无需“优先使用前端参数”的额外逻辑，现有设计已满足。

---

## 小结与实施顺序建议

| 任务 | 关键点 | 建议实施顺序 |
|------|--------|--------------|
| **任务一** | 工具层已是全量；在 moscot 循环内及必要时在其它工具末尾增加 `del` + `gc.collect()`，并补充注释/文档。 | 先做，风险低。 |
| **任务二** | 各 Tool return 中增加 `n_obs`/`n_vars`/`summary` 及已有结构的 `report_data`，不返回大 JSON。 | 与任务一可并行或紧随其后。 |
| **任务三** | 前端：`allowPreviewInteraction = isStedEcWorkflow`，仅 STED-EC 在 template 下可编辑；注入 STED_EC 默认参数；载入演示数据流程不变。 | 最后做，便于联调。 |

请架构师审查上述三项任务的思路与代码切入点，确认后再进行具体代码修改。
