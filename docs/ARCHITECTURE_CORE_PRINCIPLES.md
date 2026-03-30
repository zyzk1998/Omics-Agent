# GIBH-AGENT-V2 核心架构原则（白皮书）

本文档沉淀多轮重构后**不可破坏**的架构约束。扩展功能、修复 Bug 时须与之对齐；细节实现以代码为准。

---

## 1. 多模态资产总线（OmicsAssetManager）

**位置**：`gibh_agent/core/asset_manager.py`（`DataAsset`、`OmicsAssetManager.classify_assets`、`to_legacy_resolved_dict`）

**目标**：将用户上传的零散路径规范为**结构化资产列表**，避免在编排器、执行器、Planner 中硬编码「若文件名含 xxx 则拼接 yyy」。

**机制概要**：

- **后缀与类型嗅探**：通过 `_suffix_lower`、动态 `SUFFIX_TO_ASSET_TYPE` 等将单文件映射到 `protein_structure`、`metabolomics_csv`、`h5ad_single`（仅 `.h5ad`）、`10x_h5_matrix`（10x 原生 `.h5` 矩阵）等。
- **目录与内容特征**：例如 `_is_10x_mtx_directory`（含 `matrix.mtx`）、`_is_spatial_image_directory`（目录名或 Visium 典型文件）、Radiomics 的影像与 mask 关键词配对。
- **捆绑优先于单文件**：在同一批路径上按固定顺序尝试聚合——如 `10x_bundle`（matrix + barcodes/features）、`spatial_bundle`（矩阵端 + spatial 目录，`metadata` 含 `spatial_dir` / `visium_root`）、`radiomics_pair`（图 + mask），剩余再落单文件或 `generic_unknown`。
- **与旧管线兼容**：`to_legacy_resolved_dict` 将资产写入与 `resolve_omics_paths` 一致的桶（`tables`、`h5ad`、`10x_mtx`、`images`、`masks`、`unknown` 等），使体检与执行器推断不断裂。

**禁止**：在新工具或编排逻辑里手写「拼出 Visium 根路径」的长期方案；应扩展 `OmicsAssetManager` 或注册新 `asset_type` 并接入 legacy 映射。

---

## 2. 反射参数适配器（Executor Reflection）

**位置**：`gibh_agent/core/executor.py`（`_tool_accepted_param_names`、`_resolve_callable_for_signature`、`_inject_paths_from_asset_reflective`、`_filter_and_coerce_params_to_signature`、`_process_data_flow`）

**目标**：工具路径形参名由**真实函数签名**决定，执行器按 `DataAsset` 类型与签名**反射注入**，避免 `if tool_id == "foo"` 分支爆炸。

**机制概要**：

- 通过 `inspect.signature(inspect.unwrap(tool))` 得到允许的参数名；`VAR_KEYWORD` 时保留额外键。
- `_INJECTABLE_PATH_PARAM_KEYS` 等集合描述可注入槽位；按 `spatial_bundle` / `10x_bundle` / `radiomics_pair` 等资产类型写入 `adata_path`、`matrix_dir`、`data_path`、`image_path`、`mask_path` 等（仅当槽位为空）。
- `_coerce_value` 依据注解做类型转换（如逗号分隔字符串转为 `List[float]`），减少前端与 LLM 输出与底层 API 的摩擦。
- 占位符 `<step_id>` 解析与「提前注入」的先后顺序有严格语义（存在 `<...>` 时不应污染 `adata_path` 等），修改时需回归空间工作流链式测试。

**禁止**：在执行器中为单个工具新增「专用 if 工具名」的路径注入；应扩展反射规则、资产类型或工具签名。

---

## 3. TaaS 微服务架构（Worker 容器）

**原则**（与 `.cursor/rules/gibh-architecture-constitution.mdc` 一致）：**重型生信算子不得在主智能体进程内跑满 CPU/GPU**。

**形态**：

- 独立 **Docker / Conda** 服务（Worker），暴露 HTTP JSON API；主进程只组装参数、轮询或等待结果。
- 用户数据通过**共享挂载**（如 `/app/uploads`、`RESULTS_DIR`）以**路径字符串**传递，不在进程间拷贝大二进制。
- 典型场景：RNAFold、PyMOL、BepiPred 类工具链、长时程结构预测等。

**禁止**：在 `server.py` 或 Agent 主线程中直接调用长时间阻塞的第三方 CLI 而不做子进程/服务隔离（除非已有明确的轻量封装且符合宪法条款）。

---

## 4. 技能快车道（Skill Fast-Lane）

**位置**：`gibh_agent/core/orchestrator.py`（`[Skill_Route: tool_id]` 正则拦截）、`gibh_agent/agents/skill_agent.py`（参数抽取与工具调用）、`gibh_agent/db/seed_skills.py`（模板首行约定）

**目标**：用户或技能模板在查询中嵌入**可解析暗号**，在全局意图分类与 SOP Planner 之前**短路**到「单工具技能执行」，降低延迟与幻觉规划。

**机制概要**：

- **`[Skill_Route: registered_tool_name]`**：与 `registry.register(name=...)` **逐字一致**；Orchestrator 匹配后实例化 `SkillAgent`，流式 SSE 输出，不再走笨重全局规划。
- **`[Omics_Route: ...]`**：在技能广场列表等场景与 `Skill_Route` 一并用于**排序/识别快车道模板**（见 `gibh_agent/api/routers/skills.py`、`server.py`）；可与业务上「组学域」模板约定结合使用。
- 另有 **`target_domain`** 等 API 层快车道，跳过部分意图分类（见 Orchestrator 中 `target_domain` 分支），与 Skill 暗号互补而非重复实现业务逻辑。

**禁止**：为单个技能在 Orchestrator 顶层再写一套平行路由表；新技能应注册工具 + 模板暗号或复用 `target_domain` 契约。

---

## 5. 前端 UI 防爆与流式渲染

**主要资产**：`services/nginx/html/index.html`、`services/nginx/html/css/style.css`、`services/nginx/html/css/main.css`

**CSS 隔离**：

- 用户消息、气泡、代码块等样式**限定在** `.message-row.user`、`.message-row.assistant` 等选择器下，避免全局标签污染导致布局崩坏或主题串味。
- 大文件 `index.html` 修改须遵守 UTF-8 与**锚点式最小 diff**（见 `.cursor/rules/large-html-utf8-and-css.mdc` 与 `gibh-architecture-constitution` 前端条款）。

**`[AgenticLog]` 流式树**：

- SSE 下行中形如 `[AgenticLog] ` 前缀的 JSON 载荷由 `handleAgenticLog` 解析，挂载到 `.process-steps` 等容器。
- **微光呼吸**：步骤行使用 `agent-step__glyph--pulse` 等类名做进行态强调。
- **可折叠子步骤**：`bindAgentSubstepToggle`、chevron 旋转，避免思考与工具细节刷屏。
- 用户可见正文侧可对连续 `[AgenticLog]` 行做折叠/省略摘要，平衡「可观测」与「可读」。

**禁止**：为省事先全局覆盖 `index.html`；禁止破坏现有 Flex 布局契约；新增流式 UI 应复用上述挂载点与 BEM 风格类名。

---

## 6. 与项目宪法的关系

- 路径类工具参数命名、TaaS 黑盒、`@safe_tool_execution`、前端微创等条款见 **`.cursor/rules/gibh-architecture-constitution.mdc`**。
- **开发前**：先读本白皮书与宪法，再改代码。

---

*文档版本：与仓库同步维护；重大架构变更须更新本文并评审。*
