# 技能详情 Demo 页面规范

> **版本**：Phase 1.1（抽屉重构 + demo_visualization）  
> **适用范围**：技能广场卡片「详细描述」→ 右侧宽幅抽屉 `#skill-detail-drawer`  
> **数据源码**：`skill_detailed_specs.py`（Phase1/2 手工）+ `skill_detailed_specs_bulk.py`（批量）+ `skill_detail_demo_visualizations.py`

---

## 1. 目标与入口

| 项 | 说明 |
|----|------|
| 入口 | 技能卡片底栏：收藏 → **详细描述** → 复制 Prompt → 使用 |
| 交互 | 点击「详细描述」打开 **900px 右侧抽屉**（小屏全宽）；含 **⭐ 收藏** 与 **🚀 立即使用** |
| API | `GET /api/skills` 附带 `detailed_spec`；`GET /api/skills/detailed-specs` 返回全量索引 |

---

## 2. `detailed_spec` JSON 结构

### 2.1 字段定义

| 字段 | 类型 | 必填 | 说明 |
|------|------|------|------|
| `tool_id` | string | 是 | 与 `[Skill_Route: <tool_id>]` 一致 |
| `description_long` | string | 是 | 业务功能长描述 |
| `usage_hint` | string | 否 | 顶部蓝色提示条 |
| `inputs` | array | 是 | 输入项列表 |
| `outputs` | array | 是 | 输出项列表 |
| `parameters_table` | array | 是 | 工具参数字段表 |
| `query_examples` | array[string] **或** array[object] | 是 | ≥1 条可复制 Prompt；**旗舰组学管线**使用三阶对象（见 §2.4） |
| `workflow_highlights` | array[string] | 否 | 描述区下方能力标签（如 Harmony、MACS2） |
| **`demo_visualization`** | **string** | **否** | **输出效果预览：图片 URL 或轻量 HTML 表格/流程图** |

**`demo_visualization` 约定：**

- 以 `http(s)://` 或 `/` 开头且扩展名为图片 → 渲染 `<img>`  
- 含 HTML 标签（如 `<table>`、`<div>`）→ 渲染到右侧「输出效果预览」栏（仅信任的后端静态数据）  
- **组学旗舰管线**：使用 `skill-viz-pipeline-gallery` 包裹多段 `<section>`，每段 `<h4 class="skill-viz-pipeline-stage__title">` + `<img src="/assets/images/demos/pipelines/...">` 或真实 TSV 摘要表；图源由 `scripts/build_omics_pipeline_demo_assets.py` 从 `data/results/run_*` 同步。  
- 示例 HTML 常量见 `gibh_agent/db/skill_detail_demo_visualizations.py` 与 `skill_detail_demo_visualizations_omics.py`

### 2.4 组学旗舰管线 · 结构化 `query_examples`（`skill_pipeline_specs.py`）

每条管线 **3 个对象**，字段：

| 字段 | 说明 |
|------|------|
| `tier` | `zero_shot` / `dynamic_routing` / `post_analysis` |
| `label` | 展示标题（**一键化标准执行** / 参数精调 / 深度挖掘；禁用非专业口语） |
| `hint` | 灰色说明条 |
| `prompt` | 一键复制正文 |

广场 **name → tool_id** 见 `SKILL_NAME_TO_PIPELINE_TOOL_ID`（如 `转录组学标准全流程` → `pipeline_transcriptomics`）。

### 2.2 标准样例（含可视化）

```json
{
  "tool_id": "blueprint_drafter",
  "description_long": "根据业务流程描述生成扁平工程蓝图 HTML…",
  "usage_hint": "用箭头描述步骤顺序；生成后在工作台 iframe 预览。",
  "demo_visualization": "<div>…RNA-seq 流程节点示意 HTML…</div>",
  "inputs": [ … ],
  "outputs": [ … ],
  "parameters_table": [ … ],
  "query_examples": [ "…" ]
}
```

### 2.3 Phase 1 已覆盖的 10 个 `tool_id`

见 §4 表格。

---

## 3. 前端抽屉 UI 规范

### 3.1 DOM 结构（`index.html`）

| 元素 ID | 类名 | 职责 |
|---------|------|------|
| `#skill-detail-drawer` | `.skill-detail-drawer` | 全屏遮罩 + 滑出面板容器 |
| `#skill-detail-drawer-back` | `.skill-detail-drawer__back` | 返回 / 关闭 |
| `#skill-detail-drawer-bookmark` | `.skill-detail-drawer__btn--ghost` | 收藏（复用 `toggleSkillBookmark`） |
| `#skill-detail-drawer-use` | `.skill-detail-drawer__btn--primary` | 立即使用（关闭抽屉 + `useSkillFromApi`） |
| `#skill-detail-drawer-main` | — | 描述 / I/O / 参数表 / 调用示例 |
| `#skill-detail-drawer-preview-col` | `.skill-detail-drawer__preview-col` | 输出效果预览（有 `demo_visualization` 时显示） |

打开态：根元素加 `.skill-detail-drawer--open`；`body.skill-detail-drawer-open` 禁止背景滚动。

### 3.2 JS 入口

| 函数 | 说明 |
|------|------|
| `openSkillDetailDrawer(skill)` | 打开并渲染 |
| `closeSkillDetailDrawer()` | 关闭 |
| `useSkillFromDetailDrawer()` | 使用并注入 Prompt |
| `renderSkillDetailVisualization(viz)` | 解析 `demo_visualization` |

---

## 4. Phase 1 十技能清单

| 技能名称 | 大类 | tool_id |
|----------|------|---------|
| 科研汇报 PPT 大纲生成 | 其他技能 | `ppt_outline` |
| 分析逻辑思维导图生成 | 其他技能 | `mindmap_gen` |
| 差异表达结果解读助手 | 生物医药 | `diff_expr_interpreter` |
| 周报撰写助手 | 其他技能 | `weekly_report_writer` |
| 深度调研 | 其他技能 | `deep_research` |
| 工程蓝图制图 | 其他技能 | `blueprint_drafter` |
| GO/KEGG 富集结果叙事生成 | 生物医药 | `go_kegg_narrative` |
| 单细胞实验设计检查清单 | 生物医药 | `single_cell_checklist` |
| 学术摘要精炼 | 其他技能 | `academic_abstract_refiner` |
| ChEMBL药物检索 | 化学 | `chembl_drug_search` |

---

## 5. 批次与批量生成

| 批次 | 文件 | 数量 | 说明 |
|------|------|------|------|
| Phase 1 | `skill_detailed_specs.py` | 10 | 手工精修 + 专属 `demo_visualization` |
| Phase 2 | `skill_detailed_specs_phase2.py` | 10 | 同上 |
| **Bulk** | `skill_detailed_spec_builder.py` → `skill_detailed_specs_bulk.py` | ~72 | **除 `main_category=多模态组学` 外**，凡 seed 中带 `[Skill_Route]` 且未在 Phase1/2 登记的技能 |
| **组学旗舰管线** | `skill_detail_demo_visualizations_omics.py` + `OMICS_PIPELINE_SPECS_BY_TOOL_ID` | 7 | `pipeline_transcriptomics` 等；**代码注册表优先于 DB** |

合并顺序：`Phase1` → `Phase2` → `BULK`（后者**不覆盖**前两批 `tool_id`）。组学管线由 `resolve_detailed_spec()` 与 `/api/skills/detailed-specs` **强制 merge** `OMICS_PIPELINE_SPECS_BY_TOOL_ID`，避免 DB 旧占位覆盖新 demo 图。

### 5.2 组学管线 Demo 图运维（成功经验）

更新「输出效果预览」真实截图时，按序执行：

```bash
cd /home/ubuntu/GIBH-AGENT-V2
# 1. 从 data/results/run_* 同步 PNG 到 nginx 静态目录
PYTHONPATH=. python3 scripts/build_omics_pipeline_demo_assets.py
# 2. 可选：回写 DB detailed_specs（与代码注册表双保险）
PYTHONPATH=. python3 scripts/patch_omics_pipeline_detailed_specs.py
# 3. 重启 API 使 Gunicorn 加载新 Python
docker compose -f docker-compose.yml restart api-server
# 4. 浏览器 Ctrl+Shift+R；验收 GET /api/skills/detailed-specs 含 sample1_R1 / manhattan_plot 等字段
PYTHONPATH=. python3 -m pytest tests/test_skill_detailed_specs.py -q
```

静态资源路径：`services/nginx/html/assets/images/demos/pipelines/{genomics,transcriptomics,...}/`

### 5.3 文案与参数要求

1. 调用示例须与 `launch_skill_demos.py` / `prompt_template` 同源。  
2. 有可视化交付的技能**必须**提供 `demo_visualization`（截图路径或示意 HTML）；批量项由 builder 按技能类型生成示意 HTML。  
3. 参数名与 `skill_*.py` 的 `execute()` 一致。  
4. 变更后跑：`PYTHONPATH=. python3 -m pytest tests/test_skill_detailed_specs.py -q`（含「非组学 routed 技能全覆盖」断言）。

---

## 6. 验收清单

- [ ] 抽屉从右侧滑出，宽 900px；含收藏 + 立即使用  
- [ ] `blueprint_drafter` 右侧可见流程蓝图示意  
- [ ] 组学旗舰管线（如 `pipeline_genomics`）右侧为真实 FASTQ/曼哈顿图等，非占位 HTML  
- [ ] 点击「立即使用」关闭抽屉并填入 Prompt  
- [ ] 浏览器硬刷新后样式与脚本生效  
