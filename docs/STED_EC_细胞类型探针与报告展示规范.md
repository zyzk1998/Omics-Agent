# STED-EC 细胞类型探针与报告展示规范

> **适用范围**：`STED_EC`（四步时空轨迹）与 `SPATIOTEMPORAL_DYNAMICS`（六步时空动力学）共用 `gibh_agent/tools/sted_ec_tools.py` 中的探针与自动标注闭环。  
> **实现文件**：`gibh_agent/tools/sted_ec_tools.py`、`gibh_agent/core/workflows/sted_ec_workflow.py`、前端 `services/nginx/html/index.html`。

---

## 1. 问题与目标

单细胞时空任务依赖 `adata.obs` 中的细胞类型列。用户数据可能：

- 已有规范或别名列（如 `celltype`、`CellType`）→ **不得**重写 h5ad、不得跑自动标注；
- 仅有 `leiden` / `louvain` / `clusters` 聚类编号 → **不**视为已注释，可走自动标注；
- 完全缺失 → 触发 Scanpy 兜底聚类 + marker 映射。

---

## 2. 后端：探针与分流

### 2.1 常量

| 常量 | 含义 |
| --- | --- |
| `STED_EC_CELL_TYPE_PROBE_KEYS` | 优先精确匹配：`cell_type`、`celltype`、`CellType`、`cell_types` |
| `STED_EC_PREANNOTATED_CELL_TYPE_KEYS` | 扩展预注释列名（含 `annotation`、`label`、`level*_…` 等） |
| `STED_EC_VIZ_ONLY_CLUSTER_KEYS` | `clusters`、`louvain`、`leiden` — **仅可视化回退，不算细胞类型** |

### 2.2 核心函数链

1. **`_probe_cell_type_key(adata)`** — 在 obs 中查找具生物学含义的标签列（`_has_meaningful_cell_type_labels` 过滤空值/Unknown）。
2. **`_resolve_sted_ec_cell_types(...)`** — 探针命中 → 返回原列名、`needs_persist=False`；未命中 → **`auto_annotate_cell_types_fallback`**。
3. **`_auto_annotate_cell_types_fallback_sync`** — HVG → PCA → neighbors → Leiden（失败回退 Louvain）→ **`rna/annotation._get_fallback_markers`** 映射 → 写 `cell_type` 与 `{stem}_auto_celltype.h5ad`；经 **`emit_tool_log` / `aemit_tool_log`** 推送 SSE 过程日志。

### 2.3 接入点

| 工具 / 步骤 | 行为 |
| --- | --- |
| `sted_ec_data_validation` | 第一步校验时调用 `_resolve_sted_ec_cell_types` |
| `sted_ec_time_series_formatting` | 二次探针兜底 |
| `sted_ec_preprocess` | 遗留入口，同样接入（正式 SOP 以工作流 DAG 为准） |

两条工作流（`STED_EC` / `SPATIOTEMPORAL_DYNAMICS`）第一步均为 **`sted_ec_data_validation`**，保证双通道一致。

### 2.4 依赖

无新增 pip 包；沿用 scanpy、anndata、leidenalg（随 scanpy）、moscot 等既有栈。见 `requirements.txt` 注释。

### 2.5 测试

`tests/test_sted_ec_auto_cell_type.py` — 覆盖探针命中跳过、缺失触发标注、聚类列不视为预注释等场景。

---

## 3. 前端：报告 summary 折叠规则

**原则**：步骤手风琴「查看详情」**默认全部展开**；仅对 STED-EC 流水线 summary 中的 **Obs 透视长 Markdown 表** 做折叠。

| 函数 | 作用 |
| --- | --- |
| `isStedEcPipelineStep` | 识别 `sted_ec_*` step_id 或 STED 相关中文步骤名 |
| `splitStedEcVerboseSummary` | 自「Obs 数据透视摘要」或宽 `| batch |` 表处拆分为 headline + verbose |
| `renderStepSummaryContentHtml` | headline 正常渲染；verbose 放入 `<details>`，默认折叠，摘要文案「展开 Obs 数据透视技术详情」 |
| `truncateExecStepSublabel` | 右栏 checklist 副标题仅显示 headline 截断，不铺整段 pipe 表 |

**禁止**：为「数据与元数据校验」单步单独默认折叠手风琴（已移除 `isStedEcMetadataValidationStep` 误用）。

---

## 4. 维护清单

- [ ] 修改探针键名或自动标注逻辑 → 同步本文件与 `sted_ec_workflow.py` 步骤描述。
- [ ] 修改 summary 文本格式（Obs 表标题）→ 同步 `splitStedEcVerboseSummary` 正则。
- [ ] 新增 STED-EC 步骤 → 更新 `isStedEcPipelineStep` 识别规则（若 step_id 不含 `sted_ec`）。
