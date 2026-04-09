# 文件与 Omics 路径总线 — 与当前代码对齐

面向排障：数据从哪来、在哪被分类、与旧版 `resolve_omics_paths` 如何共存。架构原则见 `docs/ARCHITECTURE_CORE_PRINCIPLES.md` §1–§2。

---

## 1. 数据流（四层）

| 阶段 | 做什么 | 典型位置 / 形态 |
|------|--------|-----------------|
| **网关与落盘** | 环境变量 `UPLOAD_DIR`（默认 `/app/uploads`）、`RESULTS_DIR`（默认 `/app/results`）；相对路径拼到 `UPLOAD_DIR`；静态挂载 `/uploads`、`/results`（或等价）供浏览器访问 | `server.py`：`UPLOAD_DIR` / `RESULTS_DIR`、`validate_file_path`；`/api/chat` 内 `uploaded_files` 归一化 |
| **前端** | `uploaded_files`：`{ name, path }`（或 `file_name` / `file_path`）；工作流执行时常把多路径合并为**逗号串**写入步骤 `params` | `services/nginx/html/index.html`：`sendMessage`、执行分支里 `pathKeys` + `uniquePaths.join(', ')`（行号以仓库检索为准） |
| **Orchestrator（Path A）** | 多文件时 **`OmicsAssetManager.classify_assets`** 决定 **`path_for_inspect`（单一体检入口）** → `FileInspector.inspect_file` → `file_metadata`；并把嗅探摘要并入元数据供 Planner | `gibh_agent/core/orchestrator.py`：`path_for_inspect`、`_merge_routing_asset_fields`、`_classify_global_intent` 中的资产嗅探 |
| **Executor** | 对当前步骤的解析路径列表再 **`classify_assets` → `to_legacy_resolved_dict`**，得到与历史一致的六键 + **`fasta`**，并驱动**反射注入** | `gibh_agent/core/executor.py`：`OmicsAssetManager`、`to_legacy_resolved_dict`、`_inject_paths_from_asset_reflective` 等 |

---

## 2. 主总线：`OmicsAssetManager`（`gibh_agent/core/asset_manager.py`）

**定位**：将零散路径规范为 `DataAsset` 列表（**捆绑优先于单文件**），再通过 `to_legacy_resolved_dict` 与旧工具链使用的桶结构对齐。

### 2.1 `classify_assets` 固定顺序（与代码一致）

1. **`10x_bundle`**：`matrix.mtx` +（barcodes / features / genes）同批凑齐 → `primary_path` 为公共目录（规则与 `TenXPathResolver` 一致）。
2. **`spatial_bundle`**：Visium/空间矩阵（`.h5` 或 10x_mtx 目录）+ `spatial` 目录等（`_try_spatial_bundle`）。
3. **`radiomics_pair`**：影像与 mask 在未分配项中配对（同目录 1+1 或 stem 弱化匹配）。
4. **`metabolomics_csv`**：单文件表格（`_is_table_file`，含 `.tsv` 等）。
5. **`h5ad_single`**：仅 `.h5ad`。
6. **`10x_h5_matrix`**：单独出现的 10x 原生 `.h5` 矩阵（未与 spatial 成包时）。
7. **`protein_fasta`**：蛋白 FASTA。
8. **动态嗅探**：`sniff_asset_type_from_path` → `protein_structure` / `document` / `plain_text` / `generic_unknown` 等。

可选参数 **`user_intent_domain`** 时，对资产列表做**排序**（不改分类结果），用于与 `target_domain` 等快车道对齐。

### 2.2 `to_legacy_resolved_dict`（与 `resolve_omics_paths` 六键对齐 + 扩展）

返回字典含：**`tables, h5ad, 10x_mtx, images, masks, unknown`**（均为 `List[str]`），以及 **`fasta`**（蛋白序列下游）。

| `asset_type` | 写入桶（摘要） |
|--------------|----------------|
| `10x_bundle` | `10x_mtx` ← `primary_path`（目录） |
| `spatial_bundle` | 若主路径为 `.h5` 文件 → `h5ad`；否则 → `10x_mtx` |
| `10x_h5_matrix` | `h5ad` |
| `metabolomics_csv` | `tables` |
| `h5ad_single` | `h5ad` |
| `protein_fasta` | `fasta` |
| `radiomics_pair` | `images` ← 图；`masks` ← `metadata.mask_path` |
| `radiomics_image` / `radiomics_mask` | `images` / `masks` |
| `generic_unknown`、`protein_structure`、`document`、`plain_text` 等 | `unknown`（或约定类型进 `unknown`） |

**排障要点**：工具若仍直接调用 **`resolve_omics_paths`**，其 TenX/代谢/影像策略应与上述捆绑语义**一致**；若行为不一致，以 **`OmicsAssetManager` + `to_legacy_resolved_dict`** 为编排/执行主路径的准绳。

---

## 3. 兼容层：`resolve_omics_paths`（`gibh_agent/core/path_resolvers.py`）

- **仍存在**：策略链 + 门面，供未改写的工具或调试直接调用。
- **入参**：逗号/分号分隔的 `str` 或 `List[str]`。
- **返回**：六键 `tables, h5ad, 10x_mtx, images, masks, unknown`（**无** `fasta` 键；若需 FASTA 桶应走资产总线）。
- **re-export**：`gibh_agent/core/file_inspector.py` 导出 `resolve_omics_paths`（历史 import 路径仍可用）。

### 3.1 仍适用的硬规则（与总线 TenX 一致）

- **`.tsv` 表格规则**与 **TenX 候选**的先后顺序会导致：标准「三文件 `.tsv`」10x 往往**先进 `tables`**，`10x_mtx` 可能为空；与 **§2.1** 中「仅当文件名/组合满足 TenX 捆绑」才成 `10x_bundle` 的现象一致，需对照实际文件名与扩展名。
- **Radiomics**：`radiomics_tools.py` 等仍常消费 **`images[0]` / `masks[0]`**；与配对顺序、命名（mask/seg/label）相关。

---

## 4. Orchestrator：`path_for_inspect`（Path A，已以资产总线为主）

**文件**：`gibh_agent/core/orchestrator.py`（搜索 `path_for_inspect`）。

| 场景 | 行为 |
|------|------|
| **单文件** | 默认取首文件路径；若为 Spatial 相关单文件（如 `tissue_positions_list.csv`、路径段含 `spatial`），改为**父目录或 Visium 根**再 `inspect_file`，并可 `normalize_session_directory`。 |
| **多文件** | 解析为存在的绝对路径列表 → `_strip_replaced_archive_paths_static`、`_expand_visium_sibling_matrix_paths` → **`OmicsAssetManager().classify_assets(unique_paths)`**，再按**优先级**选**一个**体检根路径：`10x_bundle` → `spatial_bundle`（`visium_root` 或 `primary_path`）→ `radiomics_pair` → `metabolomics_csv` / `h5ad_single` / `10x_h5_matrix` → `radiomics_image` → 否则 `unique_paths[0]`。选中目录则 **`normalize_session_directory`** 后再 `inspect_file`。 |

**与旧文档差异**：多文件场景**不再**以「先 `resolve_omics_paths` 再看 `10x_mtx[0]`」为主路径；当前以 **`DataAsset` 类型优先级** 为准。

**仍可能出问题的点**：

- 多文件**未**形成 `10x_bundle` / `spatial_bundle` 时，会退回首路径或 `unique_paths[0]`，**commonpath** 语义弱化；目录不对则 Handler 仍可能失败。
- 影像 image/mask **跨目录**且未配成 `radiomics_pair` 时，mask 可能为空。
- `FileInspector` 对**目录**与**单文件**由不同 Handler 处理，见 `file_inspector` 与 `file_handlers/` 注册顺序。

**Planner 侧补充**：`_merge_routing_asset_fields` 将 `routing_asset_inventory` / `routing_asset_types` / `routing_asset_digest` 并入 `file_metadata`，供 LLM 规划读取（非仅原始 `files` 列表）。

---

## 5. Executor / Planner（与总线的关系）

- **Executor**：对 `resolved_file_paths` 做 **`classify_assets` → `to_legacy_resolved_dict`**，日志中的 `10x_mtx`/`h5ad`/`tables`/`images` 来自该字典；**`current_file_path`** 由 **`_seed_workflow_current_path(assets, omics_resolved, …)`** 等与资产类型联合决定；路径形参按**函数签名反射注入**（`_inject_paths_from_asset_reflective` 等），避免按工具名硬编码分支（宪法 §2）。
- **Planner**：主要消费 **`file_metadata`**（含上述 routing 字段），以及 Orchestrator 注入的摘要。
- **排障顺序**：`path_for_inspect` → `file_metadata` → Executor 日志中的 **OmicsAssetManager 桶** → 当前步 **`params`** 是否仍含未拆的逗号串。

---

## 6. 全局意图与资产（Chat vs Task）

**文件**：`orchestrator.py` 中 `_classify_global_intent`。

- 使用 **`ROUTING_WORKFLOW_NATIVE_ASSET_TYPES`** / **`ROUTING_NON_WORKFLOW_ASSET_TYPES`**（定义在 `asset_manager.py`）与上传嗅探结果，决定纯组学上传是否直接 **`task`**、蛋白/未知 + 模糊查询是否倾向 **`chat`** 等。
- 与 **快车道** `[Skill_Route:…]`、`target_domain` 等关系见 `docs/ARCHITECTURE_CORE_PRINCIPLES.md` §4。

---

## 7. 代码索引（按模块）

| 主题 | 位置 |
|------|------|
| **多模态资产总线（核心）** | `gibh_agent/core/asset_manager.py`：`DataAsset`、`OmicsAssetManager.classify_assets`、`to_legacy_resolved_dict`、`ROUTING_*` |
| 六键门面（兼容） | `gibh_agent/core/path_resolvers.py`：`resolve_omics_paths` |
| 总线 re-export | `gibh_agent/core/file_inspector.py` |
| 体检与 Handler | `gibh_agent/core/file_inspector.py`；Visium/扩展：`gibh_agent/core/file_handlers/extended_handlers.py` |
| Path A 体检路径（多文件资产优先级） | `gibh_agent/core/orchestrator.py`：`path_for_inspect`、`_merge_routing_asset_fields`、`_classify_global_intent` |
| 执行期分类与反射注入 | `gibh_agent/core/executor.py`：`OmicsAssetManager`、`to_legacy_resolved_dict`、`_seed_workflow_current_path` |
| 影像 resolved 消费 | `gibh_agent/tools/radiomics_tools.py`：`_radiomics_paths_from_resolved` 等 |
| 代谢 / RNA / 空间工具 | `gibh_agent/tools/metabolomics_tools.py`、`gibh_agent/tools/rna/*`、`gibh_agent/tools/spatial/*` |
| 上传目录与静态服务 | `server.py`：`UPLOAD_DIR`、`RESULTS_DIR`、`app.mount` |
| 前端多路径写入步骤 | `services/nginx/html/index.html`：`pathKeys`、`pathValue`（以符号搜索为准） |
| 架构白皮书 | `docs/ARCHITECTURE_CORE_PRINCIPLES.md` |

---

## 8. 建议的自检顺序（「丢文件 / 错模态」）

1. 打印 **`uploaded_files` / `file_paths`** 原始列表与是否已变为绝对路径。  
2. 对同一列表调用 **`OmicsAssetManager().classify_assets(paths)`**，核对 `asset_type` 与 **`to_legacy_resolved_dict`** 各桶。  
3. 对照 **`path_for_inspect`** 与多文件优先级（§4）。  
4. 查看 **`file_metadata`** 与 **`routing_asset_*`** 字段。  
5. 查看 Executor 日志 **「OmicsAssetManager: N 个资产」** 及当前步 **`params`**。  
6. 若某工具仍直接调 **`resolve_omics_paths`**，单独对其入参字符串做一次解析对照。

---

*与仓库代码同步维护；文中行号仅为历史参照，以符号搜索与 `docs/ARCHITECTURE_CORE_PRINCIPLES.md` 为准。*
