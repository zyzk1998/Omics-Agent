# 文件与 Omics 路径总线 — 简版（对照当前代码）

面向排障：数据从哪来、在哪被拆开、哪里容易错。不展开重构方案。

---

## 1. 数据流（三层）

| 阶段 | 做什么 | 典型形态 |
|------|--------|----------|
| 前端 | `uploaded_files` 数组 `{ name/path }`；工作流执行时多路径常合并为**一条逗号字符串**写入多个 `params` 键 | `services/nginx/html/index.html`：`sendMessage` 的 `uploaded_files`；执行逻辑里 `pathKeys` + `uniquePaths.join(', ')`（约 L5297+） |
| Orchestrator | 规范化列表 → Path A 里决定 **`path_for_inspect`** → **`FileInspector.inspect_file`** → `file_metadata` | `gibh_agent/core/orchestrator.py`：`path_for_inspect` 与 `inspect_file` |
| 工具 / Executor | 逗号串或单路径再经 **`resolve_omics_paths`** 分到 `tables` / `h5ad` / `10x_mtx` / `images` / `masks` / `unknown` | `gibh_agent/tools/*` 见下节索引 |

---

## 2. `resolve_omics_paths`（当前实现）

- **实现位置**：`gibh_agent/core/path_resolvers.py`（策略类 + 门面）。
- **对外入口**：`gibh_agent/core/file_inspector.py` 中 `from .path_resolvers import resolve_omics_paths`（旧 `import file_inspector.resolve_omics_paths` 仍有效）。
- **入参**：`str`（逗号/分号分隔）或 `List[str]`。
- **返回**（键名与顺序语义与历史一致，下游依赖此结构）：

```text
tables, h5ad, 10x_mtx, images, masks, unknown  → 均为 List[str]
```

- **策略顺序**（与旧 if-else 链一致）：`Metabolomics`（表）→ `AnnData`（h5ad/h5）→ `TenX`（候选聚合）→ `Radiomics`（影像/掩膜）→ 门面再把未分配路径丢进 `unknown`。

### 2.1 排障时优先核对的「硬规则」

1. **`.tsv` 先于 10x 候选**  
   `barcodes.tsv` / `features.tsv` 会因后缀进 **`tables`**，**不会**进入 TenX 候选。因此标准「三文件 .tsv」上传时，**`10x_mtx` 经常为空**；10x 聚合只会在「未被表规则吃掉的文件名」上触发（例如无扩展名但名中含 `barcodes`/`genes` 等，或仅 `matrix.mtx` 与其他非 `.tsv` 规则文件组合凑齐条件）。  
   **现象**：工具里 `10x_mtx` 空、路径落在 `tables` + `unknown`。先对照 `path_resolvers.TenXPathResolver` 与 `MetabolomicsPathResolver` 的先后顺序。

2. **影像 `images[0]` / `masks[0]`**  
   `gibh_agent/tools/radiomics_tools.py` 的 `_radiomics_paths_from_resolved` 仍默认 **`images[0]`、`masks[0]`**，顺序与命名（mask/seg/label 等）敏感；与 `RadiomicsPathResolver` 的关键词规则一致。

3. **单参数字符串**  
   工具若只收到 `"a,b,c"`，必须先 `resolve_omics_paths` 再取对应桶；若某处只截 **第一个路径** 再 resolve，多文件语义会丢。

---

## 3. Orchestrator：`path_for_inspect`（Path A）

**文件**：`gibh_agent/core/orchestrator.py`（搜索 `path_for_inspect`）。

| 场景 | 行为 |
|------|------|
| 单文件 | 默认首文件路径；若命中 Spatial 相关单文件（如 `tissue_positions_list.csv` / 路径含 `spatial`），改为**父目录或再上级的目录**再 `inspect_file`。 |
| 多文件 | 先解析出存在的 `resolved_paths`；**若 `len > 1`**：① **`resolve_omics_paths(各路径)`** — 若 **`10x_mtx` 非空**，则 **`path_for_inspect = 10x_mtx[0]`**（目录），并 `normalize_session_directory`；② 否则若**全部为** `.nii.gz`/`.nii`/`.dcm`，则**选一个非 mask/label 文件**单文件体检（影像）；③ 否则 **`commonpath` 会话目录** + `normalize_session_directory`。 |

**仍可能出问题的点**：

- 多文件 10x **未**进 `10x_mtx` 时，会退回 `commonpath`；若 common 根不对，TenX 目录 Handler 仍可能失败。  
- 影像若 image/mask **不在同一目录**，`RadiomicsHandler` 依赖 `pair_radiomics_files(search_dir)`，mask 可能一直为空。  
- `FileInspector` 需支持**目录**（如 TenX、Spatial、解压后的会话目录）；单文件与目录由不同 Handler 处理，见 `file_inspector` 注册优先级。

---

## 4. Executor / Planner（与总线的关系）

- **Executor**：`current_file_path` 默认取 **`file_paths[0]`**；Radiomics 等依赖步骤里已填的 **`image_path` / `mask_path`**，不会从列表自动拆双路径。  
- **Planner**：主要消费 **`file_metadata`**（来自上面那一次 `inspect_file`），不是原始 `files` 全量列表。  
**排障**：规划错先看 `path_for_inspect` 是否选对；再看 `file_metadata` 是否含 `mask_path` / `real_data_path` / `file_type`。

---

## 5. 代码索引（按模块）

| 主题 | 位置 |
|------|------|
| 路径总线实现 | `gibh_agent/core/path_resolvers.py` |
| 总线 re-export | `gibh_agent/core/file_inspector.py` |
| 体检与 Handler | `gibh_agent/core/file_inspector.py`；影像/Visium：`gibh_agent/core/file_handlers/extended_handlers.py` |
| Path A 体检路径 | `gibh_agent/core/orchestrator.py`：`path_for_inspect`、`resolve_omics_paths` |
| 影像 resolved 消费 | `gibh_agent/tools/radiomics_tools.py`：`_radiomics_paths_from_resolved` |
| 代谢 tables[0] | `gibh_agent/tools/metabolomics_tools.py` |
| RNA 10x / h5ad | `gibh_agent/tools/rna/quality_control.py`：`_resolve_adata_path` |
| 空间 h5ad / 10x | `gibh_agent/tools/spatial/analysis.py` |
| 前端多路径写入步骤 | `services/nginx/html/index.html`：`pathKeys`、`pathValue`（约 L5297、L6577 等，以实际文件为准） |

---

## 6. 建议的自检顺序（出现「丢文件 / 错模态」时）

1. 打印 **`uploaded_files` 或 `file_paths`** 原始列表是否齐全。  
2. 对同一列表调用 **`resolve_omics_paths`**，看六个桶的划分是否符合第 2.1 节预期。  
3. 看 Orchestrator **`path_for_inspect`** 最终是文件还是目录、是否走了 **`10x_mtx` 分支**。  
4. 看 **`file_metadata`** 成功字段与 `file_type` / `mask_path`。  
5. 再看 **Executor 当前步 `params`** 是否仍带逗号串、工具是否只读了第一段。

---

*文档与仓库代码同步；行号仅作近似定位，以符号搜索为准。*
