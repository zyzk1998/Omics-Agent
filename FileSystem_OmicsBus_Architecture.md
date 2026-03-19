# 通用多文件总线系统 (Omics Path Resolution) — 深度技术报告

> 面向首席系统架构师：文件上传、路径去重、resolve_omics_paths 总线、Planner/Executor 文件处理与单/多文件死循环根因分析。

---

## 1. 前端文件收集与路由层 (Frontend File Bus)

### 1.1 文件资产收集方式

**文件路径**: `services/nginx/html/index.html`

- **全局状态**: `let uploadedFiles = [];`（约 L2495），元素形态为 `{ name, path }` 或 `{ file_name, file_path }`（与后端返回一致）。
- **来源**:
  1. **POST /api/upload** 成功后，响应 `data.file_paths` 或 `data.file_path` / `data.files` 被转成上述结构并 push 到 `uploadedFiles`（见 L4709–4766 附近）。
  2. **左侧数据资产** 点击后 `addAssetAsAttachment(a.file_path, a.file_name)`，将资产路径加入当前“待发送”上下文（不一定会改写全局 `uploadedFiles`，取决于实现）。
  3. **工作流卡片执行** 时从表单参数、`uploadedFiles`、`.asset-tag`、`workflowCache`、当前消息的 `.file-chip` 等多源汇总路径（L6711–6743）。

单文件与多文件在数据结构上**无区分**：均为「对象数组」，每一项为 `{ name/file_name, path/file_path }`。

### 1.2 向后端发送请求时的 Payload 形态

**聊天/执行入口**: `sendMessage` 构建 payload（约 L5014–5025）：

```javascript
const payload = {
    message: finalMessage,
    history: chatHistory,
    uploaded_files: filesToSend,  // 即 [{ name, path }, ...]
    workflow_data: workflowData,
    stream: true,
    session_id: currentSessionId || undefined
};
```

- **单文件（如代谢组 .csv）**: `uploaded_files: [{ name: "data.csv", path: "/app/uploads/owner_id/batch/data.csv" }]`。
- **多文件（如影像组 image.nii.gz + mask.nii.gz，或 10x 三件套）**: `uploaded_files: [{ name: "a.nii.gz", path: "..." }, { name: "b_mask.nii.gz", path: "..." }]`，或 10x 多个 `{ name, path }`。

即：**单/多文件共用同一数组结构**，仅长度不同；后端通过 `files` 列表 + 后续解析区分类型。

**工作流执行**（`executeWorkflowFromForm`，L6752–6777）：

- 从多源收集路径 → `uniquePaths = [...new Set(filePaths.filter(Boolean))]`。
- **所有路径类参数统一赋值为同一字符串**：`pathValue = uniquePaths.join(', ')`，并对以下 key 写回每一步的 `params`：
  - `pathKeys = ['file_path', 'adata_path', 'data_path', 'data_path_2', 'input_path']`。
- 即：**多文件被压成一条逗号分隔字符串**，写入 `file_path`、`adata_path`、`data_path` 等，**不区分组学类型**，也不区分 image vs mask。

```javascript
// index.html L6754-6762
var pathKeys = ['file_path', 'adata_path', 'data_path', 'data_path_2', 'input_path'];
var uniquePaths = [...new Set(filePaths.filter(Boolean))];
var pathValue = uniquePaths.length > 0 ? uniquePaths.join(', ') : '';
for (var si = 0; si < updatedSteps.length; si++) {
    var sp = updatedSteps[si].params || {};
    for (var ki = 0; ki < pathKeys.length; ki++) {
        if (Object.prototype.hasOwnProperty.call(sp, pathKeys[ki]) && pathValue) {
            updatedSteps[si].params[pathKeys[ki]] = pathValue;
        }
    }
}
```

- **workflow_data** 中同时带 `file_paths: uniquePaths`（数组），供后端合并/去重使用。

### 1.3 路径去重逻辑

- **前端**: 仅 `uniquePaths = [...new Set(filePaths.filter(Boolean))]`，按引用/字符串去重，**无按模态或角色（image/mask）的去重**。
- **后端**（orchestrator）: `file_paths = list(dict.fromkeys([p for p in file_paths if p]))`（L418、L1762 等），顺序保留、去重。同样**不区分组学与语义**。

**结论**: 前端与后端都只有「扁平路径列表 + 逗号拼接」这一套逻辑；影像组 image/mask、单细胞 10x 多文件、代谢组单文件共用同一通道，差异全靠后端 `resolve_omics_paths` 与各工具对「单路字符串」的再解析，易出现「修单文件坏多文件、修多文件坏单文件」的互相干扰。

---

## 2. 后端解析总线 (Backend `resolve_omics_paths`)

### 2.1 函数位置与签名

**文件**: `gibh_agent/core/file_inspector.py`  
**入口**: `resolve_omics_paths(raw_path_string: str) -> Dict[str, List[str]]`

- **入参**: 单条字符串，可含逗号/分号分隔的多个路径（如 `"path1, path2"`）。
- **返回**:  
  `{ "tables": [], "h5ad": [], "10x_mtx": [], "images": [], "masks": [], "unknown": [] }`  
  各 key 为路径列表。

### 2.2 内部逻辑（按执行顺序）

1. **标准化与拆分为列表**  
   - `raw.replace(";", ",").strip()` → `split(",")` → 每段 strip，非空保留。  
   - 对每段做 `Path(p).resolve()`（存在则用绝对路径），否则保留原串，得到 `paths` 列表。

2. **逐路径分类（单路 if-else 链）**  
   对每个 `p` 用 `path_obj.name` / `suffix` 判断，**命中即 continue**，未命中最后进 `unknown`：

   - **表格**: `suffix` 或 name 为 `.csv`/`.tsv`/`.xlsx`/`.txt` 或 `.tsv.gz` → `out["tables"].append(p)`。
   - **h5ad**: `.h5ad` 或 `.h5` → `out["h5ad"].append(p)`。
   - **10x 相关**: 文件名含 `matrix.mtx`、`barcodes`+tsv、`features`/`genes`+tsv → 只加入 `tenx_paths`，**不立刻写入 out**。
   - **影像**: 扩展名为 `.nii`/`.nii.gz`/`.dcm`/`.tiff`/… 或 name 含 `.nii`；再按文件名是否含 `mask`/`roi`/`seg`/… 且不含 `image`/`img`/`mri`/`ct`/… → `out["masks"]` 或 `out["images"]`。
   - 其余 → `out["unknown"].append(p)`。

3. **10x_mtx 后处理**  
   - 若 `tenx_paths` 非空：检查是否存在 `matrix.mtx` 且存在 `barcodes` 或 `features`/`genes`；  
   - 若成立：`out["10x_mtx"] = [os.path.commonpath(dirs)]`（单元素：公共父目录）；  
   - 否则：将 `tenx_paths` 全部并入 `out["unknown"]`。

### 2.3 痛点暴露：「坏味道」代码与单/多文件互相干涉

**（1）单路字符串 + 统一逗号拆分**

- 所有组学共用「一条字符串、按逗号拆」：代谢组单文件为 `"path/to/data.csv"`，影像组多文件为 `"image.nii.gz, mask.nii.gz"`，10x 为多个路径或一个目录。
- 若某处只传「第一个路径」或只传「一个路径」给该总线（例如规划阶段只带 `file_path` 单值），则 10x 会得不到完整三件套、影像组可能只得到 image 或只得到 mask，**单文件修复（只传一条）会破坏多文件；多文件修复（传整串）可能破坏期望「单路径」的代谢组逻辑**。

**（2）影像 image/mask 仅靠文件名关键词**

- image vs mask 的划分完全依赖文件名是否含 `mask`/`roi`/`seg` 等（且排除 `image`/`img`/`mri`/`ct` 等）。  
- 命名不规范（如 `seg.nii.gz` 与 `t1.nii.gz`）或顺序变化（前端先传 mask 再传 image）会导致分类与下游期望不一致；且 **images/masks 之间无显式配对关系**，仅按列表顺序取 `images[0]`、`masks[0]`（见下节工具层）。

**（3）10x 的「公共父目录」与单路径假设冲突**

- 10x 正确性依赖：多路径的**公共父目录**存在且包含 matrix.mtx + barcodes/features。  
- 若上游只传一个路径（如只传 `matrix.mtx` 的路径），则 `tenx_paths` 只有一条，`dirs = [os.path.dirname(p)]`，`commonpath` 仍是该目录，可能仍能识别；但若上游误传成「多个不共父目录的路径」或只传 barcodes 不传 matrix，则 10x_mtx 为空、全部落入 unknown，**单文件与多文件共用一个入口导致边界情况复杂**。

**（4）工具层对 resolve 结果的二次假设（影像组）**

`gibh_agent/tools/radiomics_tools.py` 中 `_radiomics_paths_from_resolved`：

```python
# radiomics_tools.py L17-38
def _radiomics_paths_from_resolved(resolved: Dict[str, List[str]]) -> Tuple[List[str], Optional[str], Optional[str]]:
    images = resolved.get("images") or []
    masks = resolved.get("masks") or []
    tables = resolved.get("tables") or []
    unknown = resolved.get("unknown") or []
    all_paths: List[str] = images + masks + tables + unknown
    image_path: Optional[str] = images[0] if images else None
    mask_path: Optional[str] = masks[0] if masks else None
    if tables and not images and not masks:
        return all_paths, tables[0], None
    if images and not masks and tables:
        return all_paths, tables[0], None
    if not image_path and unknown:
        image_path = unknown[0]
    if not mask_path and len(unknown) >= 2:
        mask_path = unknown[1]
    return all_paths, image_path, mask_path
```

- 影像组工具**假定**：`images[0]` = 原图，`masks[0]` = 掩膜；无 image 时用 `unknown[0]`，无 mask 时用 `unknown[1]`。  
- 一旦 resolve 阶段把 image 和 mask 分错（或只收到一条路径），这里就会错配；且「tables 有值且 images/masks 无」时直接退回 CSV 路径，**单文件（CSV）与多文件（image+mask）共用一个出口**，分支纠缠。

---

## 3. 规划与执行层的上下文传递 (Planner & Executor Context)

### 3.1 Planning 阶段：文件路径在 Prompt 中的形态

- **Orchestrator** 在「Path A」中：从 `files` 取**第一个**元素得到 `file_path`（单值），再根据 `len(files)==1` 或 `len(files)>1` 决定 `path_for_inspect`（见下）；**仅将该单路径或目录传给 Inspector**。
- **Planner** 接收的是 **file_metadata**（来自 `file_inspector.inspect_file(path_for_inspect)`）以及**单一 `file_path`**（或 real_data_path）。  
- 规划阶段 LLM 看到的「文件路径」来自 file_metadata 的 `file_path`、`real_data_path`、`mask_path`（Radiomics）等字段，**并非完整的 `file_paths` 数组**；工作流配置里 `file_paths` 常被设为 `[file_path]`（单元素），例如：

```python
# planner.py L1983, L2481, L451
"file_paths": [file_path] if file_path else []
# 或
"file_paths": context_files  # context_files 来自上游，可能仍是单路径列表
```

因此 **Planning 阶段默认是「单路径/单目录」视角**；多文件信息仅通过 Inspector 在「单路径/目录」下自发现（如 Radiomics 在同目录找 mask、10x 在目录内找 matrix+barcodes+features）。

### 3.2 多文件时 path_for_inspect 的确定（单/多分支干涉点）

**文件**: `gibh_agent/core/orchestrator.py`（约 L1343–1419）

- **单文件** (`len(files)==1`):  
  - `path_for_inspect = file_path`；  
  - 若是 Spatial 相关（如 `tissue_positions_list.csv` 或路径含 `spatial`），则改为其父目录或更上层目录再 inspect。  
- **多文件** (`len(files)>1`):  
  - 先解析所有路径为绝对路径列表 `resolved_paths`；  
  - 若**全部为影像扩展名**（`.nii.gz`/`.nii`/`.dcm`）：  
    - **选一个「非 mask/非 label」的候选**作为 `path_for_inspect`（单文件体检）；  
  - 否则：取 `common = os.path.commonpath(...)`，对 `common` 做目录规范化后 `path_for_inspect = str(common)`（按目录体检）。

```python
# orchestrator.py L1390-1416
if all_radiomics:
    candidates = [p for p in resolved_paths if "mask" not in p.name.lower() and "label" not in p.name.lower()]
    path_for_inspect = str(candidates[0] if candidates else resolved_paths[0])
else:
    # ...
    path_for_inspect = str(common)
```

- **问题**：Radiomics 多文件时只传**一个**路径给 Inspector；Inspector（RadiomicsHandler）在该文件的**父目录**内用 `pair_radiomics_files(search_dir)` 找 image+mask。若 image 与 mask 不在同一目录（例如来自两次上传、不同 session 目录），则 **mask 永远找不到**，`file_metadata["mask_path"]` 为 None，后续 Planner/Executor 的 mask 缺失，导致单/多文件行为不一致。

### 3.3 Execution 阶段：底层工具如何接收单/多文件路径

**Executor**（`gibh_agent/core/executor.py`）：

- **输入**: `execute_workflow(workflow_data, file_paths)`，`file_paths` 为字符串列表（已去重、绝对路径）。
- **上下文**: `current_file_path = resolved_file_paths[0] if resolved_file_paths else None`（L1216），即**只保留「第一个」解析后的路径**作为后续步骤的默认输入。
- **参数映射**（L308–324、L785–798）：  
  - **scRNA-seq**: `file_param_name = "adata_path"`；若有 `file_path` 则映射为 `adata_path`。  
  - **Radiomics**: `file_param_name = "image_path"`；**不**自动注入 `file_path`，依赖步骤 params 内已有 `image_path` / `mask_path`。  
  - **其他（如代谢组）**: `file_param_name = "file_path"`，缺失时用 `current_file_path` 注入。

因此：

- **代谢组（单文件）**: 通常一个 `file_path`，Executor 用 `current_file_path` 或步骤内 `file_path` 即可；若某次修改为「多路径拼接」且工具只读第一个，可能被截断或误用。  
- **影像组（多文件）**: **不**从 `file_paths` 列表推导 image/mask，而是依赖 **Planner 在步骤 params 中写好的 `image_path` / `mask_path`**（来自 file_metadata）；若 Inspector 未返回 mask（例如跨目录），执行阶段就没有 mask。  
- **单细胞 10x（多文件/目录）**: 工具层期望的是**目录路径**或 **逗号拼接路径**，经 `resolve_omics_paths` 得到 `10x_mtx[0]` 目录；若上游只传了单文件路径，则 10x_mtx 可能为空，工具报错或回退到 unknown。

**参数映射与占位符**：  
Executor 的 `_process_data_flow` 会替换占位符（如 `<preprocess_data_output>`），并只对签名内参数做 path 解析；Radiomics 会删除 `file_path`/`output_dir` 等非签名参数，避免误注入。**image_path/mask_path 必须在 step.params 中已存在**，否则不会从 `file_paths` 自动拆成 image+mask。

### 3.4 小结：单/多文件干涉的三大代码点

| 层级 | 单文件假设 | 多文件假设 | 冲突点 |
|------|------------|------------|--------|
| 前端 | 同构数组，可只有 1 个元素 | 多元素数组 → 逗号拼接进所有 path 参数 | 同一 pathValue 写 file_path/adata_path/data_path，影像组需要 image/mask 两个语义 |
| resolve_omics_paths | 单路字符串也可 | 逗号分隔多路，按后缀/名称分类 | 10x 依赖多路或目录；影像依赖 images[0]/masks[0]，顺序/命名敏感 |
| Orchestrator / Planner | file_path 单值，path_for_inspect 单路径 | 多文件时用 commonpath 或「一个候选」 | Radiomics 只传一条路径给 Inspector，mask 依赖同目录发现 |
| Executor | current_file_path = file_paths[0] | Radiomics 不从这里拆 image/mask | 影像组完全依赖 Planner 填好的 image_path/mask_path |

---

## 4. 架构师视角的重构建议 (Refactoring Strategy)

### 4.1 缺失的抽象接口（策略 / 适配器）

当前系统**缺少**：

1. **按模态区分的「路径解析策略」接口**  
   - 应有 `OmicsPathResolver` 抽象（或协议），按域实现：  
     - `MetabolomicsResolver`: 单表路径 → 单 path；  
     - `RadiomicsResolver`: 多路径/目录 → `(image_path, mask_path)` 或 `(image_path, mask_path, tables)`；  
     - `TenXResolver`: 多路径或目录 → 单目录路径（或 10x 三件套结构）。  
   - 现有 `resolve_omics_paths` 是「一个大 if-else 函数」，应拆成上述策略的**门面**，由调用方根据 domain 或文件特征选择策略，避免单/多文件逻辑在同一分支里互相影响。

2. **「文件集 → 执行参数」的适配器**  
   - 应有 `FileSetToToolParams` 适配器（或按工具注册）：  
     - 输入：`(domain, file_paths[], file_metadata?)`；  
     - 输出：该域工具所需的参数字典（如 Radiomics → `{ image_path, mask_path }`，Metabolomics → `{ file_path }`，10x → `{ adata_path }` 或 `{ matrix_dir }`）。  
   - 这样前端/Orchestrator 只需传「路径列表 + 域」，不再写死 `pathValue = uniquePaths.join(", ")` 并塞进所有 path key。

3. **校验与组装的域专属接口**  
   - **Validation**: 应有 `ValidateFileSet(domain, paths) -> Result`，例如：  
     - Radiomics: 至少 1 个 image，可选 1 个 mask，或 1 个 CSV（已提取特征）；  
     - 10x: 目录内存在 matrix.mtx + barcodes/features，或单 h5ad；  
     - Metabolomics: 单 CSV/TSV。  
   - **Assembly**: 应有 `AssembleToolInput(domain, paths, metadata?) -> dict`，把路径列表 + 可选 metadata 组装成该域工具的参数结构，而不是在 Executor 里用 `tool_category` + 一堆 if 分支。

### 4.2 三种组学的致命缺陷归纳

| 组学 | 校验 (Validation) 缺陷 | 组装 (Assembly) 缺陷 |
|------|------------------------|----------------------|
| **影像组 (Image+Mask)** | 1）校验分散在 `radiomics_data_validation` 与 `resolve_omics_paths`，无统一「成对」校验；2）单文件影像必成对否则抛错，多文件依赖同目录发现 mask，跨目录上传即失败。3）image/mask 仅靠文件名关键词，易误判。 | 1）Planner 从 file_metadata 取 mask_path，若 inspect 只看了单路径所在目录则 mask_path 恒为 None。2）前端把所有路径拼成一条字符串赋给所有 path 参数，无法表达「image vs mask」语义。3）Executor 不根据 file_paths 组装 image_path/mask_path，完全依赖步骤 params 已填好。 |
| **单细胞 (10x 多文件/目录)** | 1）10x 有效性与「公共父目录 + matrix+barcodes+features」强绑定，若只传一个文件路径可能进 unknown 或 10x_mtx 空。2）单文件 matrix.mtx 被单独上传时，校验与解析逻辑不一致（有时当普通文件处理）。 | 1）工具期望 adata_path 或 10x 目录；若上游只传了 `file_paths[0]`（某个子文件），resolve 后可能得不到 10x_mtx。2）Executor 的 current_file_path 仅取第一个路径，10x 需目录时若列表里是多个文件路径，需依赖 resolve 得到目录并写回 params，当前依赖 Planner/Inspector 的 real_data_path。 |
| **代谢组 (单文件)** | 1）单 CSV 校验清晰，但若前端误传多路径（如多 CSV），resolve 会得到 tables 列表，工具层通常只取 tables[0]，其余被忽略，无「多表」语义。2）与影像组共用同一 pathValue 写入逻辑，易被多文件改动误伤。 | 1）Assembly 简单（单 file_path），但若为「修多文件」而改成总是传 file_paths 数组或逗号串，必须保证代谢组工具只取第一个或明确约定「多表」含义，否则单文件场景会不稳定。 |

### 4.3 重构动作建议（简要）

1. **引入 OmicsPathResolver 策略**  
   - 为 Metabolomics / Radiomics / TenX（及 Spatial）各实现一个 resolver，输入为「原始路径字符串或路径列表」，输出为该域的结构化结果（单 path / (image, mask) / directory 等）；`resolve_omics_paths` 改为根据调用方传入的 domain 或自动检测选用对应策略，避免在一个函数内用同一套 if-else 处理所有模态。

2. **Orchestrator 多文件 Radiomics**  
   - 多文件影像时，向 Inspector 传入「目录或全部路径」而非单一路径；或在 Inspector 外先调用 RadiomicsResolver，得到 image_path + mask_path 后再决定 inspect 哪一条（或两条都传），保证 file_metadata 中 mask_path 在跨目录上传时也能被正确填充（或显式校验「必须同目录」并提示用户）。

3. **前端 path 写入与 workflow 执行**  
   - 执行工作流时，不要用同一 `pathValue` 覆盖所有 path key；按**域或步骤类型**决定写入方式：  
     - Radiomics 步骤：只写 `image_path` / `mask_path`（可从后端下发的 workflow 中带出，或由前端按命名规则拆成两条）；  
     - 代谢组：只写 `file_path`（单路径或约定取第一个）；  
     - 10x/RNA：写 `adata_path` 或目录路径。  
   - 可后端在 workflow 下发时带 `path_schema: { image_path: [...], mask_path: [...] }` 或类似，前端按 schema 装配，避免盲目 join 成一条。

4. **Executor 参数注入**  
   - 用「FileSetToToolParams」适配器替代当前按 `tool_category` 的 if-else：根据 domain + file_paths + 已有 file_metadata 生成该步骤的路径参数字典，再与 LLM/Planner 填写的其他参数合并，这样 Radiomics 的 image_path/mask_path 可由适配器统一从 file_paths 或 metadata 组装，不依赖 Planner 一定填对。

5. **校验前置与统一**  
   - 在 Orchestrator 或独立层对「当前请求的 files + domain」做一次 ValidateFileSet；Radiomics 要求同目录 image+mask 或明确拒绝跨目录并提示；10x 要求目录或完整三件套路径；Metabolomics 要求单表。校验失败早退，避免进入规划/执行后再在工具层报错，从而减少「修一处、崩另一处」的回归。

---

## 5. 关键代码索引（便于重构定位）

| 模块 | 文件路径 | 行号/函数或片段 |
|------|----------|------------------|
| 前端上传结果处理 | `services/nginx/html/index.html` | L4709–4766 `file_paths` / `uploadedFiles` |
| 前端执行时路径汇总与写入 | `services/nginx/html/index.html` | L6711–6777 `executeWorkflowFromForm`，pathKeys / pathValue |
| 解析总线 | `gibh_agent/core/file_inspector.py` | L29–134 `resolve_omics_paths` |
| 影像从 resolved 取路径 | `gibh_agent/tools/radiomics_tools.py` | L17–38 `_radiomics_paths_from_resolved` |
| 代谢组路径解析 | `gibh_agent/tools/metabolomics_tools.py` | L41–46 `resolve_omics_paths(str(path_in))`，tables[0] |
| RNA 路径解析 | `gibh_agent/tools/rna/quality_control.py` | L20–32 `_resolve_adata_path`，10x_mtx / h5ad |
| Orchestrator 单/多 path_for_inspect | `gibh_agent/core/orchestrator.py` | L1343–1419，单文件 Spatial / 多文件 all_radiomics / commonpath |
| Orchestrator file_paths 去重与绝对化 | `gibh_agent/core/orchestrator.py` | L414–427，L1761–1767 |
| Planner file_paths 单值 | `gibh_agent/core/planner.py` | L1983，L2481，L451 |
| Planner Radiomics image_path/mask_path | `gibh_agent/core/planner.py` | L1607–1611 `image_path` / `mask_path` from file_metadata |
| Executor current_file_path | `gibh_agent/core/executor.py` | L1206–1216 |
| Executor 按类别 file_param_name | `gibh_agent/core/executor.py` | L308–324，L785–798，L388–403 |
| RadiomicsHandler 同目录配对 | `gibh_agent/core/file_handlers/extended_handlers.py` | L207–258，`pair_radiomics_files(search_dir)` |
| pair_radiomics_files | `gibh_agent/core/file_handlers/structure_normalizer.py` | L230–258 |

---

以上为当前「通用多文件总线」与 Omics 路径解析的架构现状、单/多文件互相干涉根因及重构方向，供架构师决策与拆分任务使用。
