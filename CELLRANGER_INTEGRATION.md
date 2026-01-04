# Cell Ranger 集成到智能体 - 完成总结

## ✅ 已完成的工作

### 1. CellRangerTool 增强
- ✅ 添加 `run_count()` 方法，支持实际执行 Cell Ranger count
- ✅ 支持 Cell Ranger 10.0.0+ 的 `--create-bam` 参数
- ✅ 自动检测 Cell Ranger 安装路径（支持 bin/cellranger 和直接路径）
- ✅ 完整的错误处理和状态返回

### 2. ScanpyTool 扩展
- ✅ 添加 `run_cellranger()` 方法，封装 Cell Ranger 调用
- ✅ 添加 `convert_cellranger_to_h5ad()` 方法，自动转换输出格式
- ✅ 在工具映射表中注册新工具

### 3. RNAAgent 更新
- ✅ 支持从 FASTQ 文件开始的分析流程
- ✅ 自动执行：Cell Ranger → 转换 → Scanpy 分析
- ✅ 更新工具描述，让 Agent 知道可以使用 Cell Ranger
- ✅ 更新 `execute_workflow()` 方法，支持完整流程

### 4. 系统提示词更新
- ✅ 添加 Cell Ranger 工具说明
- ✅ 更新工作流规则，支持从 FASTQ 开始
- ✅ 明确工具使用场景和参数

### 5. 文件类型检测增强
- ✅ 更新 `detect_file_type()` 方法，支持识别 FASTQ 目录
- ✅ 支持识别 10x MTX 目录和 Cell Ranger 输出目录

### 6. 配置更新
- ✅ 更新 `settings.yaml`，添加 Cell Ranger 配置项
- ✅ 支持环境变量配置（`CELLRANGER_PATH`, `CELLRANGER_REF`）

---

## 🚀 使用方法

### 方式 1: 通过智能体交互（推荐）

**提示词示例：**
```
我需要从 FASTQ 文件开始分析单细胞数据。

FASTQ 文件路径：/path/to/fastqs/
参考基因组路径：/path/to/refdata-gex-GRCh38-2024-A/
样本名称：sample_name

请帮我运行 Cell Ranger count，然后将输出转换为 Scanpy 格式，并执行完整的分析流程。
```

**或者更简洁：**
```
分析 /data/fastqs/ 目录下的数据，参考基因组在 /data/refdata/
```

### 方式 2: 直接调用工具

```python
from gibh_agent.tools.cellranger_tool import CellRangerTool
from gibh_agent.tools.scanpy_tool import ScanpyTool

# 初始化工具
cellranger_tool = CellRangerTool({
    "path": "/home/ubuntu/cellranger-10.0.0",
    "reference": "/path/to/refdata-gex-GRCh38-2024-A"
})

scanpy_tool = ScanpyTool(cellranger_tool=cellranger_tool)

# 1. 运行 Cell Ranger
result = scanpy_tool.run_cellranger(
    fastq_dir="/path/to/fastqs",
    sample_id="sample_name",
    output_dir="/path/to/output",
    localcores=8,
    localmem=32,
    create_bam=False
)

# 2. 转换为 .h5ad
if result["status"] == "success":
    convert_result = scanpy_tool.convert_cellranger_to_h5ad(
        cellranger_matrix_dir=result["matrix_dir"],
        output_h5ad_path="/path/to/output/sample_name.h5ad"
    )
```

---

## 📋 工作流程

### 从 FASTQ 开始的完整流程：

1. **用户上传 FASTQ 文件或提供目录路径**
2. **Agent 检测文件类型** → 识别为 `fastq`
3. **Agent 调用 `run_cellranger()`**
   - 执行 Cell Ranger count
   - 生成过滤后的表达矩阵
4. **Agent 调用 `convert_cellranger_to_h5ad()`**
   - 将 10x MTX 格式转换为 .h5ad
5. **Agent 调用 `inspect_file()`**
   - 检查转换后的数据
   - 推荐分析参数
6. **用户确认参数**
7. **Agent 执行 Scanpy 分析流程**
   - QC → Normalization → HVG → Scale → PCA → Neighbors → Clustering → UMAP → t-SNE → Markers

### 从 .h5ad 开始的流程：

1. **用户上传 .h5ad 文件**
2. **Agent 检测文件类型** → 识别为 `h5ad`
3. **Agent 调用 `inspect_file()`**
   - 检查数据特征
   - 推荐分析参数
4. **用户确认参数**
5. **Agent 执行 Scanpy 分析流程**

---

## ⚙️ 配置说明

### 环境变量

```bash
# Cell Ranger 安装路径
export CELLRANGER_PATH="/home/ubuntu/cellranger-10.0.0"

# 参考基因组路径
export CELLRANGER_REF="/path/to/refdata-gex-GRCh38-2024-A"
```

### settings.yaml

```yaml
tools:
  cellranger:
    path: "${CELLRANGER_PATH:/home/ubuntu/cellranger-10.0.0}"
    reference: "${CELLRANGER_REF:/path/to/refdata-gex-GRCh38-2024-A}"
    localcores: 8
    localmem: 32
    create_bam: false
```

---

## 🔧 工具方法说明

### CellRangerTool.run_count()

**参数：**
- `fastq_dir`: FASTQ 文件目录路径
- `sample_id`: 样本 ID（也是输出目录名）
- `output_dir`: 最终输出目录路径
- `reference`: 参考基因组路径（可选，使用配置中的默认值）
- `sample`: 样本名称（从 FASTQ 文件名提取，可选）
- `localcores`: CPU 核心数（默认：8）
- `localmem`: 内存 GB（默认：32）
- `create_bam`: 是否创建 BAM 文件（默认：False）
- `expect_cells`: 预期细胞数（可选）

**返回：**
```python
{
    "status": "success" | "error",
    "output_dir": "/path/to/output",
    "matrix_dir": "/path/to/output/filtered_feature_bc_matrix",
    "error": "错误信息（如果有）"
}
```

### ScanpyTool.convert_cellranger_to_h5ad()

**参数：**
- `cellranger_matrix_dir`: Cell Ranger 矩阵目录路径
- `output_h5ad_path`: 输出的 .h5ad 文件路径

**返回：**
```python
{
    "status": "success" | "error",
    "output_path": "/path/to/output.h5ad",
    "n_obs": 1221,
    "n_vars": 38606,
    "matrix_type": "csc_matrix",
    "file_size_mb": 35.2,
    "error": "错误信息（如果有）"
}
```

---

## 📝 注意事项

1. **Cell Ranger 版本要求**
   - 支持 Cell Ranger 10.0.0+
   - `--create-bam` 参数是必需的（设置为 `false` 可节省时间和空间）

2. **参考基因组**
   - 需要提前下载并解压参考基因组
   - 支持 10x Genomics 官方参考基因组

3. **FASTQ 文件结构**
   - 支持标准的 10x Genomics FASTQ 文件结构
   - 文件名应包含样本名称（如 `sample_S1_L001_R1_001.fastq.gz`）

4. **内存和 CPU**
   - 根据服务器配置调整 `localcores` 和 `localmem`
   - 大型数据集需要更多资源

5. **输出目录**
   - Cell Ranger 会在当前工作目录创建输出
   - 建议使用绝对路径指定输出目录

---

## 🎯 下一步优化建议

1. **支持多样本聚合**
   - 实现 `cellranger aggr` 功能
   - 支持批量处理多个样本

2. **自动参数优化**
   - 根据数据规模自动调整 Cell Ranger 参数
   - 智能推荐 `expect_cells` 值

3. **进度监控**
   - 实时显示 Cell Ranger 运行进度
   - 支持长时间运行的任务

4. **错误恢复**
   - 支持从失败点恢复
   - 自动重试机制

---

## 📚 相关文档

- `test_data/CELLRANGER_TUTORIAL.md` - Cell Ranger 使用教程
- `test_data/CELLRANGER_RESULTS.md` - Cell Ranger 运行结果示例
- `单细胞转录组测序全流程提示词.md` - 完整流程提示词

---

**最后更新：** 2024-12-19  
**提交记录：** `8789e97` - feat: 集成 Cell Ranger 到智能体

