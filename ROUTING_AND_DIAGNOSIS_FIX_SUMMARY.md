# 路由和诊断修复总结

## 修复日期
2025-01-28

## 问题分析

### 问题1：路由错误
**现象**：用户选择FASTQ文件，需求是RNA全流程分析，但系统路由到了代谢组分析工作流。

**根本原因**：
1. 在 `planner.py` 的 `_classify_intent` 方法中，文件类型检查逻辑只检查文件扩展名
2. 对于FASTQ目录（如 `/app/test_data/pbmc_1k_v3_fastqs`），这是一个目录，不是文件，所以 `file_ext` 会是空字符串或目录名
3. 没有检查 `file_type` 字段，导致FASTQ目录没有被正确识别为RNA类型

### 问题2：FASTQ文件诊断
**现象**：FASTQ文件过大，不适合直接进行数据诊断。

**根本原因**：
1. 没有专门的FASTQ目录检查器
2. 没有诊断结果缓存机制
3. 每次规划都需要重新诊断，对于大型文件效率低下

## 修复方案

### 1. 创建FASTQ目录检查器

**文件**：`gibh_agent/core/file_inspector.py`

**新增类**：`FastqDirectoryHandler`
- 优先级：8（高于其他目录检查器）
- 功能：
  - 检测包含FASTQ文件的目录（.fastq, .fq, .fastq.gz, .fq.gz）
  - 提供轻量级诊断：文件数量、总大小、文件类型（R1, R2, I1等）
  - 识别是否包含配对端数据和索引文件
  - 判断是否为10x格式

**关键代码**：
```python
@register_inspector(priority=8)
class FastqDirectoryHandler(BaseFileHandler):
    """FASTQ 目录检查器"""
    
    def can_handle(self, path: Path) -> bool:
        """检查是否为 FASTQ 目录"""
        if not path.is_dir():
            return False
        # 检查目录中是否包含 FASTQ 文件
        fastq_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
        # ... 递归搜索逻辑
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """检查 FASTQ 目录，返回轻量级诊断信息"""
        # 收集所有FASTQ文件
        # 计算总大小
        # 识别文件类型（R1, R2, I1等）
        # 返回诊断信息
```

### 2. 修复路由逻辑

**文件**：`gibh_agent/core/planner.py`

**修改**：在 `_classify_intent` 方法中，增强文件类型检查逻辑

**关键修改**：
```python
# 🔥 TASK 2 FIX: Check file_type from metadata (more reliable than extension for directories)
if file_type == "fastq":
    if domain_name == "Metabolomics":
        logger.warning(f"⚠️ LLM 将 FASTQ 目录分类为 Metabolomics，强制覆盖为 RNA")
        domain_name = "RNA"
    else:
        logger.info(f"✅ FASTQ 目录已正确分类为 {domain_name}")
elif file_type == "10x_mtx" or file_type == "anndata":
    if domain_name == "Metabolomics":
        logger.warning(f"⚠️ LLM 将 {file_type} 文件分类为 Metabolomics，强制覆盖为 RNA")
        domain_name = "RNA"
```

**改进点**：
- 不仅检查文件扩展名，还检查 `file_type` 字段
- 对于目录类型的文件（FASTQ目录、10x目录），使用 `file_type` 更可靠
- 添加了对 `10x_mtx` 和 `anndata` 类型的检查

### 3. 实现诊断结果缓存机制

**新增文件**：`gibh_agent/core/diagnosis_cache.py`

**类**：`DiagnosisCache`
- 功能：
  - 将诊断结果保存到原始数据文件同文件夹下
  - 缓存文件命名：`{原文件名}.diagnosis.json`
  - 支持缓存失效检查（基于文件哈希）
  - 在规划阶段检查并加载已保存的诊断结果

**关键方法**：
- `save_diagnosis(file_path, diagnosis_result, file_metadata)`: 保存诊断结果
- `load_diagnosis(file_path)`: 加载诊断结果
- `has_cache(file_path)`: 检查是否存在缓存
- `clear_cache(file_path)`: 清除缓存

**缓存文件格式**：
```json
{
  "file_path": "/app/test_data/pbmc_1k_v3_fastqs",
  "file_hash": "abc123...",
  "diagnosis_result": {
    "diagnosis_report": "...",
    "stats": {...},
    "recommendation": {...},
    "omics_type": "scRNA"
  },
  "file_metadata": {...},
  "cached_at": "1234567890.0"
}
```

### 4. 集成诊断缓存到BaseAgent

**文件**：`gibh_agent/agents/base_agent.py`

**修改**：`_perform_data_diagnosis` 方法

**关键修改**：
```python
# 🔥 TASK 4: 检查诊断缓存
from ..core.diagnosis_cache import get_diagnosis_cache
cache = get_diagnosis_cache()

file_path = file_metadata.get("file_path", "")
if file_path:
    cached_diagnosis = cache.load_diagnosis(file_path)
    if cached_diagnosis:
        logger.info(f"✅ [DataDiagnostician] 使用缓存的诊断结果: {file_path}")
        return cached_diagnosis.get("diagnosis_report")

# ... 执行诊断 ...

# 🔥 TASK 4: 保存诊断结果到缓存
if file_path and response:
    cache_data = {
        "diagnosis_report": response,
        "stats": stats,
        "recommendation": recommendation,
        "omics_type": omics_type
    }
    cache.save_diagnosis(file_path, cache_data, file_metadata)
```

**工作流程**：
1. 检查是否有已保存的诊断结果
2. 如果有，直接返回缓存的诊断报告
3. 如果没有，执行诊断
4. 诊断完成后，保存结果到缓存

## 修复效果

### 路由修复
- ✅ FASTQ目录现在会被正确识别为RNA类型
- ✅ 10x目录和h5ad文件也会被正确识别为RNA类型
- ✅ 文件类型检查逻辑更加健壮，不仅检查扩展名，还检查 `file_type` 字段

### 诊断缓存
- ✅ FASTQ文件等大型文件的诊断结果会被缓存
- ✅ 规划阶段如果检测到已有诊断结果，直接使用，提高效率
- ✅ 缓存文件保存在原始数据文件同文件夹下，方便管理
- ✅ 支持缓存失效检查，确保数据修改后重新诊断

## 测试建议

1. **测试路由修复**：
   - 上传FASTQ目录，请求RNA分析
   - 确认系统正确路由到RNA工作流
   - 上传CSV文件，请求代谢组分析
   - 确认系统正确路由到代谢组工作流

2. **测试诊断缓存**：
   - 首次上传FASTQ文件，执行诊断
   - 检查是否在文件同目录下生成了 `.diagnosis.json` 文件
   - 再次上传相同文件，确认使用缓存的诊断结果
   - 修改文件后，确认重新诊断

3. **测试FASTQ目录检查器**：
   - 上传包含FASTQ文件的目录
   - 确认系统正确识别为FASTQ类型
   - 确认诊断信息包含文件数量、总大小、文件类型等

## 注意事项

1. **缓存文件位置**：
   - 缓存文件保存在原始数据文件同文件夹下
   - 格式：`{原文件名}.diagnosis.json`
   - 如果原始文件是目录，缓存文件也在该目录下

2. **缓存失效**：
   - 基于文件哈希（文件大小和修改时间）
   - 如果文件被修改，缓存会自动失效
   - 可以手动删除缓存文件强制重新诊断

3. **文件类型检测优先级**：
   - FASTQ目录检查器优先级为8（最高）
   - 10x目录检查器优先级为10
   - AnnData检查器优先级为5
   - Tabular检查器优先级为1

---

**修复完成日期**: 2025-01-28  
**修复人员**: AI Assistant  
**状态**: ✅ 已完成
