# 数据诊断功能重构总结

## 📋 重构目标
将数据诊断功能从各个 Agent 中提取到统一的核心模块，使所有 Agent 自动获得数据诊断能力。

## ✅ 完成的工作

### Step 1: 创建 `gibh_agent/core/data_diagnostician.py`
- **DataDiagnostician 类**：统一的数据诊断器
- **支持的组学类型**：
  - `scRNA`：单细胞转录组（计算 n_cells, n_genes, QC 指标）
  - `Metabolomics`：代谢组学（计算 missing_rate, zero_rate, data_range）
  - `BulkRNA`：Bulk RNA-seq
  - `default`：通用统计
- **功能**：
  - 基于文件元数据计算统计事实
  - 数据质量评估
  - 参数推荐逻辑（基于数据特征）

### Step 2: 更新 `gibh_agent/agents/base_agent.py`
- **新增方法**：`_perform_data_diagnosis()`
  - 统一的诊断入口
  - 调用 `DataDiagnostician.analyze()` 计算统计
  - 使用 LLM 生成 Markdown 报告
  - 将诊断报告保存到 `self.context["diagnosis_report"]`
- **初始化**：添加 `self.diagnostician` 和 `self.context`

### Step 3: 重构 `gibh_agent/agents/specialists/rna_agent.py`
- **移除**：`_generate_diagnosis_and_recommendation()` 方法
- **更新**：`_generate_workflow_config()` 使用 `BaseAgent._perform_data_diagnosis()`
- **传递参数**：`omics_type="scRNA"`

### Step 4: 更新 `gibh_agent/agents/specialists/metabolomics_agent.py`
- **移除**：`_generate_diagnosis_and_recommendation()` 方法
- **更新**：`_generate_workflow_config()` 使用 `BaseAgent._perform_data_diagnosis()`
- **传递参数**：`omics_type="Metabolomics"`
- **修复**：确保诊断报告被包含在返回结果中（`result["diagnosis_report"]`）

## 🔄 工作流程

```
用户上传文件
    ↓
FileInspector.inspect_file() → 文件元数据
    ↓
BaseAgent._perform_data_diagnosis()
    ├─ DataDiagnostician.analyze() → 统计事实
    └─ LLM 生成 Markdown 报告
    ↓
保存到 self.context["diagnosis_report"]
    ↓
返回给前端（在 workflow_config 结果中）
```

## 📊 诊断报告格式

诊断报告是 Markdown 格式，包含：
- **数据体检报告**：数据规模、特征、质量
- **参数推荐表**：参数名、默认值、推荐值、推荐理由
- **下一步建议**

## 🎯 优势

1. **统一性**：所有 Agent 使用相同的诊断逻辑
2. **可扩展性**：新增 Agent 自动获得诊断能力
3. **可维护性**：诊断逻辑集中在一个模块
4. **数据驱动**：基于统计事实生成推荐，减少幻觉

## 🔍 验证

- ✅ 语法检查通过
- ✅ 导入测试通过
- ✅ Linter 检查通过

## 📝 注意事项

- `MetabolomicsAgent` 仍然保留 `_generate_parameter_recommendations()` 方法（用于参数提取）
- 诊断报告会自动包含在 `workflow_config` 返回结果中
- 如果诊断失败，不会阻塞工作流生成（返回 None）
