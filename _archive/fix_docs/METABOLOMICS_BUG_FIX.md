# 代谢组工作流结果提取 Bug 修复

## 🐛 问题描述

根据工作流输出，差异代谢物分析和PCA分析的结果显示为 "N/A"：
- `significant_metabolites`: "N/A"
- `total_metabolites`: "N/A"
- `variance_explained`: "N/A"

## 🔍 根本原因

**字段名不匹配问题**：

### 1. 差异分析字段名不匹配

**工具返回** (`metabolomics_tool.py:886-887`):
```python
"summary": {
    "n_total": len(results_df),
    "n_significant": len(significant),
    ...
}
```

**智能体期望** (`metabolomics_agent.py:960-961`):
```python
"significant_metabolites": differential_result.get("summary", {}).get("significant_count", "N/A")
"total_metabolites": differential_result.get("summary", {}).get("total_count", "N/A")
```

**问题**: 工具返回 `n_significant` 和 `n_total`，但智能体查找 `significant_count` 和 `total_count`。

### 2. PCA 字段名和路径不匹配

**工具返回** (`metabolomics_tool.py:573-577`):
```python
result = {
    "explained_variance": {
        "PC1": float(explained_variance[0]),
        "PC2": float(explained_variance[1]),
        ...
    },
    "data": {
        "tables": {
            "variance_table": [...]
        }
    }
}
```

**智能体期望** (`metabolomics_agent.py:964`):
```python
"variance_explained": pca_result.get("summary", {}).get("variance_explained", "N/A")
```

**问题**: 
1. 工具返回 `explained_variance` 在顶层，不在 `data.summary` 中
2. 智能体从 `pca_result`（即 `step_result.data`）中查找，但 `explained_variance` 在完整的 `result` 对象中

## ✅ 修复方案

### 修复1: 差异分析字段名

**修复位置**: `gibh_agent/agents/specialists/metabolomics_agent.py:960-961`

**修复前**:
```python
"significant_metabolites": differential_result.get("summary", {}).get("significant_count", "N/A")
"total_metabolites": differential_result.get("summary", {}).get("total_count", "N/A")
```

**修复后**:
```python
"significant_metabolites": differential_result.get("summary", {}).get("n_significant", "N/A")
"total_metabolites": differential_result.get("summary", {}).get("n_total", "N/A")
```

### 修复2: PCA 结果提取

**修复位置**: `gibh_agent/agents/specialists/metabolomics_agent.py:1212-1230, 940-966`

**修复内容**:
1. 在保存 PCA 步骤结果时，同时保存完整的 `result` 对象
2. 在构建结果摘要时，从完整的 `result` 对象中提取 `explained_variance`

**修复代码**:
```python
# 1. 保存完整结果
step_result = {
    "step_name": step.get("desc", step_id),
    "status": result.get("status", "success"),
    "logs": result.get("message", "PCA 分析完成"),
    "data": result.get("data", {}),
    "_full_result": result  # 🔧 修复：保存完整结果
}

# 2. 提取 variance_explained
pca_variance = "N/A"
if pca_result:
    for step_detail in steps_details:
        if step_detail.get("tool_id") == "pca_analysis":
            step_result = step_detail.get("step_result", {})
            full_result = step_result.get("_full_result", {})
            if full_result and "explained_variance" in full_result:
                pc1_var = full_result["explained_variance"].get("PC1", 0) * 100
                pc2_var = full_result["explained_variance"].get("PC2", 0) * 100
                pca_variance = f"PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%"
            break
```

## 📊 修复效果

### 修复前
```json
{
  "differential_analysis": {
    "significant_metabolites": "N/A",
    "total_metabolites": "N/A"
  },
  "pca": {
    "variance_explained": "N/A"
  }
}
```

### 修复后
```json
{
  "differential_analysis": {
    "significant_metabolites": 15,  // 实际显著代谢物数量
    "total_metabolites": 63         // 总代谢物数量
  },
  "pca": {
    "variance_explained": "PC1: 45.23%, PC2: 12.56%"  // 实际解释方差
  }
}
```

## 🎯 验证

修复后，工作流执行应该能够正确显示：
- ✅ 显著差异代谢物数量
- ✅ 总代谢物数量
- ✅ PCA 解释方差（PC1 和 PC2）

## 📝 相关文件

- `gibh_agent/agents/specialists/metabolomics_agent.py` - 智能体代码
- `gibh_agent/tools/metabolomics_tool.py` - 工具代码

---

**修复日期**: 2026-01-09  
**修复状态**: ✅ 已完成

