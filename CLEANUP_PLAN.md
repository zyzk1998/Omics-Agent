# 🧹 项目清理计划

## 清理策略
- **安全第一**：不确定的文件先归档到 `_archive/` 目录
- **保留核心**：保留所有生产代码和正式测试
- **清理临时**：删除日志、临时测试文件、示例文件

---

## 📋 清理清单

### ✅ 确定删除的文件（根目录）

#### 1. 日志文件
- `debug.log`
- `gibh_agent.log`
- `server.log`
- `frontend_rna_test_20260128_095800.log`
- `frontend_rna_test_20260128_100330.log`
- `debug_logs/` (整个目录)

#### 2. 临时测试文件（根目录）
- `test_api_config.py`
- `test_execution_preview_separation.py`
- `test_tool_retriever.py`
- `test_upload.py`
- `test_workflow_preparation.py`
- `test_workflow_data.csv`

#### 3. 临时输出文件
- `frontend_rna_test_output.txt`
- `frontend_rna_test_summary_20260128_095800.txt`
- `frontend_rna_test_summary_20260128_100332.txt`
- `前端RNA分析测试问题总结.txt`

#### 4. 示例/遗留文件
- `server_celery_example.py` (明确标注为示例)
- `app_visualizer.py` (V1 Streamlit遗留)
- `whotowork.html` (临时文档)

#### 5. 无扩展名文件（可能是误创建）
- `Data flow`
- `System Architecture`

#### 6. Python缓存
- 所有 `__pycache__/` 目录（会被.gitignore忽略，但可以清理）

---

### 📦 归档到 `_archive/` 的文档

以下文档是开发过程中的调试和修复总结，可以归档：

#### DEBUG文档
- `AI_REPORT_DEBUG_PROMPT.md`
- `DEBUG_AI_REPORT_FAILURE.md`
- `DEBUG_COMPLETION_SUMMARY.md`
- `DEBUG_DIAGNOSIS_FLOW.md`
- `DEBUG_PROMPT_AI_REPORT_AND_TOOLS.md`
- `DEBUG_REPORT.md`

#### FIX总结文档
- `AI_REPORT_FIX_SUMMARY.md`
- `FILE_EXPLANATION_FIX.md`
- `FINAL_FIX_SUMMARY.md`
- `FIX_DIAGNOSIS_REPORT.md`
- `FIXES_SUMMARY.md`
- `FRONTEND_UPLOAD_DEBUG.md`
- `LLM_CALL_FIX.md`
- `METABOLOMICS_BUG_FIX.md`
- `RNA_ANALYSIS_COMPLETE_FIX.md`
- `RNA_ANALYSIS_FIX_SUMMARY.md`
- `RNA_ANALYSIS_FIXES_SUMMARY.md`
- `RNA_ANALYSIS_UX_IMPROVEMENTS.md`
- `RNA_PARAMETER_EXTRACTION_FIX.md`
- `ROUTING_AND_DIAGNOSIS_FIX_SUMMARY.md`
- `SYNTAX_FIX_SUMMARY.md`
- `UPLOAD_BUG_FIX.md`
- `UPLOAD_FIX_SUMMARY.md`
- `UX_IMPROVEMENTS_SUMMARY.md`

#### 其他临时文档
- `CHANGELOG_SESSION.md`
- `FRONTEND_TEST_REPORT.md`
- `WORKFLOW_PREPARATION_TEST_REPORT.md`
- `mermaid-diagram-2026-01-09-145507.png`

---

### 🔒 保留的文件

#### 核心文档
- `README.md`
- `API.md`
- `SETUP.md`
- `ARCHITECTURE_REFACTORING.md`
- `DEPLOYMENT_SUMMARY.md`
- `PROJECT_REPORT.md`
- `QUICK_REFERENCE.md`
- `SYSTEM_AUDIT.md`

#### 正式测试
- `scripts/` 目录下的所有测试文件
- `tests/` 目录下的所有测试文件

#### 生产代码
- `server.py`
- `gibh_agent/` 目录
- `services/` 目录
- `docker-compose.yml`
- `requirements.txt`

---

## 🚀 执行命令

```bash
# 1. 创建归档目录
mkdir -p _archive/debug_docs _archive/fix_docs _archive/temp_files

# 2. 删除日志文件
rm -f *.log
rm -rf debug_logs/

# 3. 删除临时测试文件
rm -f test_*.py test_*.csv

# 4. 删除临时输出文件
rm -f frontend_rna_test_*.txt frontend_rna_test_output.txt 前端RNA分析测试问题总结.txt

# 5. 删除示例/遗留文件
rm -f server_celery_example.py app_visualizer.py whotowork.html

# 6. 删除无扩展名文件
rm -f "Data flow" "System Architecture"

# 7. 清理Python缓存
find . -type d -name "__pycache__" -exec rm -r {} + 2>/dev/null || true

# 8. 归档DEBUG文档
mv AI_REPORT_DEBUG_PROMPT.md DEBUG_*.md _archive/debug_docs/ 2>/dev/null || true

# 9. 归档FIX文档
mv *_FIX*.md *_FIXES*.md _archive/fix_docs/ 2>/dev/null || true

# 10. 归档其他临时文档
mv CHANGELOG_SESSION.md FRONTEND_TEST_REPORT.md WORKFLOW_PREPARATION_TEST_REPORT.md mermaid-diagram-*.png _archive/temp_files/ 2>/dev/null || true
```

---

## ⚠️ 注意事项

1. **备份**：执行前建议先备份整个项目
2. **Git状态**：确保所有重要更改已提交
3. **确认**：执行前请确认清单无误

---

## 📊 清理统计

- **删除文件数**：约 25+ 个文件
- **归档文件数**：约 20+ 个文档
- **清理缓存**：所有 `__pycache__/` 目录
- **预计节省空间**：取决于日志文件大小
