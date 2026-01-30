#!/bin/bash
# 🧹 项目清理脚本
# 用途：清理临时文件、日志、缓存等
# 安全：只删除明确标识的临时文件，重要文件先归档

set -e  # 遇到错误立即退出

echo "🧹 开始项目清理..."
echo ""

# 创建归档目录
echo "📦 创建归档目录..."
mkdir -p _archive/debug_docs
mkdir -p _archive/fix_docs
mkdir -p _archive/temp_files
echo "✅ 归档目录创建完成"
echo ""

# 1. 删除日志文件
echo "🗑️  删除日志文件..."
rm -f debug.log gibh_agent.log server.log 2>/dev/null || true
rm -f frontend_rna_test_*.log 2>/dev/null || true
rm -rf debug_logs/ 2>/dev/null || true
echo "✅ 日志文件清理完成"
echo ""

# 2. 删除临时测试文件（根目录）
echo "🗑️  删除临时测试文件..."
rm -f test_api_config.py 2>/dev/null || true
rm -f test_execution_preview_separation.py 2>/dev/null || true
rm -f test_tool_retriever.py 2>/dev/null || true
rm -f test_upload.py 2>/dev/null || true
rm -f test_workflow_preparation.py 2>/dev/null || true
rm -f test_workflow_data.csv 2>/dev/null || true
echo "✅ 临时测试文件清理完成"
echo ""

# 3. 删除临时输出文件
echo "🗑️  删除临时输出文件..."
rm -f frontend_rna_test_output.txt 2>/dev/null || true
rm -f frontend_rna_test_summary_*.txt 2>/dev/null || true
rm -f 前端RNA分析测试问题总结.txt 2>/dev/null || true
echo "✅ 临时输出文件清理完成"
echo ""

# 4. 删除示例/遗留文件
echo "🗑️  删除示例/遗留文件..."
rm -f server_celery_example.py 2>/dev/null || true
rm -f app_visualizer.py 2>/dev/null || true
#rm -f whotowork.html 2>/dev/null || true
echo "✅ 示例/遗留文件清理完成"
echo ""

# 5. 删除无扩展名文件（可能是误创建）
echo "🗑️  删除无扩展名文件..."
rm -f "Data flow" 2>/dev/null || true
rm -f "System Architecture" 2>/dev/null || true
echo "✅ 无扩展名文件清理完成"
echo ""

# 6. 清理Python缓存
echo "🗑️  清理Python缓存..."
find . -type d -name "__pycache__" -exec rm -r {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete 2>/dev/null || true
find . -type f -name "*.pyo" -delete 2>/dev/null || true
echo "✅ Python缓存清理完成"
echo ""

# 7. 归档DEBUG文档
echo "📦 归档DEBUG文档..."
mv AI_REPORT_DEBUG_PROMPT.md _archive/debug_docs/ 2>/dev/null || true
mv DEBUG_AI_REPORT_FAILURE.md _archive/debug_docs/ 2>/dev/null || true
mv DEBUG_COMPLETION_SUMMARY.md _archive/debug_docs/ 2>/dev/null || true
mv DEBUG_DIAGNOSIS_FLOW.md _archive/debug_docs/ 2>/dev/null || true
mv DEBUG_PROMPT_AI_REPORT_AND_TOOLS.md _archive/debug_docs/ 2>/dev/null || true
mv DEBUG_REPORT.md _archive/debug_docs/ 2>/dev/null || true
echo "✅ DEBUG文档归档完成"
echo ""

# 8. 归档FIX文档
echo "📦 归档FIX文档..."
mv AI_REPORT_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv FILE_EXPLANATION_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv FINAL_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv FIX_DIAGNOSIS_REPORT.md _archive/fix_docs/ 2>/dev/null || true
mv FIXES_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv FRONTEND_UPLOAD_DEBUG.md _archive/fix_docs/ 2>/dev/null || true
mv LLM_CALL_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv METABOLOMICS_BUG_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv RNA_ANALYSIS_COMPLETE_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv RNA_ANALYSIS_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv RNA_ANALYSIS_FIXES_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv RNA_ANALYSIS_UX_IMPROVEMENTS.md _archive/fix_docs/ 2>/dev/null || true
mv RNA_PARAMETER_EXTRACTION_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv ROUTING_AND_DIAGNOSIS_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv SYNTAX_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv UPLOAD_BUG_FIX.md _archive/fix_docs/ 2>/dev/null || true
mv UPLOAD_FIX_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
mv UX_IMPROVEMENTS_SUMMARY.md _archive/fix_docs/ 2>/dev/null || true
echo "✅ FIX文档归档完成"
echo ""

# 9. 归档其他临时文档
echo "📦 归档其他临时文档..."
mv CHANGELOG_SESSION.md _archive/temp_files/ 2>/dev/null || true
mv FRONTEND_TEST_REPORT.md _archive/temp_files/ 2>/dev/null || true
mv WORKFLOW_PREPARATION_TEST_REPORT.md _archive/temp_files/ 2>/dev/null || true
mv mermaid-diagram-*.png _archive/temp_files/ 2>/dev/null || true
echo "✅ 其他临时文档归档完成"
echo ""

echo "✨ 项目清理完成！"
echo ""
echo "📊 清理统计："
echo "  - 日志文件：已删除"
echo "  - 临时测试文件：已删除"
echo "  - 临时输出文件：已删除"
echo "  - 示例/遗留文件：已删除"
echo "  - Python缓存：已清理"
echo "  - DEBUG文档：已归档到 _archive/debug_docs/"
echo "  - FIX文档：已归档到 _archive/fix_docs/"
echo "  - 临时文档：已归档到 _archive/temp_files/"
echo ""
echo "💡 提示：归档的文件可以在 _archive/ 目录中找到，如需恢复可手动移动回来。"
