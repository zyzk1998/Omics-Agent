#!/bin/bash
# Repository Cleanup Script for Public Release
# This script moves temporary debugging files to _archive/ and removes them from git tracking

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ARCHIVE_DIR="$PROJECT_ROOT/_archive"

echo "🧹 Starting repository cleanup..."
echo "📁 Project root: $PROJECT_ROOT"
echo "📦 Archive directory: $ARCHIVE_DIR"
echo ""

# Create archive directory if it doesn't exist
mkdir -p "$ARCHIVE_DIR"
mkdir -p "$ARCHIVE_DIR/scripts"
mkdir -p "$ARCHIVE_DIR/gibh_agent"
mkdir -p "$ARCHIVE_DIR/test_data"

# Function to move and remove from git
move_to_archive() {
    local file="$1"
    local dest_dir="$2"
    
    if [ -f "$PROJECT_ROOT/$file" ]; then
        echo "Moving $file to archive..."
        mkdir -p "$ARCHIVE_DIR/$(dirname "$file")"
        mv "$PROJECT_ROOT/$file" "$ARCHIVE_DIR/$file"
        git rm --cached "$file" 2>/dev/null || true
    elif [ -d "$PROJECT_ROOT/$file" ] && [ ! -d "$ARCHIVE_DIR/$file" ]; then
        echo "Moving directory $file to archive..."
        mv "$PROJECT_ROOT/$file" "$ARCHIVE_DIR/$file"
        git rm -r --cached "$file" 2>/dev/null || true
    fi
}

# Root directory temporary documents
echo "📄 Archiving root directory temporary documents..."
move_to_archive "FEATURE_IMPLEMENTATION_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "PROMPT_OPTIMIZATION_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "SMART_PATH_RESOLUTION_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "DEPLOYMENT_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "ACCESS_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "REFACTORING_SUMMARY.md" "$ARCHIVE_DIR"
move_to_archive "gibh_agent.log" "$ARCHIVE_DIR"
move_to_archive "debug.log" "$ARCHIVE_DIR"

# Subdirectory temporary documents
echo "📄 Archiving subdirectory temporary documents..."
move_to_archive "gibh_agent/IMPLEMENTATION_SUMMARY.md" "$ARCHIVE_DIR/gibh_agent"
move_to_archive "scripts/COMPLETION_SUMMARY.md" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/README_debug_log_collector.md" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_integration_report.md" "$ARCHIVE_DIR/scripts"

# Scripts directory temporary scripts
echo "🐍 Archiving temporary test and debug scripts..."
move_to_archive "scripts/debug_llm_connection.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/debug_llm_real.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/debug_log_collector.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/simulate_frontend.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_ai_report_generation.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_backend_integrity.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_complete_flow.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_frontend_rna_analysis.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_orchestrator_sse.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_rna_ai_report_generation.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_rna_full_workflow.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/test_rna_parameter_extraction.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_full_flow.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_full_stack.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_intelligence.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_llm_flow.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_llm_flow_simple.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_metabolomics_complete.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_report_generation.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_report_quality.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_structure_only.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_ui_logic.py" "$ARCHIVE_DIR/scripts"
move_to_archive "scripts/verify_r1_flow.py" "$ARCHIVE_DIR/scripts"

# Temporary text files
echo "📝 Archiving temporary text files..."
move_to_archive "test_data/新建 Text Document.txt" "$ARCHIVE_DIR/test_data"

echo ""
echo "✅ Cleanup completed!"
echo "📦 Archived files are in: $ARCHIVE_DIR"
echo ""
echo "💡 Next steps:"
echo "   1. Review the changes: git status"
echo "   2. Commit the cleanup: git add .gitignore && git commit -m 'chore: archive temporary debugging files'"
echo "   3. The _archive/ directory is kept locally for reference but excluded from git"
