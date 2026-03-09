#!/bin/bash
# SVG 转 EMF（WPS 编辑兼容性更好）
# 需要安装 Inkscape: sudo apt install inkscape
#
# WPS 导入：插入→图片→来自文件，选择 .svg 或 .emf
# EMF 可右键「取消组合」后编辑文字、颜色、形状
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT="${1:-$PROJECT_ROOT/omics_agent_flowchart.svg}"
OUTPUT="${2:-$PROJECT_ROOT/omics_agent_flowchart.emf}"

if ! command -v inkscape &>/dev/null; then
  echo "请先安装 Inkscape: sudo apt install inkscape"
  exit 1
fi

inkscape "$INPUT" --export-type=emf --export-filename="$OUTPUT" 2>/dev/null || \
inkscape "$INPUT" --export-type=emf --export-file="$OUTPUT" 2>/dev/null || \
inkscape -z -f "$INPUT" -E "$OUTPUT" 2>/dev/null

echo "已生成: $OUTPUT"
echo "可将 EMF 文件直接插入 WPS 演示文稿，支持文字、颜色、布局修改"
