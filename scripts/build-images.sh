#!/bin/bash
# 一键生成所有图片到 output-images/
set -e
cd "$(dirname "$0")/.."

echo "生成系统拓扑图..."
node scripts/svg-to-png.js static_topology.svg output-images/static_topology.png

echo "生成市场份额饼图..."
node scripts/generate-pie-chart.js

echo "复制 Mermaid 流程图..."
cp mermaid-workflows/*.png output-images/ 2>/dev/null || true

echo "完成。图片位于 output-images/"
ls -la output-images/*.png
