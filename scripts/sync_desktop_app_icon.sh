#!/usr/bin/env bash
# 以 gibh-desktop-app/app-icon.png 为品牌源图（圆角白底 + 主体，四角透明；若只有灰底源图可先跑
#   python3 scripts/compose_squircle_brand_icon.py <源.png>），同步到：
#   - services/nginx/html/static/favicon.png（全站 favicon）
#   - services/nginx/html/downloads/client-app-icon.png（下载页展示 / 可单独缓存刷新）
# 提交前或替换图标后可执行。
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/gibh-desktop-app/app-icon.png"
cp -f "$SRC" "$ROOT/services/nginx/html/static/favicon.png"
mkdir -p "$ROOT/services/nginx/html/downloads"
cp -f "$SRC" "$ROOT/services/nginx/html/downloads/client-app-icon.png"
echo "已更新 services/nginx/html/static/favicon.png"
echo "已更新 services/nginx/html/downloads/client-app-icon.png（均来自 gibh-desktop-app/app-icon.png）"
