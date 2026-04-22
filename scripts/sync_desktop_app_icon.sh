#!/usr/bin/env bash
# 将 **仅用于瘦客户端** 的 gibh-desktop-app/app-icon.png 同步到：
#   - services/nginx/html/downloads/client-app-icon.png（与安装包/下载说明一致；与 Web UI 无强绑）
# 全站侧栏/收藏夹/标签页用图请单独维护：services/nginx/html/static/favicon.png
#（莫比乌斯主品牌，勿用本脚本从 app-icon 覆盖，避免与客户端 squircle 搞混。）
# 若需同时更新 client-release 元数据，请用：cd gibh-desktop-app && npm run sync:download-page
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/gibh-desktop-app/app-icon.png"
mkdir -p "$ROOT/services/nginx/html/downloads"
cp -f "$SRC" "$ROOT/services/nginx/html/downloads/client-app-icon.png"
echo "已更新 services/nginx/html/downloads/client-app-icon.png（仅客户端源图，未改 favicon）"
