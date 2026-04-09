#!/usr/bin/env bash
# 将网页 favicon 同步为 Electron 客户端图标源文件（提交前可执行以保持与 /static/favicon.png 一致）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cp -f "$ROOT/services/nginx/html/static/favicon.png" "$ROOT/gibh-desktop-app/app-icon.png"
echo "已更新 gibh-desktop-app/app-icon.png"
