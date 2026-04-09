#!/usr/bin/env bash
# 将 electron-builder 产物复制到 Nginx 静态目录，供 /downloads/ 页面分发。
# 用法：在仓库根目录执行  bash scripts/copy_desktop_installer_to_nginx.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/gibh-desktop-app/dist-out"
DST="$ROOT/services/nginx/html/downloads"
mkdir -p "$DST"
shopt -s nullglob
for f in "$SRC"/GIBH-Agent-Demo-Setup-*.{exe,AppImage,blockmap} "$SRC"/*.exe "$SRC"/*.AppImage; do
  [[ -f "$f" ]] || continue
  cp -v "$f" "$DST/"
done
echo "完成。请确认 Nginx 已挂载 $DST 后访问 https://你的域名/downloads/"
