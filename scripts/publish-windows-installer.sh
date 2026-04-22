#!/usr/bin/env bash
# 将本机已打好的 Windows NSIS 安装包发布到静态下载目录并刷新 client-release.json。
# 无需在服务器上「拉仓库再打包」，只需把 .exe 传到仓库所在机器（或直接在仓库根执行）。
#
# 用法：
#   ./scripts/publish-windows-installer.sh /path/to/Omics\ Agent\ Setup\ 0.12.0.exe
#   ./scripts/publish-windows-installer.sh ~/Downloads/setup.exe   # 任意文件名，会按 package.json 重命名为规范名
#
# 依赖：仓库根目录执行；需 node（用于 npm run sync:download-page）。
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${1:-}"
MIN_BYTES=$((5 * 1024 * 1024))

if [[ -z "$SRC" ]]; then
  echo "用法: $0 <安装包.exe 的路径>" >&2
  echo "示例（路径含空格务必加引号）:" >&2
  echo "  $0 \"\$HOME/Downloads/Omics Agent Setup 0.12.0.exe\"" >&2
  echo "  $0 ./dist-out/Omics\\\\ Agent\\\\ Setup\\\\ 0.12.0.exe" >&2
  exit 1
fi
if [[ ! -f "$SRC" ]]; then
  echo "错误：找不到文件（路径错误或尚未上传）：" >&2
  echo "  $SRC" >&2
  exit 1
fi

EXPECTED_NAME="$(node -e "
const p = require('$ROOT/gibh-desktop-app/package.json');
const pn = (p.build && p.build.productName) || p.name || 'Omics Agent';
const v = p.version || '0.0.0';
console.log(pn + ' Setup ' + v + '.exe');
")"

SIZE=$(stat -c%s "$SRC" 2>/dev/null || stat -f%z "$SRC" 2>/dev/null || wc -c <"$SRC" | tr -d ' ')
if [[ "${SIZE:-0}" -lt "$MIN_BYTES" ]]; then
  echo "拒绝：文件过小 (${SIZE} bytes)，完整 Electron 安装包通常为数十 MB。" >&2
  exit 1
fi

DOWNLOADS="$ROOT/services/nginx/html/downloads"
mkdir -p "$DOWNLOADS"
DEST="$DOWNLOADS/$EXPECTED_NAME"
cp -f "$SRC" "$DEST"
echo "已复制 -> $DEST (${SIZE} bytes)"

cd "$ROOT/gibh-desktop-app"
export GIBH_CLIENT_UNPIN_WIN=1
npm run sync:download-page

echo "完成。下载页 meta: $DOWNLOADS/client-release.json"
