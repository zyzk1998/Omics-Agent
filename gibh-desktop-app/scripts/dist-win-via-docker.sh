#!/usr/bin/env bash
# 在 Linux/macOS 上打 Windows NSIS 包：electron-builder 需要 Wine；用官方镜像代替本机安装 Wine。
# 用法（在仓库任意目录）： bash gibh-desktop-app/scripts/dist-win-via-docker.sh
#
# 可选环境变量（宿主机 export 后执行本脚本即可传入容器）：
#   GIBH_ELECTRON_MIRROR — 国内拉 Electron 二进制易断时，可设为：
#     https://npmmirror.com/mirrors/electron/
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ① 旧构建可能以 root 写入 dist-out（如 win-unpacked/resources/elevate.exe），非 root 跑 electron-builder 会 unlinkat: permission denied。
#    同时清理 node_modules，避免属主混杂导致 npm ci EACCES。
echo "[dist-win-via-docker] 清理旧 node_modules 与 dist-out（root 一步）…"
docker run --rm -u 0:0 \
  -v "${REPO_ROOT}:/repo" \
  electronuserland/builder:wine \
  /bin/bash -lc 'rm -rf /repo/gibh-desktop-app/node_modules /repo/gibh-desktop-app/dist-out'

ELECTRON_MIRROR_ARGS=()
if [[ -n "${GIBH_ELECTRON_MIRROR:-}" ]]; then
  ELECTRON_MIRROR_ARGS+=(-e "ELECTRON_MIRROR=${GIBH_ELECTRON_MIRROR}")
  echo "[dist-win-via-docker] 使用 ELECTRON_MIRROR=${GIBH_ELECTRON_MIRROR}"
fi

# ② 以宿主机 UID 写入 /repo；拉长 npm 拉取与 electron postinstall 超时；npm ci 失败时自动重试（缓解 socket hang up）。
echo "[dist-win-via-docker] npm ci + dist:win（宿主机 UID $(id -u):$(id -g)）…"
# 必须 -i：否则 heredoc  stdin 不会进容器，bash -s 立刻 EOF，npm ci 根本不会跑。
exec docker run --rm -i \
  -u "$(id -u):$(id -g)" \
  -v "${REPO_ROOT}:/repo" \
  -w /repo/gibh-desktop-app \
  -e ELECTRON_CACHE=/repo/gibh-desktop-app/.cache/electron \
  -e ELECTRON_BUILDER_CACHE=/repo/gibh-desktop-app/.cache/electron-builder \
  -e HOME=/repo/gibh-desktop-app/.docker-home \
  -e NPM_CONFIG_CACHE=/repo/gibh-desktop-app/.cache/npm \
  -e NPM_CONFIG_FETCH_RETRIES=8 \
  -e NPM_CONFIG_FETCH_RETRY_MINTIMEOUT=20000 \
  -e NPM_CONFIG_FETCH_RETRY_MAXTIMEOUT=120000 \
  -e NPM_CONFIG_FETCH_TIMEOUT=600000 \
  "${ELECTRON_MIRROR_ARGS[@]}" \
  electronuserland/builder:wine \
  /bin/bash -s <<'EOS'
set -euo pipefail
mkdir -p .cache/electron .cache/electron-builder .cache/npm .docker-home

echo "[dist-win-via-docker] 开始 npm ci（拉取 Electron 二进制可能需数分钟，期间终端可能较长时间无新行，属正常现象）…"

npm_ci_with_retry() {
  local max=5
  local n=1
  while true; do
    if npm ci; then
      return 0
    fi
    if [[ "$n" -ge "$max" ]]; then
      echo "[dist-win-via-docker] npm ci 已连续失败 ${max} 次，退出。" >&2
      return 1
    fi
    echo "[dist-win-via-docker] npm ci 失败（${n}/${max}），20s 后重试…（常见原因为下载 Electron 时网络中断 socket hang up）" >&2
    n=$((n + 1))
    sleep 20
  done
}

npm_ci_with_retry
npm run dist:win
EOS
