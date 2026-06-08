#!/usr/bin/env bash
# 开发环境 Label Studio 启动（Docker 镜像拉取失败时的 pip 兜底）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p ls_data

export LABEL_STUDIO_USERNAME="${LABEL_STUDIO_USERNAME:-system_admin@omics.local}"
export LABEL_STUDIO_PASSWORD="${LABEL_STUDIO_PASSWORD:-AutoBootstrap_O1m2i3c4s!}"
export LABEL_STUDIO_DISABLE_SIGNUP_WITHOUT_LINK="${LABEL_STUDIO_DISABLE_SIGNUP_WITHOUT_LINK:-true}"

if ! command -v label-studio >/dev/null 2>&1; then
  echo "安装 label-studio…"
  pip install 'label-studio>=1.12.0,<1.15' -q
fi

if python3 -c "import requests; requests.get('http://127.0.0.1:8082/user/login/',timeout=2).raise_for_status()" 2>/dev/null; then
  echo "Label Studio 已在 http://127.0.0.1:8082 运行"
  exit 0
fi

nohup label-studio start --port 8082 --data-dir "$ROOT/ls_data" > /tmp/label-studio.log 2>&1 &
echo "Label Studio 启动中… 日志: /tmp/label-studio.log"
for i in $(seq 1 30); do
  if python3 -c "import requests; requests.get('http://127.0.0.1:8082/user/login/',timeout=2).raise_for_status()" 2>/dev/null; then
    echo "OK: http://127.0.0.1:8082"
    exit 0
  fi
  sleep 2
done
echo "超时：请检查 /tmp/label-studio.log" >&2
exit 1
