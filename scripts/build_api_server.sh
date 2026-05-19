#!/usr/bin/env bash
# 构建 api-server 镜像（实时输出，避免 `| tail -N` 导致长时间无回显）
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

INSTALL_BEPIPRED3="${INSTALL_BEPIPRED3:-1}"
LOG="${LOG:-$ROOT/build-api-server.log}"

export DOCKER_BUILDKIT=1
echo "==> 构建 api-server（INSTALL_BEPIPRED3=${INSTALL_BEPIPRED3}）"
echo "==> 日志: ${LOG}"
echo "==> 勿使用 '| tail -N'：会等构建结束才显示输出"

docker compose --progress plain build \
  --build-arg "INSTALL_BEPIPRED3=${INSTALL_BEPIPRED3}" \
  api-server 2>&1 | tee "${LOG}"

echo "==> 完成。启动: docker compose up -d api-server"
