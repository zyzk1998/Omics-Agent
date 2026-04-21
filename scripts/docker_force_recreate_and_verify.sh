#!/usr/bin/env bash
# 生产抢修：down → 无缓存构建 api 服务 → force-recreate → compose ps → 容器内跑 debug_orchestrator_local.py
# 用法（在含 docker-compose.yml 的目录执行，或设置 COMPOSE_FILE）：
#   cd /path/to/stack
#   export COMPOSE_FILE=/path/to/docker-compose.yml   # 若 compose 不在当前目录
#   export API_SERVICE=api-server                      # 若服务名不是 api-server
#   bash /path/to/GIBH-AGENT-V2/scripts/docker_force_recreate_and_verify.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
API_SERVICE="${API_SERVICE:-api-server}"
COMPOSE_FILE="${COMPOSE_FILE:-}"
SCRIPT_HOST="${ROOT}/debug_orchestrator_local.py"

compose() {
  if [[ -n "${COMPOSE_FILE}" ]]; then
    docker compose -f "${COMPOSE_FILE}" "$@"
  else
    docker compose "$@"
  fi
}

if [[ -n "${COMPOSE_FILE}" ]]; then
  if [[ ! -f "${COMPOSE_FILE}" ]]; then
    echo "❌ COMPOSE_FILE 不存在: ${COMPOSE_FILE}"
    exit 2
  fi
  COMPOSE_DIR="$(cd "$(dirname "${COMPOSE_FILE}")" && pwd)"
  cd "${COMPOSE_DIR}"
else
  if [[ ! -f docker-compose.yml ]] && [[ ! -f compose.yml ]]; then
    echo "❌ 当前目录无 docker-compose.yml / compose.yml，且未设置 COMPOSE_FILE。"
    echo "   请在生产栈目录执行，或: export COMPOSE_FILE=/绝对路径/docker-compose.yml"
    exit 2
  fi
fi

echo "==> docker compose down"
compose down

echo "==> docker compose build --no-cache ${API_SERVICE}"
compose build --no-cache "${API_SERVICE}"

echo "==> docker compose up -d --force-recreate"
compose up -d --force-recreate

echo ""
echo "======== docker compose ps ========"
compose ps

echo ""
echo "======== docker ps (节选) ========"
docker ps --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" | head -35

if [[ ! -f "${SCRIPT_HOST}" ]]; then
  echo "⚠️  未找到 ${SCRIPT_HOST}，跳过容器内验证。"
  exit 0
fi

echo ""
echo "==> 将 debug 脚本拷入容器并在 /tmp 执行（工作目录按常见 /app）"
cid="$(compose ps -q "${API_SERVICE}" | head -1)"
if [[ -z "${cid}" ]]; then
  echo "❌ 未找到运行中的容器: ${API_SERVICE}"
  exit 3
fi

APP_ROOT="${APP_ROOT:-/app}"
docker cp "${SCRIPT_HOST}" "${cid}:${APP_ROOT}/debug_orchestrator_local.py"

echo "==> docker compose exec ${API_SERVICE} python ${APP_ROOT}/debug_orchestrator_local.py"
compose exec -T "${API_SERVICE}" env GIBH_AGENT_ROOT="${APP_ROOT}" python "${APP_ROOT}/debug_orchestrator_local.py"
