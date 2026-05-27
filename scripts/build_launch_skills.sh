#!/usr/bin/env bash
# 仅重建首发技能隔离栈（不触发 api-server 全量组学层）
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export DOCKER_BUILDKIT=1
LOG="${LOG:-$ROOT/build-launch-skills.log}"
echo "==> 构建 launch-skills（日志: ${LOG}）"
docker compose --progress plain build launch-skills 2>&1 | tee "${LOG}"
echo "==> 启动: docker compose up -d launch-skills"
docker compose up -d launch-skills
docker exec gibh_v2_launch_skills /app/scripts/verify_launch_skill_deps.sh
echo "==> 完成。api-server 将经 LAUNCH_SKILLS_BASE_URL 委托执行首发技能。"
