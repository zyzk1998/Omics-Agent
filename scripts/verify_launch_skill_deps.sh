#!/usr/bin/env bash
# 验收首发 10 项 BaseSkill 在 api-server 容器（或同构 venv）内的运行时依赖。
# 用法（推荐 launch-skills 隔离栈）：
#   docker exec gibh_v2_launch_skills /app/scripts/verify_launch_skill_deps.sh
#   ./scripts/build_launch_skills.sh   # 构建并自动执行本脚本
# 宿主机：PYTHON=/usr/local/bin/python3 ./scripts/verify_launch_skill_deps.sh
set -euo pipefail

# api-server 主进程为 /usr/local/bin/python3；PATH 前置 omics-conda 时裸 python3 会误用 3.11 且无 rdkit/httpx
PYTHON="${PYTHON:-/usr/local/bin/python3}"
if [[ ! -x "$PYTHON" ]]; then
  PYTHON="$(command -v python3 || echo python3)"
fi
fail=0

check_cli() {
  local name="$1"
  if command -v "$name" >/dev/null 2>&1; then
    echo "OK  CLI $name -> $(command -v "$name")"
  else
    echo "FAIL CLI $name 未在 PATH 中（见 services/api/Dockerfile apt: ncbi-blast+）"
    fail=1
  fi
}

check_py() {
  local snippet="$1"
  local label="$2"
  if "$PYTHON" -c "$snippet" >/dev/null 2>&1; then
    echo "OK  pip $label"
  else
    echo "FAIL pip $label — $snippet"
    fail=1
  fi
}

echo "==> verify_launch_skill_deps (PYTHON=$PYTHON)"
check_cli blastn
check_cli blastp
check_py "import httpx" httpx
check_py "import chembl_webresource_client" chembl_webresource_client
check_py "import rdkit; from rdkit.Chem import AllChem" rdkit-pypi

if [[ "$fail" -ne 0 ]]; then
  echo "==> 依赖验收失败；请重建 api-server: ./scripts/build_api_server.sh && docker compose up -d api-server"
  exit 1
fi
echo "==> 首发技能运行时依赖验收通过"
