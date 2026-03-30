#!/usr/bin/env bash
# 在 third_party/BepiPred-3.0 下创建 .venv 并安装官方依赖 + PyTorch（供本地运行 API、无 Docker 时使用）。
# Docker 镜像请使用 services/api/Dockerfile 内的 /opt/bepipred3-venv。
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BP="${BEPIPRED3_ROOT:-$ROOT/third_party/BepiPred-3.0}"
TORCH_INDEX="${BEPIPRED_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"

if [[ ! -f "$BP/requirements.txt" ]]; then
  echo "未找到 $BP/requirements.txt，请先克隆 BepiPred-3.0 到 third_party/BepiPred-3.0 或设置 BEPIPRED3_ROOT" >&2
  exit 1
fi

PY="${PYTHON_BIN:-python3}"
"$PY" -m venv "$BP/.venv"
# shellcheck disable=SC1090
source "$BP/.venv/bin/activate"
pip install --upgrade pip
pip install torch torchvision --index-url "$TORCH_INDEX"
pip install -r "$BP/requirements.txt"
echo "完成。请在运行 API 前设置: export BEPIPRED3_PYTHON=$BP/.venv/bin/python"
