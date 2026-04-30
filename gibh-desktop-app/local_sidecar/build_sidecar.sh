#!/usr/bin/env bash
# 在 **当前操作系统** 上将 local_server.py 打成单文件可执行体，供 Electron extraResources 使用。
# 产物：dist/local_sidecar（Linux/macOS）或 dist/local_sidecar.exe（Windows）
# 发版前：在目标平台各执行一次本脚本，再执行 npm run dist / dist:win / dist:linux。
set -euo pipefail
ROOT="$(cd "$(dirname "${0}")" && pwd)"
cd "${ROOT}"

PY="${OMICS_SIDECAR_PYTHON:-}"
if [[ -z "${PY}" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PY=python3
  elif command -v python >/dev/null 2>&1; then
    PY=python
  else
    echo "未找到 python3/python，请安装 Python 3 或设置 OMICS_SIDECAR_PYTHON" >&2
    exit 1
  fi
fi

echo "[build_sidecar] 使用解释器: ${PY}"
${PY} -m pip install -q -U pip
${PY} -m pip install -q -r "${ROOT}/requirements.txt" "pyinstaller>=6.0"

rm -rf "${ROOT}/build" "${ROOT}/dist"
# --collect-all 将依赖数据文件一并打入，避免 uvicorn/fastapi 运行期缺模块
${PY} -m PyInstaller \
  --onefile \
  --clean \
  --name local_sidecar \
  --collect-all uvicorn \
  --collect-all fastapi \
  --collect-all pydantic \
  --hidden-import uvicorn.lifespan \
  --hidden-import uvicorn.lifespan.on \
  --hidden-import uvicorn.protocols.http.auto \
  --hidden-import uvicorn.protocols.websockets.auto \
  "${ROOT}/local_server.py"

if [[ -f "${ROOT}/dist/local_sidecar.exe" ]]; then
  echo "[build_sidecar] 完成: ${ROOT}/dist/local_sidecar.exe"
else
  echo "[build_sidecar] 完成: ${ROOT}/dist/local_sidecar"
fi
