#!/usr/bin/env bash
# 在无 Windows 主机时，用 Wine + Windows 嵌入式 Python 在本机构建 local_sidecar/dist/local_sidecar.exe
# 依赖：wine64 / wine、curl、unzip；适用于 electronuserland/builder:wine 镜像内执行。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SIDECAR="${APP_ROOT}/local_sidecar"
DIST="${SIDECAR}/dist"
CACHE="${SIDECAR}/.cache-win-embed-py"
PY_VER="${GIBH_WIN_EMBED_PY_VER:-3.11.9}"
ZIP_NAME="python-${PY_VER}-embed-amd64.zip"
ZIP_URL="https://www.python.org/ftp/python/${PY_VER}/${ZIP_NAME}"

command -v wine >/dev/null 2>&1 || { echo "[build-sidecar-wine] 未找到 wine，无法交叉生成 .exe" >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { echo "[build-sidecar-wine] 未找到 curl" >&2; exit 1; }
command -v unzip >/dev/null 2>&1 || { echo "[build-sidecar-wine] 未找到 unzip" >&2; exit 1; }

export WINEARCH="${WINEARCH:-win64}"
export WINEPREFIX="${WINEPREFIX:-${SIDECAR}/.wine-sidecar-build}"

mkdir -p "${CACHE}" "${DIST}"
cd "${CACHE}"

if [[ ! -f "${ZIP_NAME}" ]]; then
  echo "[build-sidecar-wine] 下载嵌入式 Python ${PY_VER} …"
  curl -fsSL -o "${ZIP_NAME}" "${ZIP_URL}"
fi

rm -rf "${CACHE}/embed"
mkdir -p "${CACHE}/embed"
unzip -q -o "${ZIP_NAME}" -d "${CACHE}/embed"

# 启用 site-packages（嵌入式发行版默认注释掉 import site）
PY_DOT="$(echo "${PY_VER}" | cut -d. -f1-2 | tr -d .)"
PTH_FILE=""
if [[ -f "${CACHE}/embed/python${PY_DOT}._pth" ]]; then
  PTH_FILE="${CACHE}/embed/python${PY_DOT}._pth"
elif [[ -f "${CACHE}/embed/python311._pth" ]]; then
  PTH_FILE="${CACHE}/embed/python311._pth"
else
  echo "[build-sidecar-wine] 未找到 python*._pth" >&2
  ls -la "${CACHE}/embed" >&2
  exit 1
fi
if grep -q '^#import site' "${PTH_FILE}"; then
  sed -i 's/^#import site/import site/' "${PTH_FILE}"
elif ! grep -q '^import site' "${PTH_FILE}"; then
  echo 'import site' >> "${PTH_FILE}"
fi

GET_PIP="${CACHE}/get-pip.py"
if [[ ! -f "${GET_PIP}" ]]; then
  curl -fsSL -o "${GET_PIP}" https://bootstrap.pypa.io/get-pip.py
fi

echo "[build-sidecar-wine] Wine 初始化 / 安装 pip（首次可能较慢）…"
wineboot -i 2>/dev/null || true
(cd "${CACHE}/embed" && wine python.exe "${GET_PIP}" --no-warn-script-location)

echo "[build-sidecar-wine] 安装依赖与 PyInstaller …"
(cd "${CACHE}/embed" && wine python.exe -m pip install -q -U pip)
(cd "${CACHE}/embed" && wine python.exe -m pip install -q "pyinstaller>=6.0" -r "${SIDECAR}/requirements.txt")

echo "[build-sidecar-wine] PyInstaller 打包 …"
rm -rf "${CACHE}/embed/build" "${CACHE}/embed/dist" "${SIDECAR}/local_sidecar-wine.spec"
(cd "${CACHE}/embed" && wine python.exe -m PyInstaller \
  --onefile \
  --clean \
  --name local_sidecar \
  --workpath "${CACHE}/embed/build" \
  --distpath "${CACHE}/embed/dist" \
  --specpath "${SIDECAR}" \
  --collect-all uvicorn \
  --collect-all fastapi \
  --collect-all pydantic \
  --hidden-import uvicorn.lifespan \
  --hidden-import uvicorn.lifespan.on \
  --hidden-import uvicorn.protocols.http.auto \
  --hidden-import uvicorn.protocols.websockets.auto \
  "${SIDECAR}/local_server.py")

EXE_OUT="${CACHE}/embed/dist/local_sidecar.exe"
if [[ ! -f "${EXE_OUT}" ]]; then
  echo "[build-sidecar-wine] 未生成 ${EXE_OUT}" >&2
  exit 1
fi

cp -f "${EXE_OUT}" "${DIST}/local_sidecar.exe"
echo "[build-sidecar-wine] 完成: ${DIST}/local_sidecar.exe"
