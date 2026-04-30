@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

if defined OMICS_SIDECAR_PYTHON (
  set "PY=%OMICS_SIDECAR_PYTHON%"
) else (
  where python3 >nul 2>&1 && set "PY=python3" || set "PY=python"
)

echo [build_sidecar] 使用解释器: %PY%
%PY% -m pip install -q -U pip
%PY% -m pip install -q -r "%~dp0requirements.txt" "pyinstaller>=6.0"

if exist "%~dp0build" rmdir /s /q "%~dp0build"
if exist "%~dp0dist" rmdir /s /q "%~dp0dist"

%PY% -m PyInstaller ^
  --onefile ^
  --clean ^
  --name local_sidecar ^
  --collect-all uvicorn ^
  --collect-all fastapi ^
  --collect-all pydantic ^
  --hidden-import uvicorn.lifespan ^
  --hidden-import uvicorn.lifespan.on ^
  --hidden-import uvicorn.protocols.http.auto ^
  --hidden-import uvicorn.protocols.websockets.auto ^
  "%~dp0local_server.py"

echo [build_sidecar] 完成: %~dp0dist\local_sidecar.exe
endlocal
