@echo off
REM 在 gibh-desktop-app 目录双击或运行；需先 set GIBH_DEPLOY_HOST=服务器IP
cd /d "%~dp0"
if "%GIBH_DEPLOY_HOST%"=="" (
  echo [错误] 请先执行: set GIBH_DEPLOY_HOST=你的服务器IP
  exit /b 1
)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0publish-client-to-server.ps1" %*
