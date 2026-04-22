@echo off
REM 一键调用 PowerShell 发布脚本。请先在本窗口或「系统环境变量」中设置 GIBH_DEPLOY_HOST。
REM 示例（cmd）：  set GIBH_DEPLOY_HOST=192.168.x.x
cd /d "%~dp0.."
if "%GIBH_DEPLOY_HOST%"=="" (
  echo [错误] 未设置环境变量 GIBH_DEPLOY_HOST，例如：
  echo   set GIBH_DEPLOY_HOST=你的服务器IP或主机名
  exit /b 1
)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0publish-omics-client-from-windows.ps1" %*
