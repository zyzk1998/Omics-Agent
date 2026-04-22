#Requires -Version 5.1
<#
  转发到仓库根目录 scripts/publish-omics-client-from-windows.ps1。
  在 gibh-desktop-app 目录执行： .\scripts\publish-omics-client-from-windows.ps1
  先设置： $env:GIBH_DEPLOY_HOST = "服务器IP或主机名"
#>
# 本文件在 gibh-desktop-app/scripts/ ，仓库根为上两级目录
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Main = Join-Path $RepoRoot "scripts\publish-omics-client-from-windows.ps1"
if (-not (Test-Path -LiteralPath $Main)) {
    Write-Host "未找到: $Main" -ForegroundColor Red
    Write-Host "说明：本脚本只在「完整 GIBH-AGENT-V2 仓库」内有效；当前应在 gibh-desktop-app\scripts 的上两级存在 scripts 目录。" -ForegroundColor Yellow
    exit 1
}
& $Main @args
