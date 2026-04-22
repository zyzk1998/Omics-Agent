#Requires -Version 5.1
<#
  在 gibh-desktop-app 目录下也可运行本文件，将转发到仓库根目录的 scripts/publish-omics-client-from-windows.ps1。
  用法与根目录脚本相同：先设置 $env:GIBH_DEPLOY_HOST，再执行：
    .\publish-client-to-server.ps1
#>
$RepoRoot = Split-Path -Parent $PSScriptRoot
$Main = Join-Path $RepoRoot "scripts\publish-omics-client-from-windows.ps1"
if (-not (Test-Path -LiteralPath $Main)) {
    Write-Host "未找到: $Main" -ForegroundColor Red
    Write-Host '或在 gibh-desktop-app 目录执行: npm run publish:to-server（先设置 $env:GIBH_DEPLOY_HOST）' -ForegroundColor Yellow
    exit 1
}
& $Main @args
