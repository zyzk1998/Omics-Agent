# 将 gibh-desktop-app\dist-out 下的 Windows 安装包与 Linux AppImage 复制到 Nginx 静态目录 downloads\，供 /downloads/ 页面分发。
# 用法（在 PowerShell 中，仓库根目录执行）：
#   .\scripts\copy_desktop_installer_to_nginx.ps1
# 或指定仓库根路径：
#   .\scripts\copy_desktop_installer_to_nginx.ps1 -RepoRoot 'D:\omics agent'
param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)
$src = Join-Path $RepoRoot 'gibh-desktop-app\dist-out'
$dst = Join-Path $RepoRoot 'services\nginx\html\downloads'
if (-not (Test-Path -LiteralPath $src)) {
    Write-Error "未找到目录: $src"
    exit 1
}
New-Item -ItemType Directory -Force -Path $dst | Out-Null
Get-ChildItem -LiteralPath $src -File | Where-Object {
    $_.Extension -in '.exe', '.AppImage'
} | ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination $dst -Force
    Write-Host "已复制: $($_.Name) -> $dst"
}
if (-not (Get-ChildItem -LiteralPath $src -Filter '*.exe' -File -ErrorAction SilentlyContinue)) {
    Write-Warning "dist-out 下未找到 .exe（Windows 包需在 Windows 上 npm run dist）"
}
if (-not (Get-ChildItem -LiteralPath $src -Filter '*.AppImage' -File -ErrorAction SilentlyContinue)) {
    Write-Warning "dist-out 下未找到 .AppImage（Linux 包需在 Linux 上 npm run dist 或 npm run dist -- --linux）"
}
