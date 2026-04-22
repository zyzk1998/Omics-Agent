#Requires -Version 5.1
<#
.SYNOPSIS
  在 Windows 本机构建 Omics Agent 瘦客户端 NSIS 安装包，并通过 scp + ssh 上传到 Ubuntu 服务器，
  在服务器上执行 scripts/publish-windows-installer.sh，刷新下载页 client-release.json。

.DESCRIPTION
  首次请在 PowerShell 中设置部署目标（勿将含真实IP的本地配置提交到 Git，可用环境变量）：
    $env:GIBH_DEPLOY_HOST = "你的服务器主机名或IP"
  可选：
    $env:GIBH_DEPLOY_USER = "ubuntu"
    $env:GIBH_REMOTE_REPO = "/home/ubuntu/GIBH-AGENT-V2"   # 服务器上本仓库绝对路径

  执行（任选其一）：
    仓库根目录： .\scripts\publish-omics-client-from-windows.ps1
    gibh-desktop-app 子目录： .\publish-client-to-server.ps1

  依赖：Git、Node.js/npm、OpenSSH 客户端（Windows 10+ 可选功能「OpenSSH 客户端」）。

.NOTES
  请保证本机已对目标主机配置好 SSH 公钥登录，否则 scp/ssh 会提示输入密码。
#>

param(
    [string] $DeployHost = $env:GIBH_DEPLOY_HOST,
    [string] $DeployUser = $(if ($env:GIBH_DEPLOY_USER) { $env:GIBH_DEPLOY_USER } else { "ubuntu" }),
    [string] $RemoteRepo = $(if ($env:GIBH_REMOTE_REPO) { $env:GIBH_REMOTE_REPO } else { "" })
)

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..")).Path
$DesktopDir = Join-Path $RepoRoot "gibh-desktop-app"

if (-not $DeployHost) {
    Write-Host "请先设置环境变量 GIBH_DEPLOY_HOST（部署服务器的主机名或 IP），例如：" -ForegroundColor Red
    Write-Host '  $env:GIBH_DEPLOY_HOST = "192.168.x.x"' -ForegroundColor Yellow
    Write-Host "再运行本脚本。" -ForegroundColor Red
    exit 1
}

if (-not $RemoteRepo) {
    $RemoteRepo = "/home/$DeployUser/GIBH-AGENT-V2"
}

$remoteUpload = "omics-agent-setup-upload.exe"
# scp 目标使用远端用户家目录下的固定文件名（避免空格）
$remoteUploadPath = "${DeployUser}@${DeployHost}:./$remoteUpload"
$uploadAbsRemote = "/home/$DeployUser/$remoteUpload"

foreach ($cmd in @("git", "node", "npm", "scp", "ssh")) {
    if (-not (Get-Command $cmd -ErrorAction SilentlyContinue)) {
        Write-Host "未找到命令: $cmd ，请先安装并加入 PATH（含 OpenSSH 客户端）。" -ForegroundColor Red
        exit 1
    }
}

Write-Host "==> 仓库: $RepoRoot" -ForegroundColor Cyan
Set-Location $RepoRoot
Write-Host "==> git pull ..." -ForegroundColor Cyan
git pull

Set-Location $DesktopDir
Write-Host "==> npm ci ..." -ForegroundColor Cyan
npm ci
Write-Host "==> 生成 icon.ico 并打 Windows 安装包 ..." -ForegroundColor Cyan
npm run icons:win
npx electron-builder --win --x64
if ($LASTEXITCODE -ne 0) {
    throw "electron-builder 失败，退出码 $LASTEXITCODE"
}

$pkgPath = Join-Path $DesktopDir "package.json"
$pkgJson = Get-Content -Raw -Encoding UTF8 $pkgPath | ConvertFrom-Json
$ver = $pkgJson.version
$productName = $pkgJson.build.productName
if (-not $productName) { $productName = "Omics Agent" }
$exeName = "$productName Setup $ver.exe"
$exePath = Join-Path $DesktopDir "dist-out\$exeName"

if (-not (Test-Path -LiteralPath $exePath)) {
    throw "未找到产物: $exePath （请确认 package.json 中 version 与产物文件名一致）"
}

$exeItem = Get-Item -LiteralPath $exePath
if ($exeItem.Length -lt 5MB) {
    throw "安装包体积异常（小于 5MB），可能构建不完整: $($exeItem.Length) bytes"
}

Write-Host "==> 上传 $($exeItem.Length) bytes -> $remoteUploadPath" -ForegroundColor Cyan
& scp -q $exePath $remoteUploadPath
if ($LASTEXITCODE -ne 0) { throw "scp 失败" }

$publishSh = "$RemoteRepo/scripts/publish-windows-installer.sh"
# PS 5.1：避免在双引号串里嵌套 `" && "` 触发解析错误；Unix 路径无空格时用一行 bash -lc
$bashLine = ('chmod +x {0} ; {1} {2}' -f $publishSh, $publishSh, $uploadAbsRemote)

Write-Host "==> ssh 执行服务器端发布脚本 ..." -ForegroundColor Cyan
# 整条作为 ssh 远端一条命令，避免 bashLine 含空格时被拆参数
ssh "${DeployUser}@${DeployHost}" "bash -lc `"$bashLine`""
if ($LASTEXITCODE -ne 0) { throw "ssh / publish-windows-installer.sh 失败" }

Write-Host "完成。请在服务器上 reload nginx（若需要），并浏览器强刷下载页。" -ForegroundColor Green
