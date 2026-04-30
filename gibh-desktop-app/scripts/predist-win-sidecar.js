'use strict';

/**
 * 仅在「要打 Windows NSIS」前调用：
 * - 本机为 Windows：执行 PyInstaller，生成 local_sidecar/dist/local_sidecar.exe
 * - 本机为 Linux/macOS：无法在此处生成 .exe；若 dist 中已有 Windows 产物（CI 预置/从 Windows 拷入）则通过，否则失败并说明原因
 */
const fs = require('fs');
const path = require('path');
const { spawnSync } = require('child_process');

const appRoot = path.join(__dirname, '..');
const distDir = path.join(appRoot, 'local_sidecar', 'dist');
const winExe = path.join(distDir, 'local_sidecar.exe');

if (process.env.OMICS_SKIP_WIN_SIDECAR_CHECK === '1' || process.env.OMICS_SKIP_WIN_SIDECAR_CHECK === 'true') {
  console.warn('[predist-win-sidecar] 已跳过检查（OMICS_SKIP_WIN_SIDECAR_CHECK=1）');
  process.exit(0);
}

if (process.platform === 'win32') {
  const bat = path.join(appRoot, 'local_sidecar', 'build_sidecar.bat');
  if (!fs.existsSync(bat)) {
    console.error(`[predist-win-sidecar] 未找到 ${bat}`);
    process.exit(1);
  }
  const r = spawnSync(bat, [], {
    cwd: path.join(appRoot, 'local_sidecar'),
    stdio: 'inherit',
    shell: true,
    env: process.env,
  });
  process.exit(r.status === null ? 1 : r.status);
}

if (fs.existsSync(winExe)) {
  console.log('[predist-win-sidecar] 复用已有 Windows Sidecar:', winExe);
  process.exit(0);
}

const skipWine =
  process.env.OMICS_SKIP_WINE_SIDECAR === '1' || process.env.OMICS_SKIP_WINE_SIDECAR === 'true';

if (!skipWine && process.platform === 'linux') {
  const wineSh = path.join(__dirname, 'build-sidecar-wine.sh');
  if (fs.existsSync(wineSh)) {
    console.log('[predist-win-sidecar] 未找到 local_sidecar.exe，尝试 Wine + 嵌入式 Windows Python 交叉构建（可能需数分钟）…');
    const r = spawnSync('/bin/bash', [wineSh], {
      cwd: appRoot,
      stdio: 'inherit',
      env: process.env,
    });
    if (r.status === 0 && fs.existsSync(winExe)) {
      console.log('[predist-win-sidecar] Wine 交叉构建成功:', winExe);
      process.exit(0);
    }
    console.error('[predist-win-sidecar] Wine 交叉构建未完成或失败，退出码:', r.status);
  }
}

console.error(
  '[predist-win-sidecar] 无法得到 Windows 版 Sidecar（local_sidecar/dist/local_sidecar.exe）。\n' +
    '请任选其一：\n' +
    '  1) 在 Windows 上运行 local_sidecar\\\\build_sidecar.bat；\n' +
    '  2) 在装有 Wine 的 Linux 上手动执行：bash gibh-desktop-app/scripts/build-sidecar-wine.sh；\n' +
    '  3) 由 CI 产出 exe 后放到上述路径；\n' +
    '  4) 不推荐：OMICS_SKIP_WIN_SIDECAR_CHECK=1（安装包内将缺少 Sidecar）。'
);
process.exit(1);
