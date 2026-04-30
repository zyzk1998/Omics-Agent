'use strict';

/**
 * 跨平台调用 PyInstaller 脚本：Windows → build_sidecar.bat，其他 → build_sidecar.sh
 * 在非 Windows 上仅生成当前平台的 Sidecar 二进制（Linux/macOS ELF/Mach-O）。
 */
const fs = require('fs');
const path = require('path');
const { spawnSync } = require('child_process');

const appRoot = path.join(__dirname, '..');
const localSidecar = path.join(appRoot, 'local_sidecar');

function die(msg, code) {
  console.error(msg);
  process.exit(code ?? 1);
}

function main() {
  if (process.platform === 'win32') {
    const bat = path.join(localSidecar, 'build_sidecar.bat');
    if (!fs.existsSync(bat)) die(`[build-sidecar] 未找到 ${bat}`, 1);
    const r = spawnSync(bat, [], {
      cwd: localSidecar,
      stdio: 'inherit',
      shell: true,
      env: process.env,
    });
    process.exit(r.status === null ? 1 : r.status);
  }

  const sh = path.join(localSidecar, 'build_sidecar.sh');
  if (!fs.existsSync(sh)) die(`[build-sidecar] 未找到 ${sh}`, 1);
  try {
    fs.chmodSync(sh, 0o755);
  } catch (_e) {}
  const r = spawnSync('/bin/bash', [sh], {
    cwd: localSidecar,
    stdio: 'inherit',
    env: process.env,
  });
  process.exit(r.status === null ? 1 : r.status);
}

main();
