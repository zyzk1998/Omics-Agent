'use strict';

/**
 * 从仓库根 .env 读取 PROD_SERVER_URL，写入 gibh-desktop-app/dist-pack.env（仅一行）。
 * 供 electron-builder extraResources 打进安装包，避免把完整 .env（含密钥）打入 resources。
 */
const fs = require('fs');
const path = require('path');

const repoRoot = path.join(__dirname, '..', '..');
const rootEnvPath = path.join(repoRoot, '.env');
const outPath = path.join(__dirname, '..', 'dist-pack.env');

const DEFAULT_PROD = 'http://127.0.0.1:8018';

function readProdServerUrl(envPath) {
  if (!fs.existsSync(envPath)) {
    console.warn('[prepare-pack-env] 未找到', envPath, '→ 使用默认', DEFAULT_PROD);
    return DEFAULT_PROD;
  }
  const txt = fs.readFileSync(envPath, 'utf8');
  const m = txt.match(/^\s*PROD_SERVER_URL\s*=\s*(.+)$/m);
  if (!m) {
    console.warn('[prepare-pack-env]', envPath, '中无 PROD_SERVER_URL → 使用默认', DEFAULT_PROD);
    return DEFAULT_PROD;
  }
  return m[1].trim().replace(/^["']|["']$/g, '');
}

const prod = readProdServerUrl(rootEnvPath);
fs.writeFileSync(outPath, `PROD_SERVER_URL=${prod}\n`, 'utf8');
console.log('[prepare-pack-env] 已写入', outPath);
