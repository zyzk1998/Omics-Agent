'use strict';

/**
 * 根据 gibh-desktop-app/package.json：
 * 1) 写入 services/nginx/html/downloads/client-release.json（下载页读此文件）
 * 1b) 将 gibh-desktop-app/app-icon.png 复制为 downloads/client-app-icon.png（与站点图标源一致，供静态资源或外链使用）
 * 2) 将 dist-out/ 中与 electron-builder artifactName 一致的产物复制到 downloads/（与 JSON 同步）
 * 3) 当前版本文件已就位时，删除同产品旧文件名，避免磁盘上残留多个 Setup/AppImage 混淆
 *
 * 在 Windows 上执行 npm run dist 后，.exe 会出现在 dist-out，本脚本会一并拷到网页静态目录。
 * 在 Linux 上执行 npm run dist:linux 后，会同步 AppImage；若无本机构建的 .exe 则跳过并提示。
 */
const fs = require('fs');
const path = require('path');

const appRoot = path.join(__dirname, '..');
const pkgPath = path.join(appRoot, 'package.json');
const outPath = path.join(appRoot, '..', 'services', 'nginx', 'html', 'downloads', 'client-release.json');
const distOutDir = path.join(appRoot, 'dist-out');

function expandTemplate(tpl, productName, version, ext) {
  return String(tpl || '')
    .replace(/\$\{productName\}/g, productName)
    .replace(/\$\{version\}/g, version)
    .replace(/\$\{ext\}/g, ext);
}

function displayFromSemver(version) {
  const parts = String(version).split('.');
  if (parts.length >= 2) return 'V' + parts[0] + '.' + parts[1];
  return 'V' + version;
}

function escapeRe(s) {
  return String(s).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

/**
 * @returns {boolean} downloads 目录下目标文件名是否已存在（复制成功或原本就有）
 */
function ensureArtifactInDownloads(filename, downloadsDir) {
  const src = path.join(distOutDir, filename);
  const dest = path.join(downloadsDir, filename);
  if (fs.existsSync(src)) {
    fs.mkdirSync(downloadsDir, { recursive: true });
    fs.copyFileSync(src, dest);
    console.log('已复制到 downloads/', filename);
    return true;
  }
  if (fs.existsSync(dest)) {
    console.log('downloads/ 已存在（跳过复制）', filename);
    return true;
  }
  console.warn('未找到', filename, '：请先在对应平台执行 npm run dist / npm run dist:linux，或将文件放入 dist-out/ 后再运行本脚本。');
  return false;
}

function removeStaleArtifacts(downloadsDir, productName, winSetupFile, linuxAppImageFile) {
  const winRe = new RegExp('^' + escapeRe(productName) + ' Setup .+\\.exe$', 'i');
  const linuxRe = new RegExp('^' + escapeRe(productName) + '-.+\\.AppImage$', 'i');
  let names;
  try {
    names = fs.readdirSync(downloadsDir);
  } catch (e) {
    return;
  }
  const winReady = fs.existsSync(path.join(downloadsDir, winSetupFile));
  const linuxReady = fs.existsSync(path.join(downloadsDir, linuxAppImageFile));
  for (const f of names) {
    if (f === 'client-release.json' || f === 'index.html') continue;
    try {
      if (winReady && winRe.test(f) && f !== winSetupFile) {
        fs.unlinkSync(path.join(downloadsDir, f));
        console.log('已移除旧 Windows 包', f);
      } else if (linuxReady && linuxRe.test(f) && f !== linuxAppImageFile) {
        fs.unlinkSync(path.join(downloadsDir, f));
        console.log('已移除旧 Linux 包', f);
      }
    } catch (e) {
      console.warn('清理旧包时跳过', f, e.message);
    }
  }
}

const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));
const version = pkg.version || '0.0.0';
const productName = (pkg.build && pkg.build.productName) || pkg.productName || pkg.name || 'App';

const winTpl =
  (pkg.build && pkg.build.win && pkg.build.win.artifactName) || '${productName} Setup ${version}.${ext}';
const winSetupFile = expandTemplate(winTpl, productName, version, 'exe');

const linuxTpl =
  (pkg.build && pkg.build.linux && pkg.build.linux.artifactName) || '${productName}-${version}.${ext}';
const linuxAppImageFile = expandTemplate(linuxTpl, productName, version, 'AppImage');

const payload = {
  version,
  displayVersion: displayFromSemver(version),
  productName,
  winSetupFile,
  linuxAppImageFile,
  updatedAt: new Date().toISOString(),
};

const downloadsDir = path.dirname(outPath);
fs.mkdirSync(downloadsDir, { recursive: true });
fs.writeFileSync(outPath, JSON.stringify(payload, null, 2), 'utf8');
console.log('已写入', outPath, payload);

const iconSrc = path.join(appRoot, 'app-icon.png');
const iconDest = path.join(downloadsDir, 'client-app-icon.png');
if (fs.existsSync(iconSrc)) {
  fs.copyFileSync(iconSrc, iconDest);
  console.log('已复制图标', iconDest);
} else {
  console.warn('未找到', iconSrc, '：跳过 client-app-icon.png');
}

ensureArtifactInDownloads(winSetupFile, downloadsDir);
ensureArtifactInDownloads(linuxAppImageFile, downloadsDir);
removeStaleArtifacts(downloadsDir, productName, winSetupFile, linuxAppImageFile);
