'use strict';

/**
 * 根据 gibh-desktop-app/package.json：
 * 1) 写入 services/nginx/html/downloads/client-release.json（下载页读此文件）
 * 1b) 将 gibh-desktop-app/app-icon.png 复制为 downloads/client-app-icon.png（**仅客户端**安装包视觉；与 static/favicon.png 网站主 Logo 分离）
 * 2) 将 dist-out/ 中与 electron-builder artifactName 一致的产物复制到 downloads/（与 JSON 同步）
 * 3) 当前版本文件已就位时，删除同产品旧文件名，避免磁盘上残留多个 Setup/AppImage 混淆
 *
 * 在 Windows 上执行 npm run dist / npm run dist:win 后，完整 NSIS .exe（约数十 MB）会出现在 dist-out。
 * 在 Linux 上仅 npm run dist 只会打当前平台产物（多为 AppImage）；Linux 下打的 Windows .exe 若异常偏小（几百 KB），为不完整构建，勿同步到 downloads。
 * 安装包 EXE 内嵌图标来自构建时的 build/icon.ico；网页上的 client-app-icon.png 仅影响站点展示，不能单独替换资源管理器里 exe 的图标。
 * 4) 将 dist-out 中的 latest.yml / latest-linux.yml（若存在）复制到 downloads/，供 electron-updater（generic）拉取元数据。
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
  const parts = String(version)
    .split('.')
    .map((p) => p.trim())
    .filter(Boolean);
  if (parts.length >= 3) return 'V' + parts[0] + '.' + parts[1] + '.' + parts[2];
  if (parts.length === 2) return 'V' + parts[0] + '.' + parts[1];
  return 'V' + String(version).trim();
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

/** electron-updater generic：与安装包同目录的 yml 元数据 */
function copyAutoUpdateMetadataToDownloads(downloadsDir) {
  const names = ['latest.yml', 'latest-linux.yml', 'latest-mac.yml'];
  for (const name of names) {
    const src = path.join(distOutDir, name);
    const dest = path.join(downloadsDir, name);
    if (fs.existsSync(src)) {
      fs.mkdirSync(downloadsDir, { recursive: true });
      fs.copyFileSync(src, dest);
      console.log('已复制到 downloads/', name, '（自动更新元数据）');
    }
  }
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

/** 若站点已手动固定 Windows 安装包文件名（如暂提供旧版 exe），保留 winSetupFile，避免 npm run sync 覆盖。
 * 发布新包时设置环境变量 GIBH_CLIENT_UNPIN_WIN=1（见 scripts/publish-windows-installer.sh），则忽略 pin，按 package.json 版本生成文件名。
 */
let effectiveWin = winSetupFile;
let effectiveLinux = linuxAppImageFile;
const unpinWin =
  process.env.GIBH_CLIENT_UNPIN_WIN === '1' || process.env.GIBH_CLIENT_UNPIN_WIN === 'true';
try {
  if (!unpinWin && fs.existsSync(outPath)) {
    const prev = JSON.parse(fs.readFileSync(outPath, 'utf8'));
    if (prev && prev.winPinned === true && prev.winSetupFile) {
      payload.winSetupFile = prev.winSetupFile;
      payload.winPinned = true;
      effectiveWin = payload.winSetupFile;
      console.log('保留 winPinned：Windows 安装包文件名为', effectiveWin);
    }
  } else if (unpinWin) {
    console.log('GIBH_CLIENT_UNPIN_WIN：按 package.json 写入 winSetupFile，不使用 winPinned');
  }
} catch (e) {
  console.warn('读取旧 client-release.json 合并失败，使用计算值', e.message);
}

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

ensureArtifactInDownloads(effectiveWin, downloadsDir);
ensureArtifactInDownloads(effectiveLinux, downloadsDir);
copyAutoUpdateMetadataToDownloads(downloadsDir);
removeStaleArtifacts(downloadsDir, productName, effectiveWin, effectiveLinux);
