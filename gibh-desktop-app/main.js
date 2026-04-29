'use strict';
/**
 * 自动更新：electron-updater 的 generic 根地址内嵌自 package.json → build.publish[0].url。
 * 发版前请将该行中的 YOUR_SERVER_URL 换成公网可解析的域名或 IP（勿提交真实内网段入公开仓库时遵守团队 Git 规范）。
 * 也可在环境变量中设置 OMICS_AUTO_UPDATE_BASE_URL（指向 …/downloads/）或 OMICS_AGENT_WEB_URL（将自动追加 /downloads/），优先级高于内嵌 URL。
 */

const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');

(function loadProdEnvFiles() {
  try {
    const dotenv = require('dotenv');
    const candidates = [
      path.join(__dirname, '..', '.env'),
      path.join(__dirname, '..', 'dist-pack.env'),
    ];
    for (const envPath of candidates) {
      if (fs.existsSync(envPath)) {
        dotenv.config({ path: envPath });
        console.log('[Omics Agent] 已加载环境文件:', envPath);
        break;
      }
    }
  } catch (e) {
    console.warn('[Omics Agent] 加载 .env 失败（可忽略）:', e && e.message);
  }
})();

const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
let localSidecarProcess = null;

/** 打包后 extraResources 将仓库 services/nginx/html 拷至 resources/bundled-web-ui（与 electron-builder 配置一致）；当前默认仍 loadURL 远程站点，快照供运维比对或后续改为本地加载方案 */
function logBundledWebUiSnapshotIfPresent() {
  if (!app.isPackaged) return;
  try {
    const snap = path.join(process.resourcesPath, 'bundled-web-ui');
    if (fs.existsSync(snap)) {
      console.log('[Omics Agent] 安装包内附带前端静态快照目录:', snap);
    }
  } catch (_e) {}
}

/** 可选：未安装时仍可启动（仅跳过自动更新） */
let autoUpdater = null;
try {
  ({ autoUpdater } = require('electron-updater'));
} catch (e) {
  console.warn('[Omics Agent] electron-updater 未安装，自动更新已禁用:', e && e.message);
}

/**
 * 向所有 BrowserWindow 的渲染进程广播（主进程无 ipcMain.send → 渲染进程 API，须用 webContents.send）。
 * @param {string} channel
 * @param {Record<string, unknown>} payload
 */
function broadcastAutoUpdate(channel, payload) {
  BrowserWindow.getAllWindows().forEach((win) => {
    if (win.isDestroyed()) return;
    try {
      win.webContents.send(channel, payload);
    } catch (err) {
      console.warn('[Omics Agent] auto-update IPC 发送失败', channel, err && err.message);
    }
  });
}

/**
 * generic 更新包根目录（须能通过 GET 访问到 latest.yml / latest-linux.yml 及安装包与 .blockmap）。
 * 优先环境变量 OMICS_AUTO_UPDATE_BASE_URL；否则若已配置 OMICS_AGENT_WEB_URL，则回退为 {OMICS_AGENT_WEB_URL}/downloads/。
 */
function resolveGenericUpdateBaseUrl() {
  const explicit = process.env.OMICS_AUTO_UPDATE_BASE_URL && String(process.env.OMICS_AUTO_UPDATE_BASE_URL).trim();
  if (explicit) return explicit.replace(/\/?$/, '/');
  const web = process.env.OMICS_AGENT_WEB_URL && String(process.env.OMICS_AGENT_WEB_URL).trim();
  if (web) {
    const base = web.replace(/\/$/, '');
    return `${base}/downloads/`;
  }
  return '';
}

/**
 * 将 electron-updater / Node 网络错误转换为用户可读文案；不把 net::ERR_* 原文推到 UI。
 * 完整错误仅在主进程 console 记录。
 */
function friendlyAutoUpdateMessage(raw) {
  const msg = String((raw && raw.message) || raw || '');
  const lower = msg.toLowerCase();
  if (/err_name_not_resolved|enotfound|getaddrinfo|name not resolved/i.test(lower)) {
    return '无法连接更新服务器，请检查网络或稍后再试；也可在官网下载页获取最新安装包。';
  }
  if (/etimedout|timeout|econnrefused|enetunreach|ehostunreach/i.test(lower)) {
    return '连接更新服务器超时，请稍后再试或前往官网下载最新版。';
  }
  if (/certificate|ssl|tls|unable to verify/i.test(lower)) {
    return '更新通道安全校验失败，请联系管理员或手动从官网下载安装包。';
  }
  if (/404|not found|status code 404/i.test(lower)) {
    return '未找到更新描述文件（latest.yml），请确认服务端已发布对应平台的更新元数据。';
  }
  return '暂时无法检查更新，您可在官网下载页获取最新安装包。';
}

function hookAutoUpdaterEventsOnce() {
  if (!autoUpdater || !app.isPackaged) return;
  if (globalThis.__OMICS_AUTO_UPDATER_EVENTS__) return;
  globalThis.__OMICS_AUTO_UPDATER_EVENTS__ = true;

  const base = resolveGenericUpdateBaseUrl();
  if (base) {
    autoUpdater.setFeedURL({ provider: 'generic', url: base });
    console.log('[Omics Agent] 自动更新 feed（generic）:', base);
  } else {
    console.log(
      '[Omics Agent] 未设置 OMICS_AUTO_UPDATE_BASE_URL / OMICS_AGENT_WEB_URL，将使用 electron-builder 打包时写入的 generic URL（见 gibh-desktop-app/package.json → build.publish[0].url，须将 YOUR_SERVER_URL 替换为真实域名或 IP）'
    );
  }

  autoUpdater.on('checking-for-update', () => {
    broadcastAutoUpdate('app-auto-update-checking', {});
  });

  autoUpdater.on('update-available', (info) => {
    broadcastAutoUpdate('app-auto-update-available', {
      version: info.version,
      releaseDate: info.releaseDate,
      releaseNotes: info.releaseNotes,
    });
  });

  autoUpdater.on('update-not-available', (info) => {
    broadcastAutoUpdate('app-auto-update-not-available', {
      version: (info && info.version) || '',
    });
  });

  autoUpdater.on('error', (err) => {
    const raw = (err && err.message) || String(err);
    console.warn('[Omics Agent] autoUpdater error（详情仅日志）:', raw, err && err.stack);
    broadcastAutoUpdate('app-auto-update-error', {
      userMessage: friendlyAutoUpdateMessage(err),
    });
  });

  autoUpdater.on('download-progress', (p) => {
    broadcastAutoUpdate('app-auto-update-download-progress', {
      percent: p.percent,
      transferred: p.transferred,
      total: p.total,
      bytesPerSecond: p.bytesPerSecond,
    });
  });

  autoUpdater.on('update-downloaded', (info) => {
    broadcastAutoUpdate('app-auto-update-downloaded', {
      version: info.version,
      releaseDate: info.releaseDate,
      releaseNotes: info.releaseNotes,
    });
  });
}

function scheduleAutoUpdateCheck() {
  if (!autoUpdater || !app.isPackaged) return;
  setTimeout(() => {
    autoUpdater
      .checkForUpdatesAndNotify()
      .catch((err) => {
        console.warn('[Omics Agent] checkForUpdatesAndNotify', err && err.message, err && err.stack);
        broadcastAutoUpdate('app-auto-update-error', {
          userMessage: friendlyAutoUpdateMessage(err),
        });
      });
  }, 5000);
}

/** 开发模式默认；生产在 app.isPackaged 下使用 process.env.PROD_SERVER_URL（见仓库根 .env / 打包生成的 dist-pack.env） */
const DEV_SERVER_URL = 'http://127.0.0.1:8018';

const PROD_SERVER_URL_DEFAULT = 'http://127.0.0.1:8018';

/** Windows 任务栏 / 开始菜单分组图标与安装包一致，需与 package.json 的 appId 相同 */
if (process.platform === 'win32') {
  app.setAppUserModelId('com.gibh.agent.demo');
}

/**
 * 窗口 / 任务栏图标（Windows 须 .ico；打包后主进程从 app.asar 内读取，故 package.json build.files 须包含 build/icon.ico）
 */
function getWindowIconPath() {
  if (process.platform === 'win32') {
    const ico = path.join(__dirname, 'build', 'icon.ico');
    if (fs.existsSync(ico)) return ico;
  }
  const png = path.join(__dirname, 'app-icon.png');
  return fs.existsSync(png) ? png : undefined;
}

function getWebBase() {
  const envUrl =
    process.env.OMICS_AGENT_WEB_URL && String(process.env.OMICS_AGENT_WEB_URL).trim()
      ? String(process.env.OMICS_AGENT_WEB_URL).trim().replace(/\/$/, '')
      : '';
  if (envUrl) {
    console.log('[Omics Agent] Web 基址 ← 环境变量 OMICS_AGENT_WEB_URL:', envUrl);
    return envUrl;
  }
  if (app.isPackaged) {
    const raw = process.env.PROD_SERVER_URL || PROD_SERVER_URL_DEFAULT;
    const u = String(raw).trim().replace(/\/$/, '');
    console.log('[Omics Agent] Web 基址 ← 生产（process.env.PROD_SERVER_URL）:', u);
    return u;
  }
  const u = DEV_SERVER_URL.replace(/\/$/, '');
  console.log('[Omics Agent] Web 基址 ← 开发模式（常量 DEV_SERVER_URL）:', u);
  return u;
}

function loadRemote(win) {
  logBundledWebUiSnapshotIfPresent();
  const webBase = getWebBase();
  win.loadURL(webBase).catch((err) => {
    console.error('[Omics Agent] loadURL failed:', err && err.message);
  });
}

function resolvePythonCommand() {
  const explicit = process.env.OMICS_SIDECAR_PYTHON && String(process.env.OMICS_SIDECAR_PYTHON).trim();
  if (explicit) return explicit;
  if (process.platform === 'win32') return 'python';
  return 'python3';
}

function startLocalSidecar() {
  if (localSidecarProcess && !localSidecarProcess.killed) return;
  const scriptPath = path.join(__dirname, 'local_sidecar', 'local_server.py');
  if (!fs.existsSync(scriptPath)) {
    console.error('[Omics Agent] Local Sidecar 启动失败，脚本不存在:', scriptPath);
    return;
  }
  const pythonCmd = resolvePythonCommand();
  localSidecarProcess = spawn(pythonCmd, [scriptPath], {
    cwd: path.join(__dirname, 'local_sidecar'),
    stdio: ['ignore', 'pipe', 'pipe'],
    windowsHide: true,
  });
  console.log('[Omics Agent] Local Sidecar 启动命令:', pythonCmd, scriptPath);

  localSidecarProcess.stdout.on('data', (buf) => {
    console.log('[Local Sidecar]', String(buf || '').trim());
  });
  localSidecarProcess.stderr.on('data', (buf) => {
    console.error('[Local Sidecar]', String(buf || '').trim());
  });
  localSidecarProcess.on('exit', (code, signal) => {
    console.log('[Omics Agent] Local Sidecar 已退出:', { code, signal });
    localSidecarProcess = null;
  });
  localSidecarProcess.on('error', (err) => {
    console.error('[Omics Agent] Local Sidecar 进程异常:', err && err.message);
  });
}

function stopLocalSidecar() {
  if (!localSidecarProcess || localSidecarProcess.killed) return;
  try {
    localSidecarProcess.kill('SIGTERM');
  } catch (e) {
    console.warn('[Omics Agent] Local Sidecar SIGTERM 失败:', e && e.message);
  }
  localSidecarProcess = null;
}

function createMainWindow() {
  const win = new BrowserWindow({
    width: 1280,
    height: 800,
    useContentSize: true,
    autoHideMenuBar: true,
    icon: getWindowIconPath(),
    show: false,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false,
    },
  });

  win.setMenuBarVisibility(false);

  const offlineHtml = path.join(__dirname, 'offline.html');

  win.once('ready-to-show', () => {
    win.show();
  });

  /** 打包后仍拉取远程页面时，Chromium 磁盘缓存可能导致「看得见旧版 UI」；首屏前清空 session 缓存再 loadURL */
  win.webContents.session
    .clearCache()
    .then(() => {
      console.log('[Omics Agent] session.clearCache() 已完成（缓解 HTTP 缓存旧静态资源）');
      loadRemote(win);
    })
    .catch((err) => {
      console.warn('[Omics Agent] session.clearCache() 失败，仍将加载页面:', err && err.message);
      loadRemote(win);
    });

  if (process.env.OMICS_AGENT_DEBUG === '1' || process.env.OMICS_AGENT_DEBUG === 'true') {
    win.webContents.openDevTools({ mode: 'detach' });
  }

  win.webContents.on('did-fail-load', (event, errorCode, errorDescription, validatedURL, isMainFrame) => {
    if (!isMainFrame) return;
    const failed = String(validatedURL || '');
    if (failed.includes('offline.html')) return;
    console.error('[did-fail-load]', errorCode, errorDescription, validatedURL);
    try {
      win.loadFile(offlineHtml, {
        query: {
          u: encodeURIComponent(validatedURL || getWebBase()),
          c: String(errorCode),
        },
      });
    } catch (e) {
      console.error('[offline fallback]', e);
    }
  });

}

if (!globalThis.__OMICS_IPC_REGISTERED__) {
  globalThis.__OMICS_IPC_REGISTERED__ = true;
  ipcMain.on('omics-retry-load', (event) => {
    const w = BrowserWindow.fromWebContents(event.sender);
    if (w && !w.isDestroyed()) loadRemote(w);
  });
  ipcMain.on('app-print-report', async (event) => {
    const win = BrowserWindow.fromWebContents(event.sender);
    if (!win || win.isDestroyed()) return;
    const wc = win.webContents;
    if (!wc || wc.isDestroyed()) return;
    try {
      const { canceled, filePath } = await dialog.showSaveDialog(win, {
        title: '导出分析报告',
        defaultPath: 'Omics_Agent_Report.pdf',
        filters: [{ name: 'PDF 文档', extensions: ['pdf'] }],
      });
      if (canceled || !filePath) return;

      const pdfData = await wc.printToPDF({
        printBackground: true,
        landscape: false,
        pageSize: 'A4',
      });
      fs.writeFileSync(filePath, pdfData);
      try {
        shell.openPath(filePath);
      } catch (openErr) {
        console.warn('[app-print-report] openPath:', openErr && openErr.message);
      }
    } catch (error) {
      console.error('PDF 导出失败:', error);
    }
  });
  ipcMain.on('app-install-update', () => {
    if (!autoUpdater) return;
    try {
      autoUpdater.quitAndInstall(false, true);
    } catch (e) {
      console.error('[app-install-update]', e);
    }
  });
  ipcMain.handle('app-select-workspace-folder', async (event) => {
    const win = BrowserWindow.fromWebContents(event.sender);
    const targetWin = win && !win.isDestroyed() ? win : BrowserWindow.getAllWindows()[0];
    const result = await dialog.showOpenDialog(targetWin || undefined, {
      title: '选择本地项目目录',
      properties: ['openDirectory', 'createDirectory'],
    });
    if (result.canceled || !result.filePaths || result.filePaths.length === 0) {
      return { canceled: true };
    }
    const selectedPath = result.filePaths[0];
    return {
      canceled: false,
      workspace_path: selectedPath,
      workspace_name: path.basename(selectedPath),
    };
  });
}

app.whenReady().then(() => {
  startLocalSidecar();
  hookAutoUpdaterEventsOnce();
  createMainWindow();
  scheduleAutoUpdateCheck();
});

app.on('window-all-closed', () => {
  app.quit();
});

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createMainWindow();
  }
});

app.on('will-quit', () => {
  stopLocalSidecar();
});
