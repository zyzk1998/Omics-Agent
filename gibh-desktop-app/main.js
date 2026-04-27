'use strict';

const path = require('path');
const fs = require('fs');

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
      '[Omics Agent] 未设置 OMICS_AUTO_UPDATE_BASE_URL 且无法从 OMICS_AGENT_WEB_URL 推导 /downloads/，将使用打包时 publish 默认地址（见 package.json build.publish，须替换占位域名）'
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
    broadcastAutoUpdate('app-auto-update-error', {
      message: (err && err.message) || String(err),
      stack: (err && err.stack) || '',
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
        console.error('[Omics Agent] checkForUpdatesAndNotify', err);
        broadcastAutoUpdate('app-auto-update-error', {
          message: (err && err.message) || String(err),
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
  const webBase = getWebBase();
  win.loadURL(webBase).catch((err) => {
    console.error('[Omics Agent] loadURL failed:', err && err.message);
  });
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

  loadRemote(win);

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
}

app.whenReady().then(() => {
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
