'use strict';

const path = require('path');
const fs = require('fs');
const { app, BrowserWindow, ipcMain } = require('electron');

// ---------------------------------------------------------------------------
// 【打包发行前必改】生产环境（app.isPackaged === true）默认打开的 Omics Web 根地址。
// 下一行常量即为总指挥填写真实 Linux 站点（示例：http://192.168.x.x:8018）的唯一位置。
// 开发环境仍用 DEV_SERVER_URL；任意环境均可被环境变量 OMICS_AGENT_WEB_URL 覆盖。
// ---------------------------------------------------------------------------
const PROD_SERVER_URL = 'http://127.0.0.1:8018';

const DEV_SERVER_URL = 'http://127.0.0.1:8018';

/** Windows 任务栏 / 开始菜单分组图标与安装包一致，需与 package.json 的 appId 相同 */
if (process.platform === 'win32') {
  app.setAppUserModelId('com.gibh.agent.demo');
}

/** 仅客户端：gibh-desktop-app/app-icon.png（瘦客户端与 Web favicon 解耦） */
function getWindowIconPath() {
  const p = path.join(__dirname, 'app-icon.png');
  return fs.existsSync(p) ? p : undefined;
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
    const u = PROD_SERVER_URL.replace(/\/$/, '');
    console.log('[Omics Agent] Web 基址 ← 生产打包（main.js 常量 PROD_SERVER_URL）:', u);
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
}

app.whenReady().then(createMainWindow);

app.on('window-all-closed', () => {
  app.quit();
});

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createMainWindow();
  }
});
