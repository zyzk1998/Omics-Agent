'use strict';

const path = require('path');
const fs = require('fs');
const { app, BrowserWindow, ipcMain } = require('electron');

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
  const raw =
    process.env.OMICS_AGENT_WEB_URL && String(process.env.OMICS_AGENT_WEB_URL).trim()
      ? String(process.env.OMICS_AGENT_WEB_URL).trim().replace(/\/$/, '')
      : 'http://127.0.0.1:8018';
  return raw;
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
