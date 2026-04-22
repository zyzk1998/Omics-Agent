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

const { app, BrowserWindow, ipcMain } = require('electron');

/** 开发模式默认；生产在 app.isPackaged 下使用 process.env.PROD_SERVER_URL（见仓库根 .env / 打包生成的 dist-pack.env） */
const DEV_SERVER_URL = 'http://127.0.0.1:8018';

const PROD_SERVER_URL_DEFAULT = 'http://127.0.0.1:8018';

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
