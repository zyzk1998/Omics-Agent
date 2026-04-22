'use strict';

const path = require('path');
const fs = require('fs');
const { app, BrowserWindow } = require('electron');

/** Windows 任务栏 / 开始菜单分组图标与安装包一致，需与 package.json 的 appId 相同 */
if (process.platform === 'win32') {
  app.setAppUserModelId('com.gibh.agent.demo');
}

/** 仅客户端：gibh-desktop-app/app-icon.png（瘦客户端/任务栏/安装包，与 Web 的 static/favicon.png 解耦；npm run icons:win → build/icon.ico；write-download-release-meta 同步 downloads/client-app-icon.png） */
function getWindowIconPath() {
  const p = path.join(__dirname, 'app-icon.png');
  return fs.existsSync(p) ? p : undefined;
}

function createMainWindow() {
  const win = new BrowserWindow({
    width: 1280,
    height: 800,
    useContentSize: true,
    autoHideMenuBar: true,
    icon: getWindowIconPath(),
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false,
    },
  });

  win.setMenuBarVisibility(false);

  /**
   * 界面与浏览器同源：加载 Nginx 上的 `index.html` / `css/main.css` 等静态资源。
   * 仓库内修改 `services/nginx/html/` 后，部署到该 URL 所指站点即可，瘦客户端无需再拷一份前端文件。
   * 覆盖默认地址：启动前设置环境变量 OMICS_AGENT_WEB_URL（例如 http://127.0.0.1:8018）。
   */
  const webBase =
    process.env.OMICS_AGENT_WEB_URL && String(process.env.OMICS_AGENT_WEB_URL).trim()
      ? String(process.env.OMICS_AGENT_WEB_URL).trim().replace(/\/$/, '')
      : 'http://127.0.0.1:8018';
  win.loadURL(webBase);
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
