'use strict';

const path = require('path');
const fs = require('fs');
const { app, BrowserWindow } = require('electron');

/** 与网页 favicon 同源：仓库内由 services/nginx/html/static/favicon.png 复制为 app-icon.png，打包后随应用分发 */
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

  // Demo 阶段，请将此处的 IP 替换为真实的服务器 Nginx 访问地址
  win.loadURL('http://192.168.32.31:8018');
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
