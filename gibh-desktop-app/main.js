'use strict';

const { app, BrowserWindow } = require('electron');

function createMainWindow() {
  const win = new BrowserWindow({
    width: 1280,
    height: 800,
    useContentSize: true,
    autoHideMenuBar: true,
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
