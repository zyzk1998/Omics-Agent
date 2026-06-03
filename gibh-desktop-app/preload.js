'use strict';

/**
 * 桌面端预加载脚本：向渲染进程暴露安全的 IPC 封装（与 nodeIntegration 共存时挂载到 window）。
 * @see main.js ipcMain handle/on 通道名
 */
const { ipcRenderer } = require('electron');

window.getAppVersion = function getAppVersion() {
  return ipcRenderer.invoke('get-app-version');
};

window.checkForUpdates = function checkForUpdates() {
  ipcRenderer.send('check-for-updates');
};

window.downloadAppUpdate = function downloadAppUpdate() {
  ipcRenderer.send('app-download-update');
};

/** 清空 session 磁盘缓存并忽略缓存重载（「检查更新」已是最新时同步远程前端） */
window.hardReloadFrontend = function hardReloadFrontend() {
  ipcRenderer.send('omics-hard-reload');
};

/**
 * @param {(payload: {
 *   type: 'checking-for-update'|'update-available'|'update-not-available'|'download-progress'|'update-downloaded'|'error',
 *   version?: string,
 *   releaseNotes?: string|Array<{note?: string}|string>,
 *   percent?: number,
 *   devMode?: boolean,
 *   userMessage?: string
 * }) => void} callback
 * @returns {() => void} unsubscribe
 */
window.onUpdaterMessage = function onUpdaterMessage(callback) {
  if (typeof callback !== 'function') {
    return function noop() {};
  }
  const listener = function (_event, payload) {
    callback(payload);
  };
  ipcRenderer.on('updater-message', listener);
  return function unsubscribe() {
    ipcRenderer.removeListener('updater-message', listener);
  };
};
