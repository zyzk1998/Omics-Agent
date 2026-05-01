'use strict';
/**
 * 自动更新：electron-updater 的 generic 根地址内嵌自 package.json → build.publish[0].url。
 * 发版前请将该行中的 YOUR_SERVER_URL 换成公网可解析的域名或 IP（勿提交真实内网段入公开仓库时遵守团队 Git 规范）。
 * 也可在环境变量中设置 OMICS_AUTO_UPDATE_BASE_URL（指向 …/downloads/）或 OMICS_AGENT_WEB_URL（将自动追加 /downloads/），优先级高于内嵌 URL。
 *
 * 本地大文件/附件分流：业务逻辑在打包的 Web UI（services/nginx/html/index.html）的 sendMessage 前闸
 * `runLocalAttachmentSmartGate` + Sidecar `POST /api/tools/upload_to_cloud`；此处无重复实现。
 */

const path = require('path');
const fs = require('fs');
const http = require('http');
const { spawn, execSync } = require('child_process');

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

/** 禁止多开：第二实例直接退出，避免重复拉起 Sidecar 导致 8019 端口冲突（Windows: Errno 10048） */
const gotSingleInstanceLock = app.requestSingleInstanceLock();
if (!gotSingleInstanceLock) {
  app.quit();
} else {
  app.on('second-instance', () => {
    const wins = BrowserWindow.getAllWindows();
    if (wins.length && wins[0] && !wins[0].isDestroyed()) {
      const w = wins[0];
      if (w.isMinimized()) w.restore();
      w.show();
      w.focus();
    }
  });
}

let localSidecarProcess = null;

/** Sidecar 进程 stdout/stderr 环形缓冲上限（字符） */
const SIDECAR_LOG_CAP = 48000;

/** 同一启动周期内 Sidecar 致命错误只弹一次原生对话框，避免刷屏 */
let sidecarFatalDialogShown = false;

function appendSidecarLog(buf, chunk, cap) {
  const s = String(chunk || '');
  const next = buf + s;
  if (next.length <= cap) return next;
  return next.slice(-cap);
}

function truncateForErrorBox(text, maxLen) {
  const t = String(text || '').replace(/\r\n/g, '\n').trim();
  if (t.length <= maxLen) return t;
  return `${t.slice(0, maxLen)}\n\n…（输出过长已截断）`;
}

/** Windows Errno 10048 / EADDRINUSE 等：Sidecar 重复绑定 127.0.0.1:8019 */
function isPortBindConflict(stderrBuf, stdoutBuf) {
  const s = `${String(stderrBuf || '')}\n${String(stdoutBuf || '')}`;
  return /10048|EADDRINUSE|address already in use|Only one usage of each socket|error while attempting to bind/i.test(
    s
  );
}

function showSidecarFatalDialog(title, detail) {
  if (sidecarFatalDialogShown) return;
  sidecarFatalDialogShown = true;
  try {
    dialog.showErrorBox(title, truncateForErrorBox(detail, 7800));
  } catch (e) {
    console.error('[Omics Agent] showErrorBox 失败:', e && e.message);
  }
}

/**
 * 跨平台 Python 启动候选：优先 OMICS_SIDECAR_PYTHON；Windows 常见 python / py -3 / python3；类 Unix 常见 python3 / python。
 * @returns {string[][]} 每个元素为 [command, ...fixedArgs]，再与 scriptPath 拼接为 spawn 参数。
 */
function getPythonSpawnCandidates() {
  const explicit = process.env.OMICS_SIDECAR_PYTHON && String(process.env.OMICS_SIDECAR_PYTHON).trim();
  if (explicit) return [[explicit]];
  if (process.platform === 'win32') {
    return [['python'], ['py', '-3'], ['python3']];
  }
  return [['python3'], ['python']];
}

function formatSidecarSpawnLabel(spec) {
  return spec.join(' ');
}

/**
 * 生产包：PyInstaller 单文件路径（extraResources → resources/local_sidecar/）
 */
function resolvePackagedSidecarBinary(resourcesDir) {
  const base = path.join(resourcesDir, 'local_sidecar');
  if (process.platform === 'win32') {
    return path.join(base, 'local_sidecar.exe');
  }
  return path.join(base, 'local_sidecar');
}

/**
 * Sidecar 子进程日志与异常（开发态 Python 与生产态二进制共用）
 * @param {import('child_process').ChildProcess} child
 * @param {string} displayLabel 错误弹窗中展示的命令行说明
 * @param {{ onENOENT?: () => void }} options Python 多候选时传入 onENOENT；打包二进制勿传
 */
function hookSidecarProcess(child, displayLabel, options) {
  let stdoutBuf = '';
  let stderrBuf = '';
  let skipExitDialog = false;

  function attachStream(stream, kind) {
    if (!stream) return;
    stream.on('data', (chunk) => {
      if (kind === 'stdout') {
        stdoutBuf = appendSidecarLog(stdoutBuf, chunk, SIDECAR_LOG_CAP);
      } else {
        stderrBuf = appendSidecarLog(stderrBuf, chunk, SIDECAR_LOG_CAP);
      }
      const line = String(chunk || '').trim();
      if (line) {
        if (kind === 'stdout') console.log('[Local Sidecar stdout]', line);
        else console.error('[Local Sidecar stderr]', line);
      }
    });
  }
  attachStream(child.stdout, 'stdout');
  attachStream(child.stderr, 'stderr');

  child.on('error', (err) => {
    if (err && err.code === 'ENOENT' && options && typeof options.onENOENT === 'function') {
      skipExitDialog = true;
      console.warn('[Omics Agent] 未找到解释器，尝试下一候选');
      if (localSidecarProcess === child) localSidecarProcess = null;
      options.onENOENT();
      return;
    }
    const detail = `命令：${displayLabel}\n${err && err.message ? err.message : String(err)}`;
    console.error('[Omics Agent] Local Sidecar spawn error:', detail);
    showSidecarFatalDialog('本地 Sidecar 启动失败', detail);
    if (localSidecarProcess === child) localSidecarProcess = null;
  });

  child.on('close', (code, signal) => {
    console.log('[Omics Agent] Local Sidecar 已关闭:', { code, signal });
    if (localSidecarProcess === child) localSidecarProcess = null;
    if (child._omicsIntentionalKill || skipExitDialog) return;
    const abnormal = code !== 0 && code !== null;
    if (!abnormal) return;
    const detail = `命令：${displayLabel}\n退出码：${code}${signal ? `，信号：${signal}` : ''}\n\n--- stderr ---\n${stderrBuf || '(无)'}\n\n--- stdout ---\n${stdoutBuf || '(无)'}`;
    console.error('[Omics Agent] Local Sidecar 异常退出（完整日志见上）');
    if (isPortBindConflict(stderrBuf, stdoutBuf)) {
      showSidecarFatalDialog(
        '本地 Sidecar 无法绑定端口',
        `端口 8019 已被占用（常见于已有一个 Omics Agent / Sidecar 在运行，或其它程序占用该端口）。\n\n请只保留一个客户端，或在任务管理器中结束重复的 Omics Agent / local_sidecar 进程后再启动。\n\n--- 技术详情 ---\n${detail}`
      );
    } else {
      showSidecarFatalDialog('本地 Sidecar 异常退出', detail);
    }
  });

  localSidecarProcess = child;
}

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

function trySpawnSidecarCandidate(candidates, index, scriptPath, sidecarCwd, failedENOENTLabels) {
  if (index >= candidates.length) {
    const tried = failedENOENTLabels.length ? failedENOENTLabels.join('\n') : '(无)';
    showSidecarFatalDialog(
      '本地 Sidecar 无法启动',
      `未找到可用的 Python 解释器。已尝试：\n${tried}\n\n请安装 Python 3，或将环境变量 OMICS_SIDECAR_PYTHON 设为解释器的完整路径（例如 C:\\Python311\\python.exe）。`
    );
    return;
  }

  const spec = candidates[index];
  const cmd = spec[0];
  const prefixArgs = spec.slice(1);
  const spawnArgs = [...prefixArgs, scriptPath];

  console.log('[Omics Agent] Local Sidecar（开发态 Python）尝试启动:', cmd, spawnArgs);

  const child = spawn(cmd, spawnArgs, {
    cwd: sidecarCwd,
    stdio: ['ignore', 'pipe', 'pipe'],
    windowsHide: true,
  });

  const displayLabel = `${formatSidecarSpawnLabel(spec)} ${scriptPath}`;
  hookSidecarProcess(child, displayLabel, {
    onENOENT: () => {
      const label = `${formatSidecarSpawnLabel(spec)} （系统未找到该可执行文件）`;
      console.warn('[Omics Agent]', label);
      failedENOENTLabels.push(label);
      trySpawnSidecarCandidate(candidates, index + 1, scriptPath, sidecarCwd, failedENOENTLabels);
    },
  });
}

function startLocalSidecar() {
  if (localSidecarProcess && !localSidecarProcess.killed) return;
  const resourcesDir = app.isPackaged ? process.resourcesPath : __dirname;

  if (app.isPackaged) {
    const binaryPath = resolvePackagedSidecarBinary(resourcesDir);
    if (!fs.existsSync(binaryPath)) {
      showSidecarFatalDialog(
        '本地 Sidecar 无法启动',
        `未找到打包的 Sidecar 可执行文件：\n${binaryPath}\n请在 **当前目标平台** 先执行 local_sidecar/build_sidecar（生成 dist/）再打 Electron 安装包。`
      );
      return;
    }
    const sidecarCwd = path.dirname(binaryPath);
    console.log('[Omics Agent] Local Sidecar（生产 PyInstaller 二进制）:', binaryPath);
    const child = spawn(binaryPath, [], {
      cwd: sidecarCwd,
      stdio: ['ignore', 'pipe', 'pipe'],
      windowsHide: true,
    });
    hookSidecarProcess(child, binaryPath, {});
    return;
  }

  const scriptPath = path.join(__dirname, 'local_sidecar', 'local_server.py');
  const sidecarCwd = path.join(__dirname, 'local_sidecar');
  if (!fs.existsSync(scriptPath)) {
    showSidecarFatalDialog(
      '本地 Sidecar 无法启动',
      `未找到开发态入口脚本：\n${scriptPath}\n请确认 gibh-desktop-app/local_sidecar 目录存在。`
    );
    return;
  }
  const candidates = getPythonSpawnCandidates();
  trySpawnSidecarCandidate(candidates, 0, scriptPath, sidecarCwd, []);
}

function stopLocalSidecar() {
  if (!localSidecarProcess || localSidecarProcess.killed) return;
  const p = localSidecarProcess;
  p._omicsIntentionalKill = true;
  try {
    p.kill('SIGTERM');
  } catch (e) {
    console.warn('[Omics Agent] Local Sidecar SIGTERM 失败:', e && e.message);
  }
  localSidecarProcess = null;
}

/**
 * 若 8019 上已有本应用 Sidecar（如残留进程），则不再 spawn，避免端口冲突。
 * 与 local_sidecar/local_server.py 的 GET /health 约定一致。
 */
function probeExistingSidecarHealth() {
  return new Promise((resolve) => {
    const req = http.get(
      'http://127.0.0.1:8019/health',
      { timeout: 1500 },
      (res) => {
        let data = '';
        res.on('data', (c) => {
          data += String(c);
        });
        res.on('end', () => {
          if (res.statusCode !== 200) {
            resolve(false);
            return;
          }
          try {
            const j = JSON.parse(data);
            resolve(Boolean(j && j.status === 'ok' && j.service === 'omics-local-sidecar'));
          } catch (_e) {
            resolve(false);
          }
        });
      }
    );
    req.on('error', () => resolve(false));
    req.on('timeout', () => {
      try {
        req.destroy();
      } catch (_e) {}
      resolve(false);
    });
  });
}

/** 轮询直至 /health 成功或超时（spawn 后需等待 uvicorn 绑定端口） */
async function waitForSidecarReady(timeoutMs) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    if (await probeExistingSidecarHealth()) return true;
    await new Promise((r) => setTimeout(r, 280));
  }
  return false;
}

/** 退出前同步请求 Sidecar 自毁（优先 curl，避免 async before-quit 未被 Electron 等待） */
function requestSidecarShutdownBestEffortSync() {
  try {
    if (process.platform === 'win32') {
      execSync('curl.exe -s -m 2 -X POST http://127.0.0.1:8019/api/shutdown', {
        stdio: 'ignore',
        windowsHide: true,
      });
    } else {
      execSync('curl -s -m 2 -X POST http://127.0.0.1:8019/api/shutdown', { stdio: 'ignore' });
    }
  } catch (_) {}
}

/**
 * 硬核兜底：杀掉本客户端 spawn 的 Sidecar 进程树（Windows taskkill /T /F；类 Unix kill -9）。
 * 解决 Electron 退出后 Python 子进程残留导致旧 Sidecar 长期占用 8019。
 */
function forceKillSidecarProcessTree() {
  if (!localSidecarProcess || localSidecarProcess.killed || !localSidecarProcess.pid) {
    localSidecarProcess = null;
    return;
  }
  const pid = localSidecarProcess.pid;
  try {
    if (process.platform === 'win32') {
      execSync(`taskkill /PID ${pid} /T /F`, { stdio: 'ignore', windowsHide: true });
    } else {
      try {
        execSync(`kill -9 ${pid}`, { stdio: 'ignore' });
      } catch (_e) {
        try {
          process.kill(pid, 'SIGKILL');
        } catch (__e) {}
      }
    }
  } catch (_err) {
    // 进程已退出等情况忽略
  }
  localSidecarProcess = null;
}

async function startLocalSidecarIfNeeded() {
  if (localSidecarProcess && !localSidecarProcess.killed) {
    await waitForSidecarReady(2500);
    return;
  }
  const ok = await probeExistingSidecarHealth();
  if (ok) {
    console.log('[Omics Agent] 127.0.0.1:8019 已有 Sidecar 健康响应，跳过重复启动');
    return;
  }
  startLocalSidecar();
  await waitForSidecarReady(14000);
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

  /** 供渲染进程 localSidecarUrl() 使用；仅允许 127.0.0.1 / localhost，防止 OMICS_SIDECAR_URL 误指远程站导致 check_file 打到 Nginx 404。 */
  win.webContents.on('did-finish-load', () => {
    try {
      const fallback = 'http://127.0.0.1:8019';
      const raw = String(process.env.OMICS_SIDECAR_URL || fallback).trim().replace(/\/+$/, '');
      let sidecarBase = fallback;
      try {
        const withScheme = /^https?:\/\//i.test(raw) ? raw : `http://${raw}`;
        const u = new URL(withScheme);
        if (u.hostname === '127.0.0.1' || u.hostname === 'localhost') {
          sidecarBase = /^https?:\/\//i.test(raw) ? raw : `${u.protocol}//${u.host}`;
        }
      } catch (_) {
        sidecarBase = fallback;
      }
      win.webContents.executeJavaScript(
        `try{window.OMICS_SIDECAR_BASE_URL=${JSON.stringify(sidecarBase)};}catch(_){}`
      );
    } catch (e) {
      console.warn('[Omics Agent] inject OMICS_SIDECAR_BASE_URL failed:', e && e.message);
    }
  });

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

if (gotSingleInstanceLock) {
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

  app.whenReady().then(async () => {
    await startLocalSidecarIfNeeded();
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
    requestSidecarShutdownBestEffortSync();
    forceKillSidecarProcessTree();
    stopLocalSidecar();
  });
}
