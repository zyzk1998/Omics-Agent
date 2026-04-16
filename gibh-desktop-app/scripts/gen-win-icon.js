'use strict';

/**
 * 从仓库内 app-icon.png 生成 Windows 安装包 / 快捷方式用的 .ico，
 * 避免仅配 PNG 时 electron-builder 报「使用默认 Electron 图标」。
 */
const fs = require('fs');
const path = require('path');
const pngToIco = require('png-to-ico');

const root = path.join(__dirname, '..');
const src = path.join(root, 'app-icon.png');
/** electron-builder 默认从 build/ 读 Windows 图标，文件名用 icon.ico 最稳 */
const outDir = path.join(root, 'build');
const out = path.join(outDir, 'icon.ico');

(async () => {
  if (!fs.existsSync(src)) {
    console.error('缺少 app-icon.png:', src);
    process.exit(1);
  }
  fs.mkdirSync(outDir, { recursive: true });
  const buf = await pngToIco(fs.readFileSync(src));
  fs.writeFileSync(out, buf);
  console.log('已生成', out, '（供 electron-builder Windows / NSIS 使用）');
})().catch((e) => {
  console.error(e);
  process.exit(1);
});
