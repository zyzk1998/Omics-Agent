'use strict';

/**
 * 从 app-icon.png 生成 Windows 用的 build/icon.ico。
 *
 * 品牌图应为：圆角白底「卡片」+ 无穷符号主体，**四角透明**（与早期仅主体+透明底
 * 同一思路：图标外轮廓不铺满不透明方块，Explorer 才不会把安装包显示成白正方形）。
 *
 * 本脚本使用 ico-codec 将多尺寸 **PNG 原样嵌入 ICO**（Vista+ 标准），完整保留 RGBA
 * Alpha；不再使用 png-to-ico 的 BMP+DIB 路径，以免透明/半透明在资源管理器呈棋盘格
 * 或被迫整图压成不透明白底变成方块。
 */
const fs = require('fs');
const path = require('path');
const { PNG } = require('pngjs');
const { resize } = require('png-to-ico/lib/png');
const { encodeIco } = require('ico-codec');

const root = path.join(__dirname, '..');
const src = path.join(root, 'app-icon.png');
const outDir = path.join(root, 'build');
const out = path.join(outDir, 'icon.ico');

/** Windows 安装包 / 快捷方式常用尺寸（含 256 供系统缩放） */
const ICO_SIZES = [256, 48, 32, 16];

function main() {
  if (!fs.existsSync(src)) {
    console.error('缺少 app-icon.png:', src);
    process.exit(1);
  }
  fs.mkdirSync(outDir, { recursive: true });

  const png = PNG.sync.read(fs.readFileSync(src));
  if (png.width !== png.height) {
    console.error('app-icon.png 须为正方形，当前为', png.width, 'x', png.height);
    process.exit(1);
  }

  const images = ICO_SIZES.map((s) => {
    const frame = png.width === s ? png : resize(png, s, s);
    const data = PNG.sync.write(frame);
    return { size: s, data: new Uint8Array(data) };
  });

  const icoBytes = encodeIco(images);
  fs.writeFileSync(out, Buffer.from(icoBytes));
  console.log('已生成', out, '（PNG-in-ICO，保留透明外廓与圆角白底）');
}

main();
