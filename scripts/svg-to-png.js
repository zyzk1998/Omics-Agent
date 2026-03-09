#!/usr/bin/env node
/**
 * 将 SVG 转为 PNG
 * 使用 sharp 渲染（需系统安装 librsvg）
 */
const fs = require('fs');
const path = require('path');
const sharp = require('sharp');

async function main() {
  const svgPath = path.resolve(process.argv[2] || path.join(__dirname, '../static_topology.svg'));
  const outPath = process.argv[3] || svgPath.replace(/\.svg$/i, '.png');

  if (!fs.existsSync(svgPath)) {
    console.error('文件不存在:', svgPath);
    process.exit(1);
  }

  const svgBuffer = fs.readFileSync(svgPath);
  await sharp(svgBuffer)
    .png()
    .toFile(outPath);

  console.log('已生成:', outPath);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
