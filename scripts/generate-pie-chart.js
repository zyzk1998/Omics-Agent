#!/usr/bin/env node
/**
 * 生成多组学市场份额饼图 PNG
 * 纯 SVG 实现，不依赖浏览器
 */
const fs = require('fs');
const path = require('path');
const sharp = require('sharp');

const MARKET_DATA = [
  { name: '基因组学', value: 35, color: '#2563EB' },
  { name: '影像组学', value: 20, color: '#7C3AED' },
  { name: '蛋白质组学', value: 15, color: '#059669' },
  { name: '转录组学', value: 12, color: '#DB2777' },
  { name: '空间组学', value: 8, color: '#EA580C', isFastest: true },
  { name: '代谢组学', value: 6, color: '#D97706' },
  { name: '表观基因组学', value: 4, color: '#DC2626' },
];

function polarToCartesian(cx, cy, r, angle) {
  const rad = ((angle - 90) * Math.PI) / 180;
  return { x: cx + r * Math.cos(rad), y: cy + r * Math.sin(rad) };
}

function describeArc(cx, cy, r, startAngle, endAngle) {
  const start = polarToCartesian(cx, cy, r, endAngle);
  const end = polarToCartesian(cx, cy, r, startAngle);
  const largeArc = endAngle - startAngle <= 180 ? 0 : 1;
  return `M ${start.x} ${start.y} A ${r} ${r} 0 ${largeArc} 0 ${end.x} ${end.y} L ${cx} ${cy} Z`;
}

function generateSvg() {
  const w = 720;
  const h = 420;
  const cx = 200;
  const cy = 210;
  const outerR = 100;
  const innerR = 60;
  const total = MARKET_DATA.reduce((s, d) => s + d.value, 0);

  let paths = '';
  let startAngle = 0;
  const RADIAN = Math.PI / 180;

  MARKET_DATA.forEach((item) => {
    const angle = (item.value / total) * 360;
    const endAngle = startAngle + angle;
    const midAngle = startAngle + angle / 2;
    const r2 = innerR + (outerR - innerR) * 0.5;
    const tx = cx + r2 * Math.cos(-midAngle * RADIAN);
    const ty = cy + r2 * Math.sin(-midAngle * RADIAN);

    const outerStart = polarToCartesian(cx, cy, outerR, startAngle);
    const outerEnd = polarToCartesian(cx, cy, outerR, endAngle);
    const innerStart = polarToCartesian(cx, cy, innerR, endAngle);
    const innerEnd = polarToCartesian(cx, cy, innerR, startAngle);
    const largeArc = angle <= 180 ? 0 : 1;
    const pathD = `M ${outerStart.x} ${outerStart.y} A ${outerR} ${outerR} 0 ${largeArc} 0 ${outerEnd.x} ${outerEnd.y} L ${innerEnd.x} ${innerEnd.y} A ${innerR} ${innerR} 0 ${largeArc} 1 ${innerStart.x} ${innerStart.y} Z`;

    paths += `<path d="${pathD}" fill="${item.color}" stroke="white" stroke-width="2"/>`;
    paths += `<text x="${tx}" y="${ty}" text-anchor="middle" dominant-baseline="central" fill="white" font-size="12" font-weight="600">${item.value}%</text>`;
    startAngle = endAngle;
  });

  const legendX = 380;
  const legendItems = MARKET_DATA.map(
    (d, i) => {
      const y = 80 + i * 28;
      return `<g transform="translate(${legendX},${y})">
        <circle cx="6" cy="6" r="6" fill="${d.color}"/>
        <text x="20" y="10" font-size="13" fill="#475569">${d.name}${d.isFastest ? '（增长最快）' : ''}</text>
      </g>`;
    }
  ).join('');

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" font-family="PingFang SC, Microsoft YaHei, sans-serif">
  <rect width="${w}" height="${h}" fill="#ffffff"/>
  <text x="${w/2}" y="28" text-anchor="middle" font-size="18" font-weight="600" fill="#1e293b">2025 全球多组学与 AI 影像市场份额</text>
  <g transform="translate(${cx},${cy})">${paths}</g>
  <text x="${cx}" y="-30" text-anchor="middle" font-size="28" font-weight="700" fill="#1e293b">$50B+</text>
  <text x="${cx}" y="-8" text-anchor="middle" font-size="12" fill="#64748b">市场规模</text>
  <g>${legendItems}</g>
</svg>`;
}

async function main() {
  const outDir = path.resolve(process.argv[2] || path.join(__dirname, '../output-images'));
  const outPath = path.join(outDir, 'omics-market-share.png');

  const svg = generateSvg();
  fs.mkdirSync(outDir, { recursive: true });

  await sharp(Buffer.from(svg))
    .png()
    .toFile(outPath);

  console.log('已生成:', outPath);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
