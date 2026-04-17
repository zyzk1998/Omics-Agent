#!/usr/bin/env bash
# 从「设计与模块.md」提取 Mermaid，用官方 CLI 容器渲染为矢量 SVG 与高倍率 PNG（无需本机安装 Chrome）。
# 依赖：Docker。首次会拉取镜像 MERMAID_IMG（默认 minlag/mermaid-cli:11.4.0）。
# 吉卜力动画风：草木绿/天蓝描边与连线（theme）；画布与节点底色默认纯白。
# 主题 JSON：MERMAID_THEME_CONFIG（默认 scripts/mermaid-design-module-ghibli.json）。
# 画布背景：MERMAID_BG（默认 #ffffff）。
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
MD="${1:-$ROOT/设计与模块.md}"
OUT="${2:-$ROOT/docs/exports}"
MERMAID_IMG="${MERMAID_IMG:-minlag/mermaid-cli:11.4.0}"
MERMAID_THEME_CONFIG="${MERMAID_THEME_CONFIG:-$ROOT/scripts/mermaid-design-module-ghibli.json}"
MERMAID_BG="${MERMAID_BG:-#ffffff}"

mkdir -p "$OUT"
MMD="$OUT/design-module-extract.mmd"
awk '/^```mermaid$/{f=1; next} /^```$/{if(f) exit} f' "$MD" > "$MMD"
if [[ ! -s "$MMD" ]]; then
  echo "未从 $MD 解析到 mermaid 代码块" >&2
  exit 1
fi

run_mmdc() {
  # 与宿主机用户同 uid/gid 写出，避免 root 在挂载卷上写文件导致 EACCES
  # 挂载主题目录（只读），供 mmdc -c 读取 JSON
  local cfg_dir
  cfg_dir="$(dirname "$MERMAID_THEME_CONFIG")"
  local cfg_name
  cfg_name="$(basename "$MERMAID_THEME_CONFIG")"
  docker run --rm -u "$(id -u):$(id -g)" \
    -v "$OUT:/data" \
    -v "$cfg_dir:/mermaid-theme:ro" \
    "$MERMAID_IMG" \
    mmdc -c "/mermaid-theme/$cfg_name" "$@"
}

echo "==> 渲染 SVG（矢量，可无限放大）主题: $MERMAID_THEME_CONFIG"
run_mmdc -i "/data/design-module-extract.mmd" -o "/data/design-module-layered.svg" -b "$MERMAID_BG"

echo "==> 渲染 PNG（-s 倍率 -w 宽度，可按需改环境变量 SCALE、WIDTH）"
SCALE="${SCALE:-3}"
WIDTH="${WIDTH:-2800}"
run_mmdc -i "/data/design-module-extract.mmd" -o "/data/design-module-layered.png" -b "$MERMAID_BG" -s "$SCALE" -w "$WIDTH"

rm -f "$MMD"
echo "完成：$OUT/design-module-layered.svg"
echo "完成：$OUT/design-module-layered.png"
echo "提示：印刷或汇报优先用 SVG；若只要位图，可把 SCALE 改为 4 再运行。"
