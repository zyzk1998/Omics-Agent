#!/usr/bin/env bash
# 清理 data/results 下历史工作流/测试输出 run_*，释放磁盘（默认保留最近 N 次）。
#
# 用法：
#   ./scripts/cleanup_data_results.sh              # 保留最近 3 次 run
#   KEEP=5 ./scripts/cleanup_data_results.sh       # 保留最近 5 次
#   DRY_RUN=1 ./scripts/cleanup_data_results.sh    # 只打印将删除的路径
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="${ROOT}/data/results"
KEEP="${KEEP:-3}"
DRY_RUN="${DRY_RUN:-0}"

if [ ! -d "$R" ]; then
  echo "目录不存在: $R"
  exit 0
fi

mapfile -t runs < <(ls -1td "$R"/run_* 2>/dev/null || true)
n="${#runs[@]}"
if [ "$n" -le "$KEEP" ]; then
  echo "仅有 ${n} 个 run 目录，不超过 KEEP=${KEEP}，无需清理。"
  exit 0
fi

echo "==> data/results: 共 ${n} 个 run_*，保留最近 ${KEEP} 个，将删除 $((n - KEEP)) 个"
for ((i = KEEP; i < n; i++)); do
  d="${runs[i]}"
  if [ "$DRY_RUN" = "1" ]; then
    echo "DRY_RUN 将删除: $d"
  else
    echo "删除: $d"
    rm -rf "$d"
  fi
done

echo "==> 剩余:"
du -sh "$R" 2>/dev/null || true
ls -1 "$R" 2>/dev/null | wc -l | xargs echo "run 目录数:"
