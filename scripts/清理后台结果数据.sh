#!/usr/bin/env bash
# 后台分析结果数据清理（默认：所有用户、7 天前）
#
# 包含：
#   1) MySQL messages 内 execution_snapshot / 重型 state_snapshot 等
#   2) SQLite job_history 过期行
#   3) data/uploads、data/results 重型文件（StorageManager，与 API 定时任务同策略）
#   4) data/results/run_* 目录（修改时间早于 MAX_DAYS）
#
# 用法：
#   ./scripts/清理后台结果数据.sh              # 执行清理（需确认）
#   ./scripts/清理后台结果数据.sh --yes        # 跳过确认
#   DRY_RUN=1 ./scripts/清理后台结果数据.sh    # 仅统计
#   MAX_DAYS=14 ./scripts/清理后台结果数据.sh --yes
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export MAX_DAYS="${MAX_DAYS:-7}"
export UPLOAD_DIR="${UPLOAD_DIR:-$ROOT/data/uploads}"
export RESULTS_DIR="${RESULTS_DIR:-$ROOT/data/results}"
export MYSQL_HOST="${MYSQL_HOST:-127.0.0.1}"
export MYSQL_PORT="${MYSQL_PORT:-3306}"

DRY_RUN="${DRY_RUN:-0}"
YES_FLAG=()
if [[ "${1:-}" == "--yes" ]]; then
  YES_FLAG=(--yes)
  shift
fi

PY_ARGS=(scripts/purge_stale_analysis_results.py --max-days "$MAX_DAYS")
if [[ "$DRY_RUN" == "1" ]]; then
  PY_ARGS+=(--dry-run)
else
  PY_ARGS+=("${YES_FLAG[@]}")
fi

echo "==> 启动分析结果清理 MAX_DAYS=${MAX_DAYS} DRY_RUN=${DRY_RUN}"
PYTHONPATH=. python3 "${PY_ARGS[@]}"
