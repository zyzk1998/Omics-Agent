#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
临时缓存区清理脚本（cron / 手动执行）。

默认清理 RESULTS_DIR 与 UPLOAD_DIR 中超过 7 天的重型文件，
保留 .csv/.md/.png/.json 等轻量报告（与 StorageManager 策略一致）。

示例 crontab（每日 03:00）::

    0 3 * * * cd /app && python3 scripts/cleanup_temp_cache.py >> /var/log/gibh_cleanup.log 2>&1
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.utils.storage_manager import StorageManager  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="清理 GIBH 临时缓存区（/app/results 等）")
    parser.add_argument("--max-days", type=int, default=int(os.getenv("STORAGE_CLEANUP_MAX_DAYS", "7")))
    parser.add_argument("--max-size-gb", type=float, default=float(os.getenv("STORAGE_CLEANUP_MAX_SIZE_GB", "100")))
    parser.add_argument("--deep", action="store_true", help="深度清理（有效保留期缩短）")
    parser.add_argument("--dry-run", action="store_true", help="仅统计，不删除")
    args = parser.parse_args()

    results_dir = os.getenv("RESULTS_DIR", "/app/results")
    upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
    sm = StorageManager(results_dir=results_dir, upload_dir=upload_dir)
    report = sm.auto_cleanup(
        max_days=args.max_days,
        max_size_gb=args.max_size_gb,
        deep=args.deep,
        dry_run=args.dry_run,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not report.get("errors") else 1


if __name__ == "__main__":
    raise SystemExit(main())
