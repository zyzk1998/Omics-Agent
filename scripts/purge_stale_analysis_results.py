#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
清理所有用户「一周前」的分析结果数据（MySQL messages 快照 + SQLite job_history + 磁盘重型产物）。

用法（仓库根目录）:
  PYTHONPATH=. python scripts/purge_stale_analysis_results.py
  PYTHONPATH=. python scripts/purge_stale_analysis_results.py --dry-run
  MAX_DAYS=7 PYTHONPATH=. python scripts/purge_stale_analysis_results.py --yes

与 server.py StorageDaemon 策略一致：默认保留 7 日内数据；轻量报告扩展名不删。
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sqlite3
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("purge_stale_analysis")

# state_snapshot 内与「分析结果/时光机」相关的大字段（保留 text/reasoning/process_log）
_STATE_SNAPSHOT_HEAVY_KEYS = (
    "execution_snapshot",
    "execution_snapshots",
    "steps",
    "workflow",
    "report",
    "data_diagnosis_html",
    "expert_markdown",
    "active_snapshot_id",
)


def _strip_analysis_content(content: Any) -> Tuple[Any, bool]:
    if not isinstance(content, dict):
        return content, False
    changed = False
    if content.pop("execution_snapshot", None) is not None:
        changed = True
    ss = content.get("state_snapshot")
    if isinstance(ss, dict):
        for key in _STATE_SNAPSHOT_HEAVY_KEYS:
            if key in ss:
                ss.pop(key, None)
                changed = True
    tom = content.get("tool_output_memory")
    if tom and isinstance(tom, str) and len(tom) > 2000:
        content["tool_output_memory"] = ""
        changed = True
    return content, changed


def purge_mysql_messages(cutoff: datetime, *, dry_run: bool) -> Dict[str, int]:
    from gibh_agent.db.connection import SessionLocal
    from gibh_agent.db.models import Message

    stats = {"scanned": 0, "updated": 0, "skipped": 0}
    db = SessionLocal()
    try:
        rows = (
            db.query(Message)
            .filter(Message.created_at < cutoff, Message.role == "agent")
            .order_by(Message.id)
            .all()
        )
        stats["scanned"] = len(rows)
        for row in rows:
            new_content, changed = _strip_analysis_content(row.content)
            if not changed:
                stats["skipped"] += 1
                continue
            stats["updated"] += 1
            if dry_run:
                logger.info("[dry-run] 将清理 message id=%s session=%s", row.id, row.session_id)
            else:
                row.content = new_content
        if not dry_run and stats["updated"]:
            db.commit()
        elif dry_run:
            db.rollback()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
    return stats


def purge_sqlite_job_history(cutoff: datetime, db_paths: List[Path], *, dry_run: bool) -> Dict[str, int]:
    stats = {"files": 0, "deleted_rows": 0}
    cutoff_s = cutoff.strftime("%Y-%m-%d %H:%M:%S")
    for p in db_paths:
        if not p.is_file():
            logger.info("SQLite 跳过（不存在）: %s", p)
            continue
        stats["files"] += 1
        conn = sqlite3.connect(str(p))
        try:
            cur = conn.cursor()
            if dry_run:
                cur.execute(
                    "SELECT COUNT(*) FROM job_history WHERE created_at < ?",
                    (cutoff_s,),
                )
                n = cur.fetchone()[0]
                logger.info("[dry-run] %s job_history 将删除 %s 行 (< %s)", p, n, cutoff_s)
                stats["deleted_rows"] += int(n or 0)
            else:
                cur.execute("DELETE FROM job_history WHERE created_at < ?", (cutoff_s,))
                stats["deleted_rows"] += int(cur.rowcount or 0)
                conn.commit()
                logger.info("%s job_history 已删除 %s 行", p, cur.rowcount)
        finally:
            conn.close()
    return stats


def purge_filesystem(
    max_days: int,
    results_dir: Path,
    upload_dir: Path,
    *,
    dry_run: bool,
) -> Dict[str, Any]:
    from gibh_agent.utils.storage_manager import StorageManager

    sm = StorageManager(results_dir=str(results_dir), upload_dir=str(upload_dir))
    storage_result = sm.auto_cleanup(max_days=max_days, max_size_gb=100.0, deep=False, dry_run=dry_run)

    run_deleted = 0
    run_freed = 0
    cutoff_ts = time.time() - max_days * 86400
    runs_root = results_dir
    if runs_root.is_dir():
        for p in runs_root.glob("run_*"):
            if not p.is_dir():
                continue
            try:
                mtime = p.stat().st_mtime
            except OSError:
                continue
            if mtime >= cutoff_ts:
                continue
            try:
                size = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
            except OSError:
                size = 0
            run_deleted += 1
            run_freed += size
            if dry_run:
                logger.info("[dry-run] 将删除 run 目录: %s", p)
            else:
                import shutil

                shutil.rmtree(p, ignore_errors=True)
                logger.info("已删除 run 目录: %s", p)

    storage_result["run_dirs_deleted"] = run_deleted
    storage_result["run_dirs_freed_mb"] = round(run_freed / (1024 * 1024), 4)
    return storage_result


def main() -> int:
    parser = argparse.ArgumentParser(description="清理一周前所有用户的分析结果数据")
    parser.add_argument("--max-days", type=int, default=int(os.getenv("MAX_DAYS", "7")))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--yes", action="store_true", help="非 dry-run 时跳过确认")
    args = parser.parse_args()

    max_days = max(1, args.max_days)
    cutoff = datetime.utcnow() - timedelta(days=max_days)
    dry_run = args.dry_run

    upload_dir = Path(os.getenv("UPLOAD_DIR", str(ROOT / "data" / "uploads")))
    results_dir = Path(os.getenv("RESULTS_DIR", str(ROOT / "data" / "results")))
    sqlite_paths = [
        ROOT / "workflows.db",
        ROOT / "data" / "workflows.db",
    ]

    logger.info("=== 分析结果清理 max_days=%s cutoff(UTC)=%s dry_run=%s ===", max_days, cutoff.isoformat(), dry_run)
    if not dry_run and not args.yes:
        print(f"将清理早于 {cutoff.isoformat()} UTC 的分析结果（约 {max_days} 天前）。加 --yes 执行。")
        return 1

    mysql_stats = purge_mysql_messages(cutoff, dry_run=dry_run)
    sqlite_stats = purge_sqlite_job_history(cutoff, sqlite_paths, dry_run=dry_run)
    fs_stats = purge_filesystem(max_days, results_dir, upload_dir, dry_run=dry_run)

    report = {
        "cutoff_utc": cutoff.isoformat(),
        "max_days": max_days,
        "dry_run": dry_run,
        "mysql_messages": mysql_stats,
        "sqlite_job_history": sqlite_stats,
        "filesystem": fs_stats,
        "paths": {
            "upload_dir": str(upload_dir.resolve()),
            "results_dir": str(results_dir.resolve()),
        },
    }
    print("\n========== 清理汇报 ==========")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
