#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阵列二 Batch2（Batch2a 8 项 + Batch2b 6 项）写入 skills 表。

用法（仓库根目录）:
  PYTHONPATH=. python3 scripts/sync_phase2_array2_batch2_to_db.py
  docker exec gibh_v2_api /usr/local/bin/python3 /app/scripts/sync_phase2_array2_batch2_to_db.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def _run(script: str) -> None:
    cmd = [sys.executable, str(_ROOT / "scripts" / script)]
    print(f"\n>>> {' '.join(cmd)}")
    subprocess.run(
        cmd,
        check=True,
        cwd=str(_ROOT),
        env={**dict(__import__("os").environ), "PYTHONPATH": str(_ROOT)},
    )


def main() -> int:
    _run("patch_launch_skills_batch.py")
    from gibh_agent.db.connection import SessionLocal
    from gibh_agent.db.seed_skills import run_upsert_system_skills

    db = SessionLocal()
    try:
        n = run_upsert_system_skills(db)
        print(f"\n>>> run_upsert_system_skills: {n} seed rows merged")
    finally:
        db.close()
    print("\nOK: 阵列二 Batch2（14 项）已同步至 skills 表。请在广场「生物医药」「化学」分类查看。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
