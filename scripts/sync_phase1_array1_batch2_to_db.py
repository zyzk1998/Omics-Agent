#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阵列一 Batch2（剩余 13 项 Prompt 软技能）写入 skills 表。

用法（仓库根目录）:
  PYTHONPATH=. python3 scripts/sync_phase1_array1_batch2_to_db.py
  docker exec gibh_v2_api /usr/local/bin/python3 /app/scripts/sync_phase1_array1_batch2_to_db.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def main() -> int:
    cmd = [sys.executable, str(_ROOT / "scripts" / "patch_prompt_soft_skills.py")]
    print(">>>", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(_ROOT), env={**dict(__import__("os").environ), "PYTHONPATH": str(_ROOT)})

    from gibh_agent.db.connection import SessionLocal
    from gibh_agent.db.seed_skills import run_upsert_system_skills

    db = SessionLocal()
    try:
        n = run_upsert_system_skills(db)
        print(f">>> run_upsert_system_skills: {n} seed rows merged")
    finally:
        db.close()
    print("\nOK: 阵列一 Batch2（13 项）已同步。请在广场「生物医药」页签查看。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
