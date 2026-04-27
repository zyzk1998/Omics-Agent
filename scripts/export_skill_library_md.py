#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 `gibh_agent/db/seed_skills.py` 中合并后的系统技能清单导出为仓库根目录 `技能库.md`。

不含用户上传的动态插件（`dynamic_skill_plugins`）与仅 DB 存在的非种子技能；
线上完整列表以 `GET /api/skills` 为准。

用法（在仓库根目录）:
    PYTHONPATH=. python3 scripts/export_skill_library_md.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gibh_agent.db.seed_skills import get_all_system_skills_list  # noqa: E402

OUT_MD = REPO_ROOT / "技能库.md"
ROUTE_RE = re.compile(r"\[Skill_Route:\s*(\w+)\s*\]")


def _extract_route(prompt: str) -> str | None:
    if not prompt:
        return None
    m = ROUTE_RE.search(prompt)
    return m.group(1) if m else None


def main() -> None:
    raw = get_all_system_skills_list()
    slim = []
    n_route = 0
    for s in raw:
        pt = s.get("prompt_template") or ""
        route = _extract_route(pt)
        if route:
            n_route += 1
        slim.append(
            {
                "name": (s.get("name") or "").strip(),
                "main_category": (s.get("main_category") or "").strip(),
                "sub_category": (s.get("sub_category") or "").strip(),
                "description": (s.get("description") or "").strip(),
                "skill_route_tool": route,
            }
        )

    n = len(slim)
    body = json.dumps(slim, ensure_ascii=False, indent=2)

    lines = [
        "# 技能库（系统种子技能索引）",
        "",
        "本文档由 `scripts/export_skill_library_md.py` 生成，反映 **`seed_skills.get_all_system_skills_list()`** 合并后的系统预置技能（`author_id=system` 幂等 Upsert 的数据源）。",
        "",
        "**说明：**",
        "",
        "1. **不含**用户 ZIP 上传的 **动态插件技能**（`dynamic_skill_plugins`，广场 `id` 形如 `dyn_*`）；也不含仅存在于数据库、未在种子文件维护的条目。",
        "2. **线上完整列表**（含收藏态、动态插件合并）以 **`GET /api/skills`** 为准（见 `gibh_agent/api/routers/skills.py` / `server.py`）。",
        "3. **`skill_route_tool`**：从 `prompt_template` 中解析的 **`[Skill_Route: …]`** 工具名；与 `@registry.register(name=...)` 对齐时可走语义路由快车道；`null` 表示未声明直连原子工具（多为占位或通用编排）。",
        "4. **扩容规范**见 `docs/技能扩展规范文档.md`；**工具 JSON 全量**见根目录 `工具库.md`（`scripts/export_tool_library_md.py`）。",
        "",
        f"**系统种子技能条数：** {n}",
        f"**其中已声明 [Skill_Route] 的条数：** {n_route}",
        "",
        "---",
        "",
        "## 索引 JSON（name / 分类 / 简介 / skill_route_tool）",
        "",
        "```json",
        body,
        "```",
        "",
        "*生成命令：`PYTHONPATH=. python3 scripts/export_skill_library_md.py`*",
        "",
    ]

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_MD} (skills={n}, with_skill_route={n_route})")


if __name__ == "__main__":
    main()
