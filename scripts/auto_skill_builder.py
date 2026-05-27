#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skill Factory — 从 CSV 批量调用 LLM 生成 BaseSkill 技能模块。

用法（仓库根目录）:
  PYTHONPATH=. python scripts/auto_skill_builder.py \\
    --csv docs/生物信息分析工具.csv --limit 5 --dry-run

  PYTHONPATH=. python scripts/auto_skill_builder.py \\
    --csv docs/生物信息分析工具.csv --start 0 --limit 50

环境变量:
  SKILL_FACTORY_MODEL   默认 deepseek-chat（走 MODEL_ROUTING_TABLE）
  DEEPSEEK_API_KEY / ZHIPU_API_KEY / MOONSHOT_API_KEY  至少其一
  SKILL_FACTORY_OUT_DIR 默认 gibh_agent/skills
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.skill_factory_common import (  # noqa: E402
    DEFAULT_FACTORY_MODEL,
    extract_dependencies_from_source,
    extract_python_from_llm,
    factory_llm_client,
    output_filename_for_skill_id,
    parse_csv_row,
    row_to_prompt_block,
    skill_id_from_row,
    write_skill_dependency_report,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

SKILL_CODEGEN_SYSTEM_PROMPT = """你是一个资深 Python 生信 MLOps 工程师与 GIBH-AGENT 技能代码生成器。

请根据用户提供的工具元数据，生成**单个完整可导入**的 Python 模块，要求：

1. **必须**定义一个继承自 `BaseSkill` 的类（从 `gibh_agent.skills.base_skill` 导入）。
2. 类属性（ClassVar）：
   - `skill_id`：与用户给出的 skill_id 完全一致（蛇形，仅字母数字下划线）。
   - `display_name`：中文展示名。
   - `description`：用于 LLM 路由的自然语言描述（含适用场景）。
   - `category`：学科大类（如「化学与分子信息学」）。
   - `sub_category`：子分类。
   - `aliases`：列表，含 tool_chain_key 及 2–5 个英文/中文别名。
   - `tool_chain_key`：磐石 ToolChain key（若有）。
   - `required_parameters`：人类可读的必填项说明列表（供意图路由 MISSING_PARAM）。
   - `__dependencies__`：**必须**在类级别声明（见下条）；无额外依赖时写 `[]`。
   - `__abstractskill__ = False`
2b. **依赖自报（`__dependencies__`）** — 类属性 `ClassVar[List[str]]` 或普通类属性均可：
   - 调用了第三方 **Python 库** → 写入 `pip:包名`（如 `pip:scanpy`、`pip:rdkit-pypi`）。
   - 通过 **subprocess** 调用底层 **Linux CLI** → 结合生信常识写入系统包，如 `apt:samtools`、`apt:openbabel`、`conda:bwa`、`conda:samtools`。
   - 仅标准库、无外部 CLI/包 → `__dependencies__ = []`。
   - 示例：`__dependencies__ = ["pip:rdkit-pypi", "apt:openbabel"]`
3. **路径类参数命名铁律**（仅允许下列名称，禁止自造 input_path / fasta_file 等）：
   `file_path`, `data_path`, `input_dir`, `matrix_dir`, `image_path`, `mask_path`, `sequence_or_path`
4. 实现 `execute(self, **kwargs) -> dict`：
   - 返回值**必须**为 `{"status": "success"|"error", "message": str, ...}`。
   - 重型计算须用 `asyncio.create_subprocess_exec` 或注释标明应部署为 Worker/TaaS，**禁止**在主进程做大规模 CPU/GPU 循环。
   - 对未知/占位实现：返回 status=error 并说明「待接入底层算子」，或返回 success 且 message 说明为脚手架。
5. **不要**在类上使用 `@registry.register`；BaseSkill 会在导入时自动注册。
6. **不要**手写 Pydantic 模型；参数类型通过 `execute` 的函数签名类型注解表达（ToolRegistry 会生成 schema）。
7. 文件顶部：`# -*- coding: utf-8 -*-` 与必要 imports（含 `from __future__ import annotations`）。
8. 只输出一个 ```python 代码块，不要解释文字。

参考骨架（须替换为真实逻辑）：

```python
# -*- coding: utf-8 -*-
from __future__ import annotations

import asyncio
import logging
from typing import Any, ClassVar, Dict, List

from gibh_agent.skills.base_skill import BaseSkill

logger = logging.getLogger(__name__)


class ExampleSkill(BaseSkill):
    __abstractskill__ = False
    skill_id: ClassVar[str] = "example_tool"
    display_name: ClassVar[str] = "示例工具"
    description: ClassVar[str] = "……"
    category: ClassVar[str] = "……"
    sub_category: ClassVar[str] = "……"
    aliases: ClassVar[List[str]] = ["example", "示例"]
    tool_chain_key: ClassVar[str] = "ExampleTool"
    required_parameters: ClassVar[List[str]] = ["输入文件 file_path 或序列 sequence_or_path"]
    __dependencies__: ClassVar[List[str]] = []  # 例: ["pip:rdkit-pypi", "apt:samtools"]

    def execute(self, file_path: str = "", sequence_or_path: str = "") -> Dict[str, Any]:
        # TODO: subprocess 调用 CLI 或 HTTP Worker
        return {"status": "error", "message": "底层算子待接入", "skill_id": self.skill_id}
```
"""


def generate_skill_code(
    client,
    skill_id: str,
    parsed: Dict[str, str],
    *,
    model: str,
) -> str:
    user_block = row_to_prompt_block(parsed, skill_id)
    messages = [
        {"role": "system", "content": SKILL_CODEGEN_SYSTEM_PROMPT},
        {
            "role": "user",
            "content": (
                "请为以下工具生成完整 Python 技能模块。\n\n"
                f"{user_block}\n"
                "注意：类名请用 PascalCase，且 skill_id 必须与上文一致。"
            ),
        },
    ]
    resp = client.chat(messages=messages, model=model, temperature=0.2, max_tokens=8192)
    content = (resp.choices[0].message.content or "").strip()
    code = extract_python_from_llm(content)
    if "BaseSkill" not in code or "skill_id" not in code:
        raise ValueError(f"LLM 输出未包含有效 BaseSkill 代码 (skill_id={skill_id})")
    if "__dependencies__" not in code:
        logger.warning("LLM 输出缺少 __dependencies__，将注入空列表占位 (skill_id=%s)", skill_id)
        code = code.replace(
            "__abstractskill__ = False",
            "__abstractskill__ = False\n    __dependencies__: ClassVar[List[str]] = []",
            1,
        )
    return code


def load_csv_rows(csv_path: Path) -> List[Dict[str, Any]]:
    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    df = df.fillna("")
    return df.to_dict(orient="records")


def main() -> int:
    parser = argparse.ArgumentParser(description="Skill Factory: CSV → BaseSkill Python 模块")
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "docs" / "生物信息分析工具.csv",
        help="工具清单 CSV",
    )
    parser.add_argument("--out-dir", type=Path, default=ROOT / "gibh_agent" / "skills")
    parser.add_argument("--model", default=DEFAULT_FACTORY_MODEL)
    parser.add_argument("--start", type=int, default=0, help="起始行（0-based）")
    parser.add_argument("--limit", type=int, default=0, help="最多生成条数，0=到文件末尾")
    parser.add_argument("--dry-run", action="store_true", help="只打印将生成的 skill_id，不调 LLM")
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", action="store_false", dest="skip_existing")
    parser.add_argument("--sleep", type=float, default=1.0, help="每条 LLM 调用间隔秒数")
    parser.add_argument("--manifest", type=Path, default=None, help="写入生成清单 JSON")
    parser.add_argument(
        "--deps-report",
        type=Path,
        default=ROOT / "docs" / "skill_missing_dependencies.txt",
        help="汇总 __dependencies__ 的报告路径",
    )
    parser.add_argument(
        "--deps-report-scope",
        choices=("all", "batch"),
        default="all",
        help="all=扫描 out-dir 下全部 skill_*.py；batch=仅本次 written 的文件",
    )
    parser.add_argument("--no-deps-report", action="store_true", help="跳过依赖汇总报告")
    args = parser.parse_args()

    if not args.csv.is_file():
        logger.error("CSV 不存在: %s", args.csv)
        return 1

    rows = load_csv_rows(args.csv)
    end = len(rows) if args.limit <= 0 else min(len(rows), args.start + args.limit)
    slice_rows = rows[args.start : end]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest: List[Dict[str, Any]] = []
    written_paths: List[Path] = []
    client = None if args.dry_run else factory_llm_client(args.model)

    for i, row in enumerate(slice_rows):
        parsed = parse_csv_row(row)
        if not parsed["display_name"]:
            continue
        sid = skill_id_from_row(row)
        out_name = output_filename_for_skill_id(sid)
        out_path = args.out_dir / out_name

        if args.skip_existing and out_path.is_file():
            logger.info("[%d] 跳过已存在: %s -> %s", args.start + i, sid, out_name)
            manifest.append({"skill_id": sid, "path": str(out_path), "status": "skipped"})
            continue

        if args.dry_run:
            logger.info("[%d] dry-run %s -> %s", args.start + i, sid, out_name)
            manifest.append({"skill_id": sid, "path": str(out_path), "status": "dry_run"})
            continue

        try:
            code = generate_skill_code(client, sid, parsed, model=args.model)
            out_path.write_text(code + "\n", encoding="utf-8")
            deps = extract_dependencies_from_source(code)
            logger.info("[%d] ✅ 已写入 %s (deps=%s)", args.start + i, out_path, deps or "[]")
            written_paths.append(out_path)
            manifest.append(
                {
                    "skill_id": sid,
                    "path": str(out_path),
                    "status": "written",
                    "dependencies": deps,
                }
            )
        except Exception as e:
            logger.error("[%d] ❌ %s: %s", args.start + i, sid, e)
            manifest.append({"skill_id": sid, "path": str(out_path), "status": "error", "error": str(e)})

        if args.sleep > 0:
            time.sleep(args.sleep)

    manifest_path = args.manifest or (args.out_dir / "_skill_factory_manifest.json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    logger.info("清单已写入: %s (%d 条)", manifest_path, len(manifest))

    if not args.no_deps_report:
        scope_files = written_paths if args.deps_report_scope == "batch" else None
        n_deps, n_files = write_skill_dependency_report(
            args.out_dir,
            args.deps_report,
            only_files=scope_files,
        )
        logger.info(
            "依赖报告已写入: %s（%d 个技能文件声明依赖，%d 条去重条目）",
            args.deps_report,
            n_files,
            n_deps,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
