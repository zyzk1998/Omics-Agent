"""
ZIP / Markdown 解析：安全解压到 dynamic_skills 根目录下唯一子目录，并生成 OpenAI Function 风格 parameters schema。
"""
from __future__ import annotations

import json
import logging
import os
import re
import shutil
import uuid
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple

SkillType = Literal["prompt", "script"]

import yaml

logger = logging.getLogger(__name__)

DEFAULT_DYNAMIC_ROOT = "/app/uploads/dynamic_skills"


def _dynamic_root() -> Path:
    root = Path(os.path.abspath(os.environ.get("DYNAMIC_SKILLS_DIR", DEFAULT_DYNAMIC_ROOT)))
    root.mkdir(parents=True, exist_ok=True)
    return root


def get_dynamic_skills_root() -> Path:
    """对外：动态技能解压根目录（绝对路径）。"""
    return _dynamic_root()


def _is_under_parent(child: Path, parent: Path) -> bool:
    try:
        child = child.resolve()
        parent = parent.resolve()
        return os.path.commonpath([str(child), str(parent)]) == str(parent)
    except (OSError, ValueError):
        return False


def _safe_extract_zip(zip_path: Path, dest_dir: Path) -> None:
    """解压 zip，拒绝路径穿越（Zip Slip）。"""
    dest_dir = dest_dir.resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        for member in zf.infolist():
            if member.is_dir():
                continue
            raw = member.filename
            if raw.startswith("/") or ".." in Path(raw).parts:
                raise ValueError(f"非法 zip 成员路径: {raw!r}")
            target = (dest_dir / raw).resolve()
            if not _is_under_parent(target, dest_dir):
                raise ValueError(f"路径穿越被拒绝: {raw!r} -> {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(member, "r") as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)


def _section(md: str, heading: str) -> str:
    """取 ## Heading 与下一同级 ## 之间的正文。"""
    pat = re.compile(
        rf"(?ms)^##\s+{re.escape(heading)}\s*$\n(.*?)(?=^##\s+|\Z)",
    )
    m = pat.search(md)
    return (m.group(1).strip() if m else "") or ""


def _slug_name(name: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_]+", "_", (name or "").strip().lower()).strip("_")
    return s or "dynamic_skill"


def parse_skill_md(md_content: str) -> Dict[str, Any]:
    """
    从 Markdown 提取 Name / Description / Parameters，生成 OpenAI tools.function.parameters 风格 JSON Schema。
    Parameters 区支持行：`- param_name (string): 说明` 或 `- param_name: 说明`（默认 string）。
    """
    text = md_content or ""
    name_block = _section(text, "Name") or _section(text, "技能名称")
    if not name_block:
        m = re.search(r"^\s*#\s+(.+)$", text, re.MULTILINE)
        name_block = m.group(1).strip() if m else "unnamed_skill"
    else:
        name_block = name_block.splitlines()[0].strip()

    desc = _section(text, "Description") or _section(text, "描述") or ""
    if not desc:
        desc = (text[:500] + "…") if len(text) > 500 else text

    params_body = _section(text, "Parameters") or _section(text, "参数")
    properties: Dict[str, Any] = {}
    required: List[str] = []

    line_re = re.compile(
        r"^\s*[-*]\s*`?([a-zA-Z_][a-zA-Z0-9_]*)`?\s*(?:\(\s*([a-zA-Z]+)\s*\))?\s*:?\s*(.*)$",
    )
    for line in params_body.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = line_re.match(line)
        if not m:
            continue
        pname, ptype, pdesc = m.group(1), (m.group(2) or "string").lower(), (m.group(3) or "").strip()
        json_type = "string"
        if ptype in ("int", "integer"):
            json_type = "integer"
        elif ptype in ("float", "number"):
            json_type = "number"
        elif ptype in ("bool", "boolean"):
            json_type = "boolean"
        elif ptype in ("object", "dict", "array", "list"):
            json_type = "object" if ptype in ("object", "dict") else "array"
        properties[pname] = {"type": json_type, "description": pdesc or pname}
        if "(required)" in line.lower() or "**required**" in line.lower():
            required.append(pname)

    parameters_schema: Dict[str, Any] = {
        "type": "object",
        "properties": properties,
        "required": required,
    }

    return {
        "name": _slug_name(name_block),
        "display_name": name_block[:255],
        "description": desc.strip()[:8000],
        "parameters_schema": parameters_schema,
    }


def parse_skill_markdown_upload(md_content: str, dest_dir: Path, skill_uid: str) -> Dict[str, Any]:
    """
    纯 .md 上架：仅写入磁盘元数据，**不**生成 main.py。skill_type 恒为 prompt。
    """
    dest_dir = dest_dir.resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    skill_md = dest_dir / "SKILL.md"
    skill_md.write_text(md_content or "", encoding="utf-8")
    meta = parse_skill_md(md_content)
    return {
        **meta,
        "skill_type": "prompt",
        "main_py_path": None,
        "extract_dir": str(dest_dir),
        "skill_uid": skill_uid,
        "raw_yaml": None,
    }


def parse_skill_zip(zip_path: str) -> Dict[str, Any]:
    """
    解压到 dynamic_skills/{uuid}/。
    - 含唯一 main.py → skill_type=script，可经 worker 执行。
    - 不含 main.py 但含 SKILL.md/skill.md/skill.yaml → skill_type=prompt，仅提示词技能。
    """
    zp = Path(os.path.abspath(zip_path))
    if not zp.is_file():
        raise FileNotFoundError(f"zip 不存在: {zp}")

    root = _dynamic_root()
    skill_uid = str(uuid.uuid4())
    dest = (root / skill_uid).resolve()
    if not _is_under_parent(dest, root.resolve()):
        raise ValueError("目标目录非法")
    dest.mkdir(parents=True, exist_ok=True)

    _safe_extract_zip(zp, dest)

    found_mains = [p.resolve() for p in dest.rglob("main.py") if p.is_file()]
    skill_type: SkillType = "script"
    main_py: Optional[Path] = None
    if len(found_mains) == 1:
        main_py = found_mains[0]
        if not _is_under_parent(main_py, dest):
            raise ValueError("main.py 路径非法")
    elif len(found_mains) > 1:
        raise FileNotFoundError("压缩包内只能包含一个 main.py")
    else:
        skill_type = "prompt"

    skill_md = dest / "SKILL.md"
    skill_md_alt = dest / "skill.md"
    yaml_path = dest / "skill.yaml"

    result: Dict[str, Any] = {
        "extract_dir": str(dest),
        "main_py_path": str(main_py) if main_py else None,
        "skill_uid": skill_uid,
        "skill_type": skill_type,
    }

    if skill_md.is_file():
        raw = skill_md.read_text(encoding="utf-8", errors="replace")
        parsed = parse_skill_md(raw)
        result.update(parsed)
        result["skill_type"] = skill_type
        result["raw_yaml"] = None
    elif skill_md_alt.is_file():
        raw = skill_md_alt.read_text(encoding="utf-8", errors="replace")
        parsed = parse_skill_md(raw)
        result.update(parsed)
        result["skill_type"] = skill_type
        result["raw_yaml"] = None
    elif yaml_path.is_file():
        raw_yaml = yaml.safe_load(yaml_path.read_text(encoding="utf-8", errors="replace"))
        if not isinstance(raw_yaml, dict):
            raise ValueError("skill.yaml 必须是 YAML mapping")
        name = str(raw_yaml.get("name") or "unnamed_skill")
        desc = str(raw_yaml.get("description") or "")
        params = raw_yaml.get("parameters") or raw_yaml.get("parameters_schema")
        if isinstance(params, dict) and "type" in params:
            parameters_schema = params
        else:
            properties = {}
            req: List[str] = []
            if isinstance(params, list):
                for item in params:
                    if not isinstance(item, dict):
                        continue
                    pn = item.get("name")
                    if not pn:
                        continue
                    properties[str(pn)] = {
                        "type": item.get("type", "string"),
                        "description": item.get("description", ""),
                    }
                    if item.get("required"):
                        req.append(str(pn))
            parameters_schema = {
                "type": "object",
                "properties": properties,
                "required": req,
            }
        result.update(
            {
                "name": _slug_name(name),
                "display_name": name[:255],
                "description": desc[:8000],
                "parameters_schema": parameters_schema,
                "raw_yaml": raw_yaml,
                "skill_type": skill_type,
            }
        )
    else:
        raise FileNotFoundError("ZIP 内须包含 SKILL.md、skill.md 或 skill.yaml 之一")

    return result
