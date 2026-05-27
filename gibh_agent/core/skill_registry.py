# -*- coding: utf-8 -*-
"""
技能动态发现与热挂载 — Skill Factory 注册中心。

启动时（及显式 reload）扫描 ``gibh_agent/skills/`` 下所有 ``.py``，
导入继承自 ``BaseSkill`` 的类，自动注入 ToolRegistry 与 Intent Router。
"""
from __future__ import annotations

import importlib
import importlib.util
import logging
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Dict, List, Optional, Set, Type

from gibh_agent.skills.base_skill import BaseSkill, SkillMeta, is_concrete_skill_class

logger = logging.getLogger(__name__)

_SKILLS_PACKAGE = "gibh_agent.skills"
_SKIP_MODULES = frozenset({"__init__", "base_skill", "_generated_readme"})

# 已加载模块名 → 上次 mtime（热重载）
_loaded_modules: Dict[str, float] = {}
_discovered_metas: List[SkillMeta] = []


def get_skills_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "skills"


def _iter_skill_py_files(skills_dir: Optional[Path] = None) -> List[Path]:
    root = skills_dir or get_skills_dir()
    if not root.is_dir():
        return []
    out: List[Path] = []
    for p in sorted(root.glob("*.py")):
        if p.stem in _SKIP_MODULES or p.name.startswith("_"):
            continue
        out.append(p)
    return out


def _module_name_for_path(py_path: Path) -> str:
    return f"{_SKILLS_PACKAGE}.{py_path.stem}"


def _collect_skill_classes(module: ModuleType) -> List[Type[BaseSkill]]:
    found: List[Type[BaseSkill]] = []
    for attr_name in dir(module):
        if attr_name.startswith("_"):
            continue
        obj = getattr(module, attr_name, None)
        if is_concrete_skill_class(obj):
            found.append(obj)
    return found


def _import_skill_module(py_path: Path, *, force_reload: bool = False) -> Optional[ModuleType]:
    mod_name = _module_name_for_path(py_path)
    mtime = py_path.stat().st_mtime
    if not force_reload and mod_name in sys.modules and _loaded_modules.get(mod_name) == mtime:
        return sys.modules[mod_name]

    if mod_name in sys.modules:
        try:
            importlib.reload(sys.modules[mod_name])
        except Exception as e:
            logger.warning("热重载技能模块失败 %s: %s", mod_name, e)
            return None
        _loaded_modules[mod_name] = mtime
        return sys.modules[mod_name]

    spec = importlib.util.spec_from_file_location(mod_name, py_path)
    if spec is None or spec.loader is None:
        logger.error("无法为技能文件创建 spec: %s", py_path)
        return None
    module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        logger.error("加载技能模块失败 %s: %s", py_path.name, e, exc_info=True)
        sys.modules.pop(mod_name, None)
        return None
    _loaded_modules[mod_name] = mtime
    return module


def discover_and_register_skills(
    *,
    skills_dir: Optional[Path] = None,
    force_reload: bool = False,
) -> Dict[str, int]:
    """
    扫描技能目录、导入 BaseSkill 子类（导入即注册 ToolRegistry）。

    Returns:
        {"loaded": n, "failed": m, "skills": k}
    """
    global _discovered_metas
    root = skills_dir or get_skills_dir()
    root.mkdir(parents=True, exist_ok=True)

    loaded = 0
    failed = 0
    metas: List[SkillMeta] = []
    seen_ids: Set[str] = set()

    for py_path in _iter_skill_py_files(root):
        module = _import_skill_module(py_path, force_reload=force_reload)
        if module is None:
            failed += 1
            continue
        classes = _collect_skill_classes(module)
        if not classes:
            logger.debug("技能文件无 BaseSkill 子类: %s", py_path.name)
            loaded += 1
            continue
        for cls in classes:
            meta = cls.meta()
            if meta.skill_id in seen_ids:
                logger.warning("重复 skill_id 已跳过: %s (%s)", meta.skill_id, py_path.name)
                continue
            seen_ids.add(meta.skill_id)
            metas.append(meta)
        loaded += 1
        logger.info("已发现技能模块: %s (%d 类)", py_path.name, len(classes))

    _discovered_metas = metas
    return {"loaded": loaded, "failed": failed, "skills": len(metas)}


def reload_skills() -> Dict[str, int]:
    """热重载：强制重新导入 skills 目录下所有模块。"""
    return discover_and_register_skills(force_reload=True)


def get_discovered_skill_metas() -> List[SkillMeta]:
    return list(_discovered_metas)


def list_registered_skill_ids() -> List[str]:
    return [m.skill_id for m in _discovered_metas]
