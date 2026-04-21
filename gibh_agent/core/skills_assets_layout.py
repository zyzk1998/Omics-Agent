# -*- coding: utf-8 -*-
"""
Skills 底层资产目录布局（服务器磁盘约定）。

三类并列目录（均在 SKILLS_ASSETS_ROOT 下）：
- prompt/     纯 Prompt：SKILL.md、说明文档、上架前校验材料等
- script/     脚本类：inbound（指挥上传 zip 暂存）、runtime（动态解压与 worker 沙箱根）
- models/     权重与大文件（不进 Git，由运维按需挂载）
"""
from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "get_skills_assets_root",
    "get_prompt_skills_dir",
    "get_script_skills_inbound_dir",
    "get_script_skills_runtime_dir",
    "get_models_skills_dir",
    "get_dynamic_skills_dir",
    "get_dynamic_skills_root",
]


def _mkdir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_skills_assets_root() -> Path:
    """
    技能资产根目录。
    默认：${UPLOAD_DIR}/skills_assets（容器内常为 /app/uploads/skills_assets，
    对应 compose 宿主机 ./data/uploads/skills_assets）。
    可整体覆盖：环境变量 SKILLS_ASSETS_ROOT。
    """
    explicit = os.environ.get("SKILLS_ASSETS_ROOT")
    if explicit:
        return _mkdir(Path(os.path.abspath(explicit.strip())))
    upload = os.environ.get("UPLOAD_DIR", "/app/uploads").strip()
    root = Path(os.path.abspath(upload)) / "skills_assets"
    return _mkdir(root)


def get_prompt_skills_dir() -> Path:
    """纯 Prompt / 文档类材料。"""
    return _mkdir(get_skills_assets_root() / "prompt")


def get_script_skills_inbound_dir() -> Path:
    """脚本类压缩包上传暂存（入库前人工/脚本校验）。"""
    return _mkdir(get_skills_assets_root() / "script" / "inbound")


def get_script_skills_runtime_dir() -> Path:
    """
    脚本类运行时根：ZIP 安全解压目录（uuid 子目录）、main.py 所在沙箱根。
    与 DYNAMIC_SKILLS_DIR 默认值一致（除非单独覆盖 DYNAMIC_SKILLS_DIR）。
    """
    return _mkdir(get_skills_assets_root() / "script" / "runtime")


def get_models_skills_dir() -> Path:
    """模型权重与大体量二进制（路径引用由工具参数或配置给出）。"""
    return _mkdir(get_skills_assets_root() / "models")


def get_dynamic_skills_dir() -> Path:
    """
    动态技能解压与 Worker 校验共用的根目录。
    若设置 DYNAMIC_SKILLS_DIR，优先使用（兼容旧部署）；
    否则为 script/runtime。
    """
    override = os.environ.get("DYNAMIC_SKILLS_DIR")
    if override and str(override).strip():
        return _mkdir(Path(os.path.abspath(override.strip())))
    return get_script_skills_runtime_dir()


# plugin_system/parser 等场所沿用的对外名称
def get_dynamic_skills_root() -> Path:
    return get_dynamic_skills_dir()
