# -*- coding: utf-8 -*-
"""
技能基类 — Skill Factory 产物须继承本类。

- 子类在模块导入时自动注册到 ToolRegistry（供 Skill Fast-Lane / LLM 填参）
- 元数据同步供 Master Intent Router 注入
"""
from __future__ import annotations

import inspect
import logging
from abc import ABC, abstractmethod
from typing import Any, ClassVar, Dict, List

from pydantic import BaseModel, Field

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution

logger = logging.getLogger(__name__)

# 架构宪法：路径类参数仅允许下列名称（见 .cursor/rules/gibh-architecture-constitution.mdc）
ALLOWED_PATH_PARAM_NAMES = frozenset(
    {
        "file_path",
        "data_path",
        "input_dir",
        "matrix_dir",
        "image_path",
        "mask_path",
        "sequence_or_path",
    }
)


class SkillMeta(BaseModel):
    """技能广场 / 意图路由用元数据（与 IntentRegistryEntry 字段对齐）。"""

    skill_id: str
    display_name: str
    description: str = ""
    category: str = "General"
    sub_category: str = ""
    aliases: List[str] = Field(default_factory=list)
    required_parameters: List[str] = Field(default_factory=list)
    tool_chain_key: str = ""
    dependencies: List[str] = Field(default_factory=list)


class BaseSkill(ABC):
    """
    可热挂载技能基类。

    子类须定义 ``skill_id``、``display_name``、``description``，并实现 ``execute``。
    模块被 skill_registry 导入后自动调用 ``_register_tool``。
    """

    __abstractskill__: ClassVar[bool] = True

    skill_id: ClassVar[str] = ""
    display_name: ClassVar[str] = ""
    description: ClassVar[str] = ""
    category: ClassVar[str] = "General"
    sub_category: ClassVar[str] = ""
    aliases: ClassVar[List[str]] = []
    required_parameters: ClassVar[List[str]] = []
    tool_chain_key: ClassVar[str] = ""
    output_type: ClassVar[str] = "json"
    # Skill Factory 依赖自报：pip:pkg / apt:deb / conda:pkg（供 docs/skill_missing_dependencies.txt 汇总）
    __dependencies__: ClassVar[List[str]] = []

    _registered_ids: ClassVar[set[str]] = set()

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        if getattr(cls, "__abstractskill__", False):
            return
        if not (cls.skill_id or "").strip():
            logger.warning("跳过未定义 skill_id 的 BaseSkill 子类: %s", cls.__name__)
            return
        cls._register_tool()

    @classmethod
    def meta(cls) -> SkillMeta:
        return SkillMeta(
            skill_id=(cls.skill_id or "").strip(),
            display_name=(cls.display_name or cls.skill_id or "").strip(),
            description=(cls.description or "").strip(),
            category=(cls.category or "General").strip(),
            sub_category=(cls.sub_category or "").strip(),
            aliases=list(cls.aliases or []),
            required_parameters=list(cls.required_parameters or []),
            tool_chain_key=(cls.tool_chain_key or "").strip(),
            dependencies=list(cls.__dependencies__ or []),
        )

    @abstractmethod
    def execute(self, **kwargs: Any) -> Dict[str, Any]:
        """业务逻辑；须返回含 status/message 的字典。"""

    @classmethod
    def _register_tool(cls) -> None:
        sid = (cls.skill_id or "").strip()
        if not sid:
            return
        if sid in cls._registered_ids:
            logger.debug("技能 %s 已注册，跳过重复", sid)
            return

        desc = (cls.description or cls.display_name or sid).strip()
        cat = (cls.category or "General").strip()
        otype = (cls.output_type or "json").strip()

        if inspect.iscoroutinefunction(cls.execute):

            @registry.register(name=sid, description=desc, category=cat, output_type=otype)
            @safe_tool_execution
            async def _async_skill_entrypoint(**kwargs: Any) -> Dict[str, Any]:
                return await cls().execute(**kwargs)

        else:

            @registry.register(name=sid, description=desc, category=cat, output_type=otype)
            @safe_tool_execution
            def _skill_entrypoint(**kwargs: Any) -> Dict[str, Any]:
                from gibh_agent.skills.launch_skill_isolated import (
                    delegate_launch_skill_execute,
                    should_delegate_to_launch_worker,
                )

                if should_delegate_to_launch_worker(sid):
                    return delegate_launch_skill_execute(sid, kwargs)
                return cls().execute(**kwargs)

        cls._registered_ids.add(sid)
        logger.info("✅ BaseSkill 已注册工具: %s (%s)", sid, cls.__name__)

    @classmethod
    def validate_path_param_names(cls, param_names: List[str]) -> List[str]:
        """生成器/审查用：返回不在白名单中的路径参数名。"""
        bad = [n for n in param_names if n.endswith("_path") or n.endswith("_dir")]
        return [n for n in bad if n not in ALLOWED_PATH_PARAM_NAMES]


def is_concrete_skill_class(obj: Any) -> bool:
    return (
        inspect.isclass(obj)
        and issubclass(obj, BaseSkill)
        and obj is not BaseSkill
        and not getattr(obj, "__abstractskill__", False)
    )
