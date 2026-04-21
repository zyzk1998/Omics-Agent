# -*- coding: utf-8 -*-
"""
云端 LLM — 模型注册中心（唯一权威路由表）。

禁止按前缀/斜杠「猜测」厂商：仅 MODEL_ROUTING_TABLE 中登记的 model id 可解析，
确保 base_url、api_key 与 model 三者绝对一致。

官方文档：
- 智谱 GLM：https://docs.bigmodel.cn/cn/guide/start/introduction
- DeepSeek：https://api-docs.deepseek.com/zh-cn/
- Kimi：https://platform.kimi.com/docs/overview
"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# 未设置或填写已弃用厂商名时，默认走 DeepSeek 官方（硅基流动 SiliconFlow 已从本仓库运维路径移除）
_DEFAULT_CLOUD_PROVIDER_ID = "deepseek"

# 厂商级网关（不含具体 model id）
# 注：硅基流动（SiliconFlow）聚合网关已弃用，不再在 PROVIDERS 中维护；历史 .env 中 LLM_CLOUD_PROVIDER=siliconflow 将在 normalize 中回落到 deepseek。
PROVIDERS: Dict[str, Dict[str, str]] = {
    "deepseek": {
        "base_url": "https://api.deepseek.com/v1",
        "api_key_env": "DEEPSEEK_API_KEY",
        "model_env": "DEEPSEEK_MODEL",
        "default_model": "deepseek-reasoner",
        "label": "DeepSeek",
    },
    "glm": {
        "base_url": "https://open.bigmodel.cn/api/paas/v4/",
        "api_key_env": "ZHIPU_API_KEY",
        "model_env": "GLM_MODEL",
        "default_model": "glm-5.1",
        "label": "ZhipuGLM",
    },
    "kimi": {
        "base_url": "https://api.moonshot.cn/v1",
        "api_key_env": "MOONSHOT_API_KEY",
        "model_env": "KIMI_MODEL",
        "default_model": "kimi-k2.5",
        "label": "MoonshotKimi",
    },
}

# ---------------------------------------------------------------------------
# 手术点 1：绝对权威的模型 → 厂商 → 密钥环境变量（显式登记，拒绝猜测）
# 前端 services/nginx/html/index.html 下拉 value 须与本表 key 完全一致（如 deepseek-reasoner），
# 不要使用硅基目录名如 deepseek-ai/DeepSeek-R1（已移除）；否则 validate_and_resolve_model_name 抛 unsupported_model。
# ---------------------------------------------------------------------------
MODEL_ROUTING_TABLE: Dict[str, Dict[str, str]] = {
    # DeepSeek 官方
    "deepseek-reasoner": {"provider": "deepseek", "env_key": "DEEPSEEK_API_KEY"},
    "deepseek-chat": {"provider": "deepseek", "env_key": "DEEPSEEK_API_KEY"},
    # 智谱 GLM 官方
    "glm-5.1": {"provider": "glm", "env_key": "ZHIPU_API_KEY"},
    # Kimi / Moonshot 官方
    "kimi-k2.5": {"provider": "kimi", "env_key": "MOONSHOT_API_KEY"},
    # 以下曾为硅基流动目录型 id，已随 SiliconFlow 弃用从注册表移除；请改用左侧官方直连 id。
}


def list_registered_models() -> List[str]:
    return sorted(MODEL_ROUTING_TABLE.keys())


def chat_mode_skip_openai_tools(model_id: Optional[str]) -> bool:
    """
    闲聊（非 DeepReAct）下，部分推理模型与「首轮 tools + MCP」组合在官方网关侧易 400/不兼容，
    故对这些 model id 跳过 tools，仅走纯文本流式（需联网时请改用 deepseek-chat 等）。
    """
    mid = (model_id or "").strip()
    return mid == "deepseek-reasoner"


def normalize_provider_id(raw: Optional[str]) -> str:
    if not raw:
        return _DEFAULT_CLOUD_PROVIDER_ID
    s = str(raw).strip().lower()
    aliases = {
        "zhipu": "glm",
        "bigmodel": "glm",
        "moonshot": "kimi",
        # 历史别名：硅基已弃用，映射到默认官方直连
        "sf": _DEFAULT_CLOUD_PROVIDER_ID,
        "siliconflow": _DEFAULT_CLOUD_PROVIDER_ID,
    }
    s = aliases.get(s, s)
    return s if s in PROVIDERS else _DEFAULT_CLOUD_PROVIDER_ID


def get_active_provider_id() -> str:
    return normalize_provider_id(os.getenv("LLM_CLOUD_PROVIDER"))


def get_default_chat_model() -> str:
    """
    未传 model_name 时的默认 id：优先 DEFAULT_CHAT_MODEL（须在注册表中），
    否则使用当前 LLM_CLOUD_PROVIDER 对应厂商的 default_model（亦须在注册表中）。
    """
    direct = (os.getenv("DEFAULT_CHAT_MODEL") or "").strip()
    if direct:
        if direct not in MODEL_ROUTING_TABLE:
            raise ValueError(
                f"环境变量 DEFAULT_CHAT_MODEL={direct!r} 未在 MODEL_ROUTING_TABLE 中登记。"
                f"已登记: {', '.join(list_registered_models())}"
            )
        return direct
    pid = get_active_provider_id()
    mid = (os.getenv(PROVIDERS[pid]["model_env"]) or "").strip() or PROVIDERS[pid]["default_model"]
    if mid not in MODEL_ROUTING_TABLE:
        raise ValueError(
            f"当前 LLM_CLOUD_PROVIDER={pid} 的默认模型 {mid!r} 未在 MODEL_ROUTING_TABLE 中登记。"
            f"请设置 DEFAULT_CHAT_MODEL 为已登记 id，或扩展 MODEL_ROUTING_TABLE / 修正 {PROVIDERS[pid]['model_env']}。"
        )
    return mid


def validate_and_resolve_model_name(raw: Optional[str]) -> str:
    """
    编排器入口：将请求中的 model_name 规范为注册表中的绝对 id；非法则抛错（由上层转 SSE）。
    raw 为空则使用 get_default_chat_model()。
    """
    s = (raw or "").strip()
    if not s:
        return get_default_chat_model()
    if s not in MODEL_ROUTING_TABLE:
        raise ValueError(
            f"不支持的模型: {s!r}。已登记模型: {', '.join(list_registered_models())}"
        )
    return s


def get_client_config(model_id: str, api_key_override: Optional[str] = None) -> Dict[str, Any]:
    """
    根据注册表中的 model id 返回本次 HTTP 调用所需的完整配置（权威、无猜测）。

    Args:
        model_id: MODEL_ROUTING_TABLE 中登记的模型 id。
        api_key_override: 非空时优先使用该密钥（兼容脚本/测试显式传入）；否则读 env_key 对应环境变量。

    Returns:
        dict: base_url, api_key, model（与请求体一致）, provider（日志标签）, provider_id, env_key
    """
    mid = (model_id or "").strip()
    if mid not in MODEL_ROUTING_TABLE:
        raise ValueError(
            f"不支持的模型: {mid!r}。已登记: {', '.join(list_registered_models())}"
        )
    row = MODEL_ROUTING_TABLE[mid]
    pid = row["provider"]
    env_key = row["env_key"]
    if pid not in PROVIDERS:
        raise ValueError(f"内部错误: provider {pid!r} 未在 PROVIDERS 中定义")
    meta = PROVIDERS[pid]
    if api_key_override is not None and str(api_key_override).strip():
        api_key = str(api_key_override).strip()
    else:
        api_key = (os.getenv(env_key) or "").strip()
    if not api_key:
        raise ValueError(
            f"模型 {mid!r} 需要环境变量 {env_key}（厂商 {meta['label']}）。请在 .env 中配置。"
        )
    base_url = meta["base_url"]
    logger.info(
        "🚀 [ModelRegistry] base_url=%s model=%s provider=%s env_key=%s",
        base_url,
        mid,
        meta["label"],
        env_key,
    )
    return {
        "base_url": base_url,
        "api_key": api_key,
        "model": mid,
        "provider": meta["label"],
        "provider_id": pid,
        "env_key": env_key,
    }


def assert_default_model_configurable() -> None:
    """启动自检：默认模型必须可解析（避免进程起来后首聊才爆）。"""
    get_client_config(get_default_chat_model())


def build_cloud_llm_client_dict() -> Dict[str, Any]:
    """
    兼容旧调用：等价于「当前默认模型」的 create_from_config 字典。
    新架构下领域智能体不应依赖此方法；仅保留供极少数过渡路径使用。
    """
    cfg = get_client_config(get_default_chat_model())
    return {
        "base_url": cfg["base_url"],
        "api_key": cfg["api_key"],
        "model": cfg["model"],
        "temperature": 0.1,
        "max_tokens": 4096,
        "timeout": None,
        "provider": cfg["provider"],
        "_provider_id": cfg["provider_id"],
        "_provider_label": cfg["provider"],
    }


def resolve_provider_tuple(model_name: str) -> Tuple[str, str, str]:
    """供 llm_client._resolve_provider_for_model 使用。"""
    c = get_client_config(model_name)
    return (c["base_url"], c["api_key"], c["provider"])
