"""
DNA 序列续写生成：调用外部 HTTP `/generate` 服务（与其它 JSON POST 预测类工具同接入范式）。

环境变量 `DNA_GENERATE_API_URL` 可覆盖默认端点。
"""
from __future__ import annotations

import asyncio
import logging
import os
from typing import Any, Dict, List, Optional

import httpx

from gibh_agent.core.tool_registry import registry

logger = logging.getLogger(__name__)

DEFAULT_DNA_GENERATE_URL = "http://8.130.86.61:38008/generate"
_RETRYABLE_HTTP = frozenset({500, 502, 503, 504})
_MAX_RETRIES = 3
_RETRY_DELAY_SEC = 3.0
_MAX_PROBS_IN_RESPONSE = 32
# 读响应可长；连接阶段单独限时，避免不可达主机每次卡满 300s
_CONNECT_TIMEOUT_SEC = 20.0
_READ_TIMEOUT_SEC = 300.0


def _humanize_dna_api_request_error(url: str, exc: Exception) -> str:
    """远程 DNA 生成 API 不可达时给用户可操作的说明（非 LLM 问题）。"""
    detail = str(exc).strip() or type(exc).__name__
    return (
        "无法连接到 DNA 生成服务（对端节点不可达、本机防火墙/安全组未放行、或服务过载）。\n"
        f"当前请求端点：`{url}`\n\n"
        "**建议排查：**\n"
        "1. 在**部署 Omics Agent 的服务器**上测试能否访问上述 URL（curl/浏览器）。\n"
        "2. 若服务提供方文档中的基地址不同，请设置环境变量 `DNA_GENERATE_API_URL` 为可达地址。\n"
        "3. 与基础设施/服务运维确认该端口是否对您的出口 IP 开放。\n\n"
        f"底层错误：`{detail}`"
    )


def _clamp_num_tokens(n: int) -> int:
    if n < 1:
        return 1
    if n > 200:
        return 200
    return n


def _clamp_top_k(k: int) -> int:
    if k < 0:
        return 0
    if k > 10:
        return 10
    return k


def _clamp_temperature(t: float) -> float:
    if t < 0.0:
        return 0.0
    if t > 1.0:
        return 1.0
    return t


def _clamp_top_p(p: float) -> float:
    if p < 0.0:
        return 0.0
    if p > 1.0:
        return 1.0
    return p


@registry.register(
    name="dna_sequence_generate",
    description=(
        "根据给定 DNA 提示序列，调用外部 `/generate` 接口续写 num_tokens 长度。"
        "输入为文本 DNA（IUPAC 大写）；建议提示序列≤500 nt，可选分类学前缀 D__…;P__…;S__… 再接序列。"
        "可调 num_tokens(1–200)、temperature、top_k、top_p、enable_logits、enable_sampled_probs。"
    ),
    category="Genomics",
    output_type="json",
)
async def dna_sequence_generate(
    sequence: str,
    num_tokens: int = 100,
    temperature: float = 0.7,
    top_k: int = 3,
    top_p: float = 1.0,
    enable_logits: bool = False,
    enable_sampled_probs: bool = True,
) -> Dict[str, Any]:
    """
    POST JSON 至 `/generate`，返回生成序列与耗时、可选采样概率摘要。

    参数与官方 API 文档对齐；超出范围的标量会自动钳制到文档允许区间。
    """
    seq = (sequence or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if not seq:
        return {"status": "error", "message": "sequence 为空，请提供 DNA 提示序列。"}
    if len(seq) > 500:
        logger.warning(
            "DNA 提示序列长度=%s，超过服务文档建议的 500 nt，API 可能拒绝或截断。",
            len(seq),
        )

    url = (os.getenv("DNA_GENERATE_API_URL") or "").strip() or DEFAULT_DNA_GENERATE_URL
    nt = _clamp_num_tokens(int(num_tokens))
    payload: Dict[str, Any] = {
        "sequence": seq,
        "num_tokens": nt,
        "temperature": _clamp_temperature(float(temperature)),
        "top_k": _clamp_top_k(int(top_k)),
        "top_p": _clamp_top_p(float(top_p)),
        "enable_logits": bool(enable_logits),
        "enable_sampled_probs": bool(enable_sampled_probs),
    }

    timeout = httpx.Timeout(_READ_TIMEOUT_SEC, connect=_CONNECT_TIMEOUT_SEC)
    async with httpx.AsyncClient(timeout=timeout) as client:
        for attempt in range(1, _MAX_RETRIES + 1):
            try:
                response = await client.post(url, json=payload, headers={"Content-Type": "application/json"})
            except httpx.TimeoutException as e:
                logger.error("DNA generate 超时 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e)
                if attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": (
                        f"DNA 生成请求超时（已重试 {_MAX_RETRIES} 次）。"
                        "若对端服务繁忙或生成过长，可稍后再试或减小 num_tokens。\n"
                        f"详情: {e}"
                    ),
                }
            except httpx.RequestError as e:
                logger.error("DNA generate 网络错误 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e, exc_info=True)
                if attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": _humanize_dna_api_request_error(url, e),
                }

            raw_text = response.text or ""
            logger.info(
                "DNA generate API | attempt=%s | status=%s | body_len=%s",
                attempt,
                response.status_code,
                len(raw_text),
            )

            if response.status_code != 200:
                if response.status_code in _RETRYABLE_HTTP and attempt < _MAX_RETRIES:
                    logger.warning("DNA generate HTTP %s，%s 秒后重试", response.status_code, _RETRY_DELAY_SEC)
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": f"HTTP {response.status_code}，详情: {raw_text[:2000]}",
                }

            try:
                result: Dict[str, Any] = response.json()
            except Exception as e:
                return {
                    "status": "error",
                    "message": f"响应非合法 JSON: {e}，原文: {raw_text[:1500]}",
                }

            out_seq = result.get("sequence")
            if not isinstance(out_seq, str) or not out_seq.strip():
                return {
                    "status": "error",
                    "message": f"响应缺少有效 sequence 字段: {result!r}"[:4000],
                }

            elapsed = result.get("elapsed_ms")
            gen_params = result.get("generation_params")
            probs: Optional[List[Any]] = None
            if isinstance(result.get("sampled_probs"), list):
                probs = result["sampled_probs"]

            probs_preview: Optional[List[Any]] = None
            probs_len = len(probs) if probs is not None else 0
            if probs is not None and probs_len > 0:
                probs_preview = probs[:_MAX_PROBS_IN_RESPONSE]

            return {
                "status": "success",
                "message": f"已生成 {len(out_seq)} nt（请求续写约 {nt} token），耗时约 {elapsed} ms。",
                "generated_sequence": out_seq.strip(),
                "elapsed_ms": elapsed,
                "generation_params": gen_params if isinstance(gen_params, dict) else {},
                "sampled_probs_count": probs_len,
                "sampled_probs_preview": probs_preview,
                "logits": result.get("logits"),
            }

    return {"status": "error", "message": "DNA 生成失败（未预期退出重试循环）"}
