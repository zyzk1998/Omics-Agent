"""
蛋白质序列层工具：B 细胞表位预测（BepiPred3）等。

与 `proteomics_tools.py`（质谱矩阵统计）互补；本模块聚焦 FASTA/序列级任务。
"""
from __future__ import annotations

import asyncio
import logging
import os
from typing import Any, Dict, Literal, Optional, Tuple

import httpx

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)

TopEpitopePercentageCutoff = Literal["top_20", "top_50", "top_70", "all"]

DEFAULT_BEPIPRED3_PREDICT_URL = "http://120.220.102.26:38089/predict"

_RETRYABLE_HTTP = frozenset({500, 502, 503, 504})
_MAX_RETRIES = 3
_RETRY_DELAY_SEC = 5.0

_GPU_BUSY_USER_MESSAGE = (
    "磐石服务器当前计算资源繁忙（可能由于排队人数过多或 GPU 满载）。请稍后重试。"
)


def _panshi_error_suggests_gpu_transient(error_message: str) -> bool:
    """磐石 data.error_message 是否像 GPU/CUDA 瞬时繁忙，可重试。"""
    if not error_message:
        return False
    low = error_message.lower()
    return "cuda" in low or "busy" in low or "unavailable" in low


def _resolve_fasta_content(sequence_or_path: str) -> str:
    """路径存在则读文件，否则视为 FASTA 纯文本。"""
    raw = (sequence_or_path or "").strip()
    if not raw:
        raise ValueError("sequence_or_path 为空")
    if os.path.isfile(raw):
        with open(raw, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    return raw


def _extract_output_urls(data: Dict[str, Any]) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """从磐石 `data.output_files` 解析 html / csv / zip 的 download_url。"""
    if not isinstance(data, dict):
        return None, None, None
    of = data.get("output_files")
    html_url: Optional[str] = None
    csv_url: Optional[str] = None
    zip_url: Optional[str] = None

    def _pick_url(entry: Any) -> Optional[str]:
        if isinstance(entry, dict):
            return entry.get("download_url") or entry.get("url")
        if isinstance(entry, str):
            return entry
        return None

    if isinstance(of, dict):
        html_url = _pick_url(of.get("html"))
        csv_url = _pick_url(of.get("csv"))
        zip_url = _pick_url(of.get("zip"))
    elif isinstance(of, list):
        for item in of:
            if not isinstance(item, dict):
                continue
            u = item.get("download_url") or item.get("url")
            if not u:
                continue
            name = str(item.get("name") or item.get("type") or item.get("format") or "").lower()
            if "html" in name or name.endswith(".html"):
                html_url = html_url or u
            elif "csv" in name or name.endswith(".csv"):
                csv_url = csv_url or u
            elif "zip" in name or name.endswith(".zip"):
                zip_url = zip_url or u
    return html_url, csv_url, zip_url


@registry.register(
    name="bepipred3_prediction",
    description=(
        "BepiPred-3.0：基于蛋白语言模型的 B 细胞线性/构象表位预测。"
        "当用户需要 B 细胞表位预测、BepiPred3、免疫表位分析时调用。"
        "输入可为 FASTA 文本或已上传的 .fasta/.fa/.faa/.ffn 文件路径。"
    ),
    category="Proteomics",
    output_type="json",
)
async def bepipred3_prediction(
    sequence_or_path: str,
    top_epitope_percentage_cutoff: TopEpitopePercentageCutoff = "top_20",
    use_sequential_smoothing: bool = False,
) -> Dict[str, Any]:
    """
    调用磐石 HTTP API 对蛋白质序列进行 BepiPred-3.0 B 细胞表位预测。

    **何时调用（给大模型）**：
    当用户要求进行 **B 细胞表位预测**、**BepiPred3 分析**、**免疫表位/抗原表位预测**（且语境为蛋白序列而非质谱定量矩阵）时，**应直接调用本工具** `bepipred3_prediction`。

    **参数**：
    - `sequence_or_path`：FASTA 文本，或服务器上存在的 FASTA 文件路径。
    - `top_epitope_percentage_cutoff`：`top_20` | `top_50` | `top_70` | `all`。
    - `use_sequential_smoothing`：是否顺序平滑。

    环境变量 `BEPIPRED3_PREDICT_URL` 可覆盖默认预测端点。
    """
    try:
        fasta_content = _resolve_fasta_content(sequence_or_path)
    except Exception as e:
        logger.exception("BepiPred3 输入解析失败: %s", e)
        return {"status": "error", "message": f"输入解析失败: {e}"}

    preview = fasta_content[:50] + "..." if len(fasta_content) > 50 else fasta_content
    logger.info(
        "正在执行 BepiPred3 预测（磐石 API），输入: %s, 阈值: %s, 平滑: %s",
        preview,
        top_epitope_percentage_cutoff,
        use_sequential_smoothing,
    )

    url = (os.getenv("BEPIPRED3_PREDICT_URL") or "").strip() or DEFAULT_BEPIPRED3_PREDICT_URL
    payload: Dict[str, Any] = {
        "url_or_content": fasta_content,
        "top_epitope_percentage_cutoff": top_epitope_percentage_cutoff,
        "use_sequential_smoothing": use_sequential_smoothing,
        "prediction_mode": "vt_pred",
    }

    async with httpx.AsyncClient(timeout=300.0) as client:
        for attempt in range(1, _MAX_RETRIES + 1):
            try:
                response = await client.post(url, json=payload)
            except httpx.TimeoutException as e:
                logger.error("BepiPred3 API 超时 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e)
                if attempt < _MAX_RETRIES:
                    logger.warning(
                        "⚠️ [BepiPred3] 磐石服务器繁忙，等待 5 秒后进行第 %s 次重试 (%s/%s)...",
                        attempt + 1,
                        attempt,
                        _MAX_RETRIES,
                    )
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {"status": "error", "message": f"网络超时（>300s，已重试 {_MAX_RETRIES} 次）: {e}"}
            except httpx.RequestError as e:
                logger.error(
                    "BepiPred3 API 网络异常 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e, exc_info=True
                )
                if attempt < _MAX_RETRIES:
                    logger.warning(
                        "⚠️ [BepiPred3] 磐石服务器繁忙，等待 5 秒后进行第 %s 次重试 (%s/%s)...",
                        attempt + 1,
                        attempt,
                        _MAX_RETRIES,
                    )
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": f"网络错误（连接失败、DNS、TLS 或中断等，已重试 {_MAX_RETRIES} 次）: {e}",
                }
            except Exception as e:
                logger.exception("BepiPred3 API 请求未预期异常 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e)
                if attempt < _MAX_RETRIES:
                    logger.warning(
                        "⚠️ [BepiPred3] 请求异常，等待 5 秒后进行第 %s 次重试 (%s/%s)...",
                        attempt + 1,
                        attempt,
                        _MAX_RETRIES,
                    )
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {"status": "error", "message": f"请求异常: {e}"}

            raw_text = response.text or ""
            logger.info(
                "BepiPred3 磐石 API 原始响应 | attempt=%s | status_code=%s | response.text=%s",
                attempt,
                response.status_code,
                raw_text,
            )

            if response.status_code != 200:
                logger.error(
                    "BepiPred3 磐石 API HTTP 错误 | status_code=%s | response.text=%s",
                    response.status_code,
                    raw_text,
                )
                if response.status_code in _RETRYABLE_HTTP and attempt < _MAX_RETRIES:
                    logger.warning(
                        "⚠️ [BepiPred3] 磐石服务器繁忙，等待 5 秒后进行第 %s 次重试 (%s/%s)...",
                        attempt + 1,
                        attempt,
                        _MAX_RETRIES,
                    )
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": f"磐石 API 报错: HTTP {response.status_code}, 详情: {raw_text}",
                }

            try:
                result: Dict[str, Any] = response.json()
            except Exception as e:
                logger.error(
                    "BepiPred3 响应非 JSON | 解析异常=%s | response.text=%s",
                    e,
                    raw_text,
                    exc_info=True,
                )
                return {
                    "status": "error",
                    "message": f"磐石 API 报错: 响应非合法 JSON（{e}）, 原文: {raw_text}",
                }

            if not result.get("success") is True:
                detail = result.get("message") or result.get("error") or str(result)
                logger.error(
                    "BepiPred3 磐石 API 业务失败 | success!=True | 详情=%s | 完整 JSON=%s",
                    detail,
                    result,
                )
                return {
                    "status": "error",
                    "message": f"磐石 API 报错: success=False, 详情: {detail}；完整响应: {result}",
                }

            data = result.get("data")
            if not isinstance(data, dict):
                return {"status": "error", "message": "响应缺少 data 对象"}

            inner_status = data.get("status")
            if inner_status == "failed":
                error_msg = (
                    data.get("error_message")
                    or "磐石服务器内部计算失败（可能 CUDA 资源不足或设备忙）"
                )
                if _panshi_error_suggests_gpu_transient(error_msg) and attempt < _MAX_RETRIES:
                    logger.warning(
                        "⚠️ [BepiPred3] 磐石服务器繁忙，等待 5 秒后进行第 %s 次重试 (%s/%s)...",
                        attempt + 1,
                        attempt,
                        _MAX_RETRIES,
                    )
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                if _panshi_error_suggests_gpu_transient(error_msg):
                    logger.error("❌ 磐石真实 API 节点 GPU 资源耗尽，重试失败。")
                    return {
                        "status": "error",
                        "message": _GPU_BUSY_USER_MESSAGE,
                    }
                logger.error("❌ [BepiPred3] 磐石计算失败 (data.status=failed): %s", error_msg)
                return {
                    "status": "error",
                    "message": f"磐石服务器计算失败: {error_msg}",
                }

            if inner_status == "completed":
                html_url, csv_url, zip_url = _extract_output_urls(data)
                if not (html_url or csv_url or zip_url):
                    logger.warning(
                        "BepiPred3 data.status=completed 但未解析到 output_files 链接: data=%s", data
                    )
                    return {
                        "status": "error",
                        "message": (
                            "磐石返回计算完成，但未提供可下载的 html/csv/zip 链接（output_files 为空或格式异常）。"
                        ),
                    }
                return {
                    "status": "success",
                    "html_url": html_url or "",
                    "csv_url": csv_url or "",
                    "zip_url": zip_url or "",
                    "message": "预测成功",
                }

            if inner_status == "running":
                return {
                    "status": "error",
                    "message": f"任务仍在运行中，状态: {inner_status}",
                }

            if inner_status in (None, ""):
                html_url, csv_url, zip_url = _extract_output_urls(data)
                if html_url or csv_url or zip_url:
                    logger.warning(
                        "BepiPred3 响应缺少 data.status，但存在 output_files，按成功降级处理"
                    )
                    return {
                        "status": "success",
                        "html_url": html_url or "",
                        "csv_url": csv_url or "",
                        "zip_url": zip_url or "",
                        "message": "预测成功",
                    }

            return {
                "status": "error",
                "message": f"磐石任务状态异常或未就绪: {inner_status!r}",
            }

