"""
蛋白质序列层工具：B 细胞表位预测（BepiPred3）等。

与 `proteomics_tools.py`（质谱矩阵统计）互补；本模块聚焦 FASTA/序列级任务。

默认在独立子进程中调用已部署的 BepiPred-3.0 CLI，将产物写入 RESULTS_DIR 并返回与历史远程 API 一致的
html_url / csv_url / zip_url 契约。若设置 BEPIPRED3_USE_REMOTE_API=1，则回退为磐石 HTTP 预测服务。
"""
from __future__ import annotations

import asyncio
import logging
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple

import httpx

from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution

logger = logging.getLogger(__name__)

TopEpitopePercentageCutoff = Literal["top_20", "top_50", "top_70", "all"]

DEFAULT_BEPIPRED3_PREDICT_URL = "http://120.220.102.26:38089/predict"

_RETRYABLE_HTTP = frozenset({500, 502, 503, 504})
_MAX_RETRIES = 3
_RETRY_DELAY_SEC = 5.0

_GPU_BUSY_USER_MESSAGE = (
    "远程计算服务当前资源繁忙（可能由于排队人数过多或 GPU 满载）。请稍后重试。"
)

# 解释器探测较耗时，成功后按 bp_root 缓存（Worker 重启后失效）
_bepipred_python_cache: Dict[str, str] = {}
_bootstrap_lock = asyncio.Lock()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_bepipred_root() -> Path:
    override = (os.getenv("BEPIPRED3_ROOT") or "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return (_repo_root() / "third_party" / "BepiPred-3.0").resolve()


def _which_python_executable(cmd: str) -> Optional[str]:
    p = Path(cmd)
    if p.is_file():
        return str(p.resolve())
    found = shutil.which(cmd)
    return found


def _python_can_run_bepipred(exe: str, bp_root: Path) -> bool:
    """子进程需能 import torch、esm（fair-esm），且 cwd 下可加载本地 bp3。"""
    try:
        r = subprocess.run(
            [exe, "-c", "import torch; import esm"],
            cwd=str(bp_root),
            capture_output=True,
            text=True,
            timeout=int(os.getenv("BEPIPRED3_PYTHON_PROBE_TIMEOUT", "120")),
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        if r.returncode != 0:
            logger.debug("BepiPred python 探测失败 %s: %s", exe, (r.stderr or r.stdout or "")[:500])
        return r.returncode == 0
    except Exception as exc:
        logger.debug("BepiPred python 探测异常 %s: %s", exe, exc)
        return False


def _pick_bepipred_python(bp_root: Path) -> Tuple[Optional[str], str]:
    """
    按顺序探测可运行 BepiPred 的解释器，避免容器内误用无 torch 的 python3。
    候选：BEPIPRED3_PYTHON → <ROOT>/.venv/bin/python → venv/bin/python → 当前进程解释器 → PATH 中 python3/python。
    """
    cache_off = (os.getenv("BEPIPRED3_PYTHON_PROBE_CACHE") or "").strip().lower() in ("0", "false", "no")
    key = str(bp_root.resolve())
    if not cache_off and key in _bepipred_python_cache:
        return _bepipred_python_cache[key], ""

    candidates: List[str] = []

    env_py = (os.getenv("BEPIPRED3_PYTHON") or "").strip()
    if env_py:
        w = _which_python_executable(env_py)
        if w:
            candidates.append(w)

    _opt_py = Path("/opt/bepipred3-venv/bin/python")
    if _opt_py.is_file():
        s = str(_opt_py.resolve())
        if s not in candidates:
            candidates.append(s)

    for rel in (bp_root / ".venv" / "bin" / "python", bp_root / "venv" / "bin" / "python"):
        if rel.is_file():
            s = str(rel.resolve())
            if s not in candidates:
                candidates.append(s)

    main_py = (sys.executable or "").strip()
    if main_py and main_py not in candidates:
        candidates.append(main_py)

    fb = (os.getenv("BEPIPRED3_FALLBACK_PYTHON") or "").strip()
    if fb:
        w = _which_python_executable(fb)
        if w and w not in candidates:
            candidates.append(w)

    for name in ("python3", "python"):
        w = shutil.which(name)
        if w and w not in candidates:
            candidates.append(w)

    tried: List[str] = []
    for exe in candidates:
        if not exe:
            continue
        tried.append(exe)
        if _python_can_run_bepipred(exe, bp_root):
            logger.info("BepiPred-3.0 选用解释器: %s", exe)
            if not cache_off:
                _bepipred_python_cache[key] = exe
            return exe, ""

    detail = "; ".join(tried[:10]) if tried else "(无候选)"
    hints: List[str] = []
    if env_py and not _which_python_executable(env_py):
        hints.append(f"环境变量 BEPIPRED3_PYTHON={env_py} 指向的文件不存在或不可执行（该路径未参与探测）")
    if not _opt_py.is_file():
        hints.append("镜像内 /opt/bepipred3-venv 不存在（需重新 docker compose build api-server，或等待首次请求惰性安装）")
    hint_txt = (" " + "；".join(hints)) if hints else ""
    msg = (
        "未找到已安装 PyTorch 与 fair-esm（import torch; import esm）的 Python。"
        f"已尝试: {detail}。{hint_txt}"
        f"可在 {bp_root} 下执行: python3 -m venv .venv && .venv/bin/pip install torch torchvision -r requirements.txt "
        "（torch 需按官方 index），或设置 BEPIPRED3_PYTHON。默认会在目录可写时自动创建 .venv（BEPIPRED3_LAZY_BOOTSTRAP=0 可关闭）。"
    )
    return None, msg


def _lazy_bootstrap_enabled() -> bool:
    return (os.getenv("BEPIPRED3_LAZY_BOOTSTRAP") or "1").strip().lower() not in ("0", "false", "no", "off")


def _bootstrap_bepipred_under_bp_root(bp_root: Path) -> Tuple[bool, str]:
    """
    在 bp_root 下创建/重建 .venv 并 pip 安装 torch + requirements.txt。
    要求 bp_root 对当前进程可写；在 Docker 挂载 .:/app 且目录属主为 root 时可能失败。
    """
    req = bp_root / "requirements.txt"
    if not req.is_file():
        return False, f"缺少 {req}"
    if not os.access(bp_root, os.W_OK):
        return False, f"目录不可写，无法自动创建 .venv: {bp_root}"

    vdir = bp_root / ".venv"
    vpy = vdir / "bin" / "python"
    pip = vdir / "bin" / "pip"
    if vpy.is_file() and _python_can_run_bepipred(str(vpy), bp_root):
        return True, ""

    if vdir.exists():
        logger.warning("BepiPred 现有 .venv 无法 import torch/esm，将删除后重建: %s", vdir)
        try:
            shutil.rmtree(vdir)
        except OSError as exc:
            return False, f"无法删除旧 .venv: {exc}"

    base_py = sys.executable or shutil.which("python3") or "python3"
    try:
        subprocess.run(
            [base_py, "-m", "venv", str(vdir)],
            cwd=str(bp_root),
            check=True,
            timeout=180,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        err = getattr(exc, "stderr", None) or getattr(exc, "stdout", None) or str(exc)
        return False, f"python -m venv 失败: {err[:2000]}"

    if not pip.is_file():
        return False, "venv 创建后未找到 bin/pip"

    torch_idx = (os.getenv("BEPIPRED_TORCH_INDEX") or "https://download.pytorch.org/whl/cpu").strip()
    steps = [
        ([str(pip), "install", "-U", "pip"], 600),
        ([str(pip), "install", "torch", "torchvision", "--index-url", torch_idx], 3600),
        ([str(pip), "install", "-r", str(req)], 3600),
    ]
    for cmd, tmo in steps:
        try:
            r = subprocess.run(
                cmd,
                cwd=str(bp_root),
                capture_output=True,
                text=True,
                timeout=tmo,
                env={**os.environ, "PYTHONUNBUFFERED": "1"},
            )
            if r.returncode != 0:
                tail = (r.stderr or r.stdout or "")[-4000:]
                return False, f"命令失败 {' '.join(cmd[:6])}… : {tail}"
        except subprocess.TimeoutExpired:
            return False, f"超时（>{tmo}s）: {' '.join(cmd[:6])}"

    if not vpy.is_file() or not _python_can_run_bepipred(str(vpy), bp_root):
        return False, "安装完成但仍无法 import torch; import esm"
    logger.info("BepiPred 惰性安装完成: %s", vpy)
    return True, ""


def _results_dir() -> Path:
    raw = (os.getenv("RESULTS_DIR") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    # 容器内与 server.py 一致；裸机未配置时写入仓库 results/，避免默认 /app/results 无权限
    docker_results = Path("/app/results")
    if docker_results.parent.is_dir():
        return docker_results.resolve()
    return (_repo_root() / "results").resolve()


def _public_path_url(rel_path: str) -> str:
    """返回前端可访问路径：默认 /results/...，可选绝对前缀。"""
    rel = rel_path if rel_path.startswith("/") else "/" + rel_path.lstrip("/")
    base = (os.getenv("PUBLIC_RESULTS_BASE_URL") or "").rstrip("/")
    return f"{base}{rel}" if base else rel


def _cutoff_to_cli_top(cutoff: TopEpitopePercentageCutoff) -> float:
    return {"top_20": 0.2, "top_50": 0.5, "top_70": 0.7, "all": 1.0}[cutoff]


def _resolve_fasta_content(sequence_or_path: str) -> str:
    raw = (sequence_or_path or "").strip()
    if not raw:
        raise ValueError("sequence_or_path 为空")
    if os.path.isfile(raw):
        with open(raw, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    return raw


def _extract_output_urls(data: Dict[str, Any]) -> Tuple[Optional[str], Optional[str], Optional[str]]:
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


def _remote_gpu_error_suggests_transient(error_message: str) -> bool:
    if not error_message:
        return False
    low = error_message.lower()
    return "cuda" in low or "busy" in low or "unavailable" in low


async def _bepipred_via_local_cli(
    seq_input: str,
    top_epitope_percentage_cutoff: TopEpitopePercentageCutoff,
    use_sequential_smoothing: bool,
) -> Dict[str, Any]:
    bp_root = _default_bepipred_root()
    cli = bp_root / "bepipred3_CLI.py"
    if not cli.is_file():
        return {
            "status": "error",
            "message": (
                f"未找到本地 BepiPred-3.0 入口脚本: {cli}。"
                "请克隆官方仓库至 third_party/BepiPred-3.0 或设置环境变量 BEPIPRED3_ROOT。"
            ),
        }

    py, py_err = await asyncio.to_thread(_pick_bepipred_python, bp_root)
    if not py and _lazy_bootstrap_enabled():
        async with _bootstrap_lock:
            py, py_err = await asyncio.to_thread(_pick_bepipred_python, bp_root)
            if not py:
                logger.warning("BepiPred 无可用解释器，惰性安装到 %s/.venv …", bp_root)
                ok, berr = await asyncio.to_thread(_bootstrap_bepipred_under_bp_root, bp_root)
                if ok:
                    _bepipred_python_cache.pop(str(bp_root.resolve()), None)
                    py, py_err = await asyncio.to_thread(_pick_bepipred_python, bp_root)
                else:
                    py_err = f"{py_err}\n惰性安装失败: {berr}"

    if not py:
        return {"status": "error", "message": py_err}
    run_id = uuid.uuid4().hex[:16]
    base = _results_dir() / "bepipred3" / run_id
    out_dir = base / "output"
    esm_dir = base / "esm_encodings"
    input_fasta = base / "input.fasta"
    try:
        base.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建结果目录: {exc}"}

    input_fasta.write_text(seq_input.strip() + "\n", encoding="utf-8")

    top = _cutoff_to_cli_top(top_epitope_percentage_cutoff)
    cmd: list[str] = [
        py,
        str(cli),
        "-i",
        str(input_fasta),
        "-o",
        str(out_dir),
        "-pred",
        "vt_pred",
        "-esm_dir",
        str(esm_dir),
        "-top",
        str(top),
    ]
    if use_sequential_smoothing:
        cmd.append("-plot_linear_epitope_scores")

    timeout = float(os.getenv("BEPIPRED3_SUBPROCESS_TIMEOUT", "3600"))
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    cwd = str(bp_root)

    logger.info("BepiPred3 本地子进程: %s", " ".join(cmd[:8]) + " ...")

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=cwd,
        env=env,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    try:
        stdout_b, stderr_b = await asyncio.wait_for(proc.communicate(), timeout=timeout)
    except asyncio.TimeoutError:
        proc.kill()
        await proc.communicate()
        return {"status": "error", "message": f"BepiPred-3.0 子进程超时（>{timeout}s）"}

    stderr_text = (stderr_b or b"").decode("utf-8", errors="replace").strip()
    if stderr_text:
        logger.warning("BepiPred3 stderr: %s", stderr_text[:8000])

    if proc.returncode != 0:
        return {
            "status": "error",
            "message": f"BepiPred-3.0 退出码 {proc.returncode}；stderr: {stderr_text[:4000]}",
        }

    html_p = out_dir / "output_interactive_figures.html"
    csv_p = out_dir / "raw_output.csv"
    if not html_p.is_file():
        return {
            "status": "error",
            "message": "本地预测已完成但未生成 output_interactive_figures.html，请检查输入与依赖环境。",
        }

    zip_base = base / "bepipred3_results"
    try:
        shutil.make_archive(str(zip_base), "zip", root_dir=str(out_dir))
    except OSError as exc:
        logger.warning("打包 zip 失败（仍返回 html/csv）: %s", exc)

    zip_p = Path(str(zip_base) + ".zip")
    rel = f"/results/bepipred3/{run_id}"
    html_url = _public_path_url(f"{rel}/output/{html_p.name}")
    csv_url = _public_path_url(f"{rel}/output/{csv_p.name}") if csv_p.is_file() else ""
    zip_url = _public_path_url(f"{rel}/{zip_p.name}") if zip_p.is_file() else ""

    if not (html_url or csv_url or zip_url):
        return {"status": "error", "message": "未能构造结果下载路径"}

    return {
        "status": "success",
        "html_url": html_url,
        "csv_url": csv_url or "",
        "zip_url": zip_url or "",
        "message": "预测成功（本地 BepiPred-3.0）",
    }


async def _bepipred_via_remote_api(
    seq_input: str,
    top_epitope_percentage_cutoff: TopEpitopePercentageCutoff,
    use_sequential_smoothing: bool,
) -> Dict[str, Any]:
    url = (os.getenv("BEPIPRED3_PREDICT_URL") or "").strip() or DEFAULT_BEPIPRED3_PREDICT_URL
    payload: Dict[str, Any] = {
        "url_or_content": seq_input,
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
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {"status": "error", "message": f"网络超时（>300s，已重试 {_MAX_RETRIES} 次）: {e}"}
            except httpx.RequestError as e:
                logger.error(
                    "BepiPred3 API 网络异常 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e, exc_info=True
                )
                if attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": f"网络错误（已重试 {_MAX_RETRIES} 次）: {e}",
                }
            except Exception as e:
                logger.exception("BepiPred3 API 请求未预期异常 (尝试 %s/%s): %s", attempt, _MAX_RETRIES, e)
                if attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {"status": "error", "message": f"请求异常: {e}"}

            raw_text = response.text or ""

            if response.status_code != 200:
                if response.status_code in _RETRYABLE_HTTP and attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                return {
                    "status": "error",
                    "message": f"预测服务报错: HTTP {response.status_code}, 详情: {raw_text}",
                }

            try:
                result: Dict[str, Any] = response.json()
            except Exception as e:
                return {
                    "status": "error",
                    "message": f"预测服务报错: 响应非合法 JSON（{e}）, 原文: {raw_text}",
                }

            if not result.get("success") is True:
                detail = result.get("message") or result.get("error") or str(result)
                return {
                    "status": "error",
                    "message": f"预测服务报错: success=False, 详情: {detail}；完整响应: {result}",
                }

            data = result.get("data")
            if not isinstance(data, dict):
                return {"status": "error", "message": "响应缺少 data 对象"}

            inner_status = data.get("status")
            if inner_status == "failed":
                error_msg = (
                    data.get("error_message")
                    or "远程服务内部计算失败（可能 CUDA 资源不足或设备忙）"
                )
                if _remote_gpu_error_suggests_transient(error_msg) and attempt < _MAX_RETRIES:
                    await asyncio.sleep(_RETRY_DELAY_SEC)
                    continue
                if _remote_gpu_error_suggests_transient(error_msg):
                    return {"status": "error", "message": _GPU_BUSY_USER_MESSAGE}
                return {"status": "error", "message": f"远程计算服务失败: {error_msg}"}

            if inner_status == "completed":
                html_url, csv_url, zip_url = _extract_output_urls(data)
                if not (html_url or csv_url or zip_url):
                    return {
                        "status": "error",
                        "message": (
                            "服务返回计算完成，但未提供可下载的 html/csv/zip 链接（output_files 为空或格式异常）。"
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
                return {"status": "error", "message": f"任务仍在运行中，状态: {inner_status}"}

            if inner_status in (None, ""):
                html_url, csv_url, zip_url = _extract_output_urls(data)
                if html_url or csv_url or zip_url:
                    return {
                        "status": "success",
                        "html_url": html_url or "",
                        "csv_url": csv_url or "",
                        "zip_url": zip_url or "",
                        "message": "预测成功",
                    }

            return {
                "status": "error",
                "message": f"预测任务状态异常或未就绪: {inner_status!r}",
            }

    return {"status": "error", "message": "预测失败（未知原因）"}


def _use_remote_api() -> bool:
    return (os.getenv("BEPIPRED3_USE_REMOTE_API") or "").strip().lower() in ("1", "true", "yes", "on")


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
@safe_tool_execution
async def bepipred3_prediction(
    sequence_or_path: str,
    top_epitope_percentage_cutoff: TopEpitopePercentageCutoff = "top_20",
    use_sequential_smoothing: bool = False,
) -> Dict[str, Any]:
    """
    默认：在独立子进程中执行本地已部署的 BepiPred-3.0 CLI，结果写入 RESULTS_DIR 并通过静态 /results 暴露。

    设置环境变量 BEPIPRED3_USE_REMOTE_API=1 时，改为调用远程 HTTP 预测服务（BEPIPRED3_PREDICT_URL）。

    环境变量（本地模式）：
    - BEPIPRED3_ROOT：BepiPred-3.0 根目录（默认仓库内 third_party/BepiPred-3.0）
    - BEPIPRED3_PYTHON：解释器路径（默认 <ROOT>/.venv/bin/python 或 python3）
    - BEPIPRED3_SUBPROCESS_TIMEOUT：秒，默认 3600
    - RESULTS_DIR：与 FastAPI 静态挂载一致，默认 /app/results
    - PUBLIC_RESULTS_BASE_URL：可选，为链接加绝对前缀
    """
    try:
        fasta_content = _resolve_fasta_content(sequence_or_path)
    except Exception as e:
        logger.exception("BepiPred3 输入解析失败: %s", e)
        return {"status": "error", "message": f"输入解析失败: {e}"}

    seq_input = fasta_content.strip()
    if seq_input and not seq_input.startswith("http://") and not seq_input.startswith("https://"):
        if not seq_input.startswith(">"):
            logger.info("检测到序列缺失 FASTA 头部 '>'，自动补全。")
            seq_input = f">Sequence\n{seq_input}"

    preview = seq_input[:50] + "..." if len(seq_input) > 50 else seq_input
    mode = "远程 API" if _use_remote_api() else "本地 CLI"
    logger.info(
        "BepiPred3 预测（%s），输入: %s, 阈值: %s, 平滑: %s",
        mode,
        preview,
        top_epitope_percentage_cutoff,
        use_sequential_smoothing,
    )

    try:
        _results_dir().mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return {"status": "error", "message": f"结果目录不可用: {exc}"}

    if _use_remote_api():
        return await _bepipred_via_remote_api(
            seq_input, top_epitope_percentage_cutoff, use_sequential_smoothing
        )
    return await _bepipred_via_local_cli(
        seq_input, top_epitope_percentage_cutoff, use_sequential_smoothing
    )
