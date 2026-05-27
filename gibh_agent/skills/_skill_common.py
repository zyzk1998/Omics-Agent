# -*- coding: utf-8 -*-
"""技能脚本共用辅助（非 BaseSkill，不会被 skill_registry 注册）。"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote

import httpx

DEFAULT_HTTP_TIMEOUT = float(os.environ.get("GIBH_SKILL_HTTP_TIMEOUT", "60"))


def ok(message: str, **data: Any) -> Dict[str, Any]:
    payload = dict(data)
    md = _format_data_markdown(message, payload)
    out: Dict[str, Any] = {"status": "success", "message": message, "data": payload}
    if md:
        out["markdown"] = md
    return out


def _format_data_markdown(message: str, data: Dict[str, Any]) -> str:
    """右栏工作台简要 Markdown（与 DRT/chem 范式对齐的轻量版）。"""
    import json

    lines = [f"### {message}\n"]
    _table_row_keys = ("hits", "results", "rows", "records", "matches", "experiments", "alignments")
    for key in ("hits", "cids", "experiment", "result", "mol_block", "sdf_text"):
        if key in data and data[key] is not None:
            val = data[key]
            if (
                key in _table_row_keys
                and isinstance(val, list)
                and val
                and isinstance(val[0], dict)
            ):
                continue
            if isinstance(val, (list, dict)) and len(str(val)) > 1200:
                lines.append(f"#### {key}\n\n```json\n{json.dumps(val, ensure_ascii=False, indent=2)[:4000]}\n```\n")
            else:
                lines.append(f"**{key}**: `{val}`\n")
    if len(lines) <= 1:
        lines.append(f"```json\n{json.dumps(data, ensure_ascii=False, indent=2)[:3000]}\n```\n")
    return "\n".join(lines)


def err(message: str, **data: Any) -> Dict[str, Any]:
    out: Dict[str, Any] = {"status": "error", "message": message}
    if data:
        out["data"] = dict(data)
    return out


def http_get_json(url: str, *, params: Optional[Dict[str, Any]] = None, timeout: Optional[float] = None) -> Any:
    with httpx.Client(timeout=timeout or DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        return resp.json()


def resolve_sequence_or_fasta(
    sequence_or_path: str = "",
    *,
    sequence_text: str = "",
) -> Tuple[Optional[str], Optional[str]]:
    """
    从纯文本序列或 FASTA 文件路径解析单条查询序列。

    Returns:
        (sequence, error_message)
    """
    raw = (sequence_text or sequence_or_path or "").strip()
    if not raw:
        return None, "请提供 sequence_text 或 sequence_or_path（FASTA 文件路径）"

    path = Path(raw)
    if path.is_file():
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            return None, f"无法读取 FASTA 文件: {exc}"
        lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith(";")]
        seq_lines: List[str] = []
        for ln in lines:
            if ln.startswith(">"):
                if seq_lines:
                    break
                continue
            seq_lines.append(re.sub(r"\s+", "", ln))
        if not seq_lines:
            return None, "FASTA 文件中未找到有效序列"
        return "".join(seq_lines).upper(), None

    if path.suffix.lower() in {".fa", ".fasta", ".fna", ".faa"} and not path.exists():
        return None, f"FASTA 文件不存在: {raw}"

    cleaned = re.sub(r"\s+", "", raw)
    if cleaned.startswith(">"):
        parts = cleaned.split(">", 1)
        if len(parts) > 1:
            body = parts[1]
            if "\n" in body or " " in body:
                body = body.split(None, 1)[-1] if " " in body else body.split("\n", 1)[-1]
            cleaned = re.sub(r"[^A-Za-z*]", "", body)
    else:
        cleaned = re.sub(r"[^A-Za-z*]", "", cleaned)
    if len(cleaned) < 5:
        return None, "序列过短或格式无效（至少 5 个残基/碱基）"
    return cleaned.upper(), None


def resolve_remote_blastn_options(
    seq_len: int,
    database: str,
    blast_task: str = "",
    *,
    use_remote: bool,
) -> Tuple[str, str, int]:
    """
    远程 blastn 参数：短序列用 core_nt + blastn-short，避免默认 nt 全库排队数分钟。
    Returns: (database, task, timeout_seconds)
    """
    db = (database or "nt").strip() or "nt"
    task = (blast_task or "").strip()
    if not use_remote:
        return db, task, 600
    if db == "nt" and seq_len < 200:
        db = "core_nt"
    if not task:
        if seq_len < 80:
            task = "blastn-short"
        elif seq_len < 500:
            task = "megablast"
    # 跨境访问 NCBI 远程队列：短序列 core_nt+blastn-short 通常 1–3 分钟，高峰更长
    if seq_len < 80:
        timeout_s = 300
    elif seq_len < 500:
        timeout_s = 300
    else:
        timeout_s = 600
    return db, task, timeout_s


def resolve_remote_blastp_options(
    seq_len: int,
    database: str,
    blast_task: str = "",
    *,
    use_remote: bool,
) -> Tuple[str, str, int]:
    """远程 blastp：短肽用 swissprot + blastp-short（比默认 nr 全库快）。"""
    db = (database or "nr").strip() or "nr"
    task = (blast_task or "").strip()
    if not use_remote:
        return db, task, 600
    if db == "nr" and seq_len < 50:
        db = "swissprot"
    if not task and seq_len < 30:
        task = "blastp-short"
    if seq_len < 30:
        timeout_s = 300
    elif seq_len < 200:
        timeout_s = 300
    else:
        timeout_s = 600
    return db, task, timeout_s


def write_temp_fasta(sequence: str, *, prefix: str = "query") -> str:
    fd, path = tempfile.mkstemp(suffix=".fasta", prefix=prefix)
    os.close(fd)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{prefix}\n{sequence}\n")
    return path


def find_cli(name: str) -> Optional[str]:
    return shutil.which(name)


def run_cli(
    cmd: List[str],
    *,
    timeout_s: int = 300,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[int, str, str]:
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
            env=env,
        )
        return proc.returncode, proc.stdout or "", proc.stderr or ""
    except subprocess.TimeoutExpired as exc:
        partial_out = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
        partial_err = (exc.stderr or "") if isinstance(exc.stderr, str) else ""
        return (
            -1,
            partial_out,
            (partial_err or "")
            + f"\n[timeout] 命令超过 {timeout_s}s 未结束: {' '.join(cmd[:6])}...",
        )


def run_cli_with_progress(
    cmd: List[str],
    *,
    timeout_s: int = 300,
    env: Optional[Dict[str, str]] = None,
    progress_interval_s: float = 30.0,
    progress_message: str = "⏳ 仍在执行中，请稍候...",
) -> Tuple[int, str, str]:
    """
    阻塞 CLI 执行，并按间隔触发 ``emit_tool_log``（需调用方已绑定 tool_log_sink）。
    用于 NCBI ``-remote`` 等长耗时、无 stdout 心跳的外部进程。
    """
    from gibh_agent.core.tool_stream_log import emit_tool_log

    deadline = time.monotonic() + max(1, int(timeout_s))
    last_progress = time.monotonic()
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except OSError as exc:
        return -1, "", str(exc)

    stdout_chunks: List[str] = []
    stderr_chunks: List[str] = []
    try:
        while proc.poll() is None:
            now = time.monotonic()
            if now >= deadline:
                proc.kill()
                proc.wait(timeout=5)
                return (
                    -1,
                    "".join(stdout_chunks),
                    "".join(stderr_chunks)
                    + f"\n[timeout] 命令超过 {timeout_s}s 未结束: {' '.join(cmd[:6])}...",
                )
            if progress_interval_s > 0 and (now - last_progress) >= progress_interval_s:
                emit_tool_log(progress_message, state="running")
                last_progress = now
            time.sleep(0.5)
        out, err = proc.communicate(timeout=10)
        if out:
            stdout_chunks.append(out)
        if err:
            stderr_chunks.append(err)
        return proc.returncode or 0, "".join(stdout_chunks), "".join(stderr_chunks)
    except Exception as exc:
        try:
            proc.kill()
        except OSError:
            pass
        return -1, "".join(stdout_chunks), "".join(stderr_chunks) + f"\n{exc}"


def pubchem_smiles_url(smiles: str, suffix: str) -> str:
    encoded = quote(smiles, safe="")
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/{suffix}"
