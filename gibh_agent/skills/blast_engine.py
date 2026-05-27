# -*- coding: utf-8 -*-
"""
BLAST 双引擎路由：本地库（极速） / NCBI 远程（排队）。

本地库检测顺序：
  1. 环境变量 LOCAL_BLASTDB_PATH
  2. /app/gibh_agent/blast_db（容器默认）
  3. 仓库内 gibh_agent/blast_db
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.skills._skill_common import resolve_remote_blastn_options

# 面向用户的产品级文案（禁止向前端暴露 TimeoutError / Stack Trace）
MSG_REMOTE_TIMEOUT = (
    "NCBI 官方服务器当前排队拥挤或网络存在波动，远程比对已超时。"
    "建议您：1. 稍后重试；2. 换用更短的序列；"
    "3. 若需要秒级响应，请联系管理员在服务器挂载本地 BLAST 数据库"
    "（设置 LOCAL_BLASTDB_PATH 或将库文件放入 gibh_agent/blast_db/）。"
)

MSG_REMOTE_NETWORK = (
    "无法稳定连接 NCBI 远程 BLAST 服务（网络波动或防火墙限制）。"
    "建议您稍后重试；若需稳定快速比对，请由管理员挂载本地 nt/core_nt 数据库。"
)

MSG_LOCAL_MISSING = (
    "当前服务器未检测到本地核酸 BLAST 数据库，且已关闭远程比对模式。"
    "请联系管理员配置 LOCAL_BLASTDB_PATH 或 gibh_agent/blast_db/ 下的 nt 库，"
    "或允许使用远程 NCBI（use_remote=true）。"
)

MSG_LOCAL_TIMEOUT = (
    "本地核酸比对耗时过长，可能数据库未完整挂载或磁盘负载过高。"
    "请联系管理员检查 gibh_agent/blast_db/ 是否包含完整的 nt.nhr/.nin/.nsq 文件。"
)


def blast_db_search_roots() -> List[Path]:
    roots: List[Path] = []
    env_root = (os.getenv("LOCAL_BLASTDB_PATH") or "").strip()
    if env_root:
        roots.append(Path(env_root))
    roots.extend(
        [
            Path("/app/gibh_agent/blast_db"),
            Path(__file__).resolve().parents[1] / "blast_db",
        ]
    )
    seen: set[str] = set()
    out: List[Path] = []
    for r in roots:
        key = str(r.resolve())
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def local_blast_database(db_name: str) -> Optional[Dict[str, str]]:
    """
    检测本地 BLAST 库是否可用（存在 .nhr/.nin/.nsq 或 .phr/.pin/.psq 之一组合）。

    Returns:
        db_arg: 传给 blastn/blastp -db 的前缀路径
        blastdb_dir: BLASTDB 环境变量目录
        db_name: 实际库名
    """
    name = (db_name or "").strip()
    if not name:
        return None
    for root in blast_db_search_roots():
        if not root.is_dir():
            continue
        prefix = root / name
        nuc = (root / f"{name}.nhr").is_file() or prefix.with_suffix(".nhr").is_file()
        prot = (root / f"{name}.phr").is_file() or prefix.with_suffix(".phr").is_file()
        if nuc or prot:
            return {
                "db_arg": str(prefix),
                "blastdb_dir": str(root.resolve()),
                "db_name": name,
            }
    return None


def _local_db_candidates(requested: str) -> List[str]:
    """按优先级尝试本地库名（请求 core_nt 时可回退 nt）。"""
    db = (requested or "nt").strip() or "nt"
    order = [db]
    if db == "core_nt":
        order.append("nt")
    elif db == "nt":
        order.append("core_nt")
    seen: set[str] = set()
    out: List[str] = []
    for d in order:
        if d not in seen:
            seen.add(d)
            out.append(d)
    return out


def plan_blastn_execution(
    seq_len: int,
    database: str,
    blast_task: str,
    *,
    use_remote: bool,
) -> Dict[str, Any]:
    """
    双引擎路由计划。

    Returns dict keys:
        ok: bool
        engine: local | remote | unavailable
        database, blast_task, timeout_s, db_arg, use_remote
        blastdb_dir (local only)
        message (when ok=False)
    """
    requested = (database or "nt").strip() or "nt"
    task_in = (blast_task or "").strip()

    for cand in _local_db_candidates(requested):
        loc = local_blast_database(cand)
        if loc:
            task = task_in
            if not task and seq_len < 80:
                task = "blastn-short"
            timeout_s = 180 if seq_len < 500 else 600
            return {
                "ok": True,
                "engine": "local",
                "database": loc["db_name"],
                "blast_task": task,
                "timeout_s": timeout_s,
                "db_arg": loc["db_arg"],
                "blastdb_dir": loc["blastdb_dir"],
                "use_remote": False,
            }

    if not use_remote:
        return {"ok": False, "engine": "unavailable", "message": MSG_LOCAL_MISSING}

    db, task, _ = resolve_remote_blastn_options(
        seq_len, requested, task_in, use_remote=True
    )
    return {
        "ok": True,
        "engine": "remote",
        "database": db,
        "blast_task": task,
        "timeout_s": 600,
        "db_arg": db,
        "blastdb_dir": "",
        "use_remote": True,
    }


def user_facing_blastn_error(
    *,
    code: int,
    stderr: str,
    stdout: str,
    engine: str,
    timeout_s: int,
) -> str:
    """将底层 blast 退出信息转为产品级 message（不含堆栈）。"""
    blob = f"{stderr or ''}\n{stdout or ''}"
    lower = blob.lower()

    if code == -1 and "[timeout]" in blob:
        return MSG_REMOTE_TIMEOUT if engine == "remote" else MSG_LOCAL_TIMEOUT

    network_markers = (
        "network",
        "connection",
        "could not connect",
        "failed to connect",
        "name or service not known",
        "temporary failure",
        "ssl",
        "http",
    )
    if engine == "remote" and any(m in lower for m in network_markers):
        return MSG_REMOTE_NETWORK

    if "traceback" in lower or "exception" in lower:
        return (
            "核酸比对引擎内部异常，任务已安全终止。"
            "请稍后重试；若反复失败，请联系管理员检查 BLAST 安装或本地数据库挂载。"
        )

    tail = blob.strip()
    if len(tail) > 240:
        tail = tail[-240:]
    if not tail:
        tail = "blastn 进程异常退出"
    return (
        f"核酸比对未能完成（{tail}）。"
        "建议您稍后重试；远程模式受 NCBI 排队影响，"
        "稳定场景请由管理员配置本地 BLAST 数据库。"
    )
