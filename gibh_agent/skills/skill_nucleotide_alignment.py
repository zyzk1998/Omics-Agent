# -*- coding: utf-8 -*-
"""核酸序列 BLAST 比对 — 首发核心技能（本地/远程双引擎）。"""
from __future__ import annotations

import os
from typing import Any, Dict, List

from gibh_agent.core.tool_stream_log import emit_tool_log
from gibh_agent.skills._skill_common import (
    err,
    find_cli,
    ok,
    resolve_sequence_or_fasta,
    run_cli,
    run_cli_with_progress,
    write_temp_fasta,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.blast_engine import plan_blastn_execution, user_facing_blastn_error
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class NucleotideSequenceAlignmentSkill(BaseSkill):
    __abstractskill__ = False

    """
    核酸序列比对（NCBI BLAST blastn）— 本地/远程双引擎自适应。

    **执行策略（助手须在调用前向用户说明）**
    - 若服务器已挂载本地 nt/core_nt 库（``LOCAL_BLASTDB_PATH`` 或 ``gibh_agent/blast_db/``），
      系统自动**本地极速模式**，通常秒级～数分钟完成，耗时可预期。
    - 若无本地库，则走 **NCBI 远程排队**，耗时 1–10 分钟不等，受跨境网络与 NCBI 负载影响；
      请提前安抚用户「正在连接 NCBI 官方库排队，请耐心等待」，勿承诺秒回。

  参数:
        sequence_or_path: 核酸序列或 FASTA 文件路径。
        sequence_text: 序列文本（与 sequence_or_path 二选一）。
        database: 库名，默认 nt；远程短序列会自动优化为 core_nt。
        max_target_seqs: 最大命中条数（1–50）。
        use_remote: 无本地库时是否允许远程（默认 true）；有本地库时忽略此项并优先本地。
        blast_task: 可选，如 blastn-short；留空时按序列长度自动选择。
    """

    skill_id = "nucleotide_sequence_blast"
    display_name = "核酸序列比对"
    description = (
        "核酸 BLAST（blastn）：优先本地 nt 库（秒级稳定）；无本地库时走 NCBI 远程（可能排队数分钟）。"
        "调用前请向用户说明远程模式需耐心等待。"
    )
    category = "生物医药"
    sub_category = "数据分析"
    aliases = ["blastn", "核酸BLAST", "核酸序列比对", "BLAST核酸"]
    required_parameters = []
    tool_chain_key = ""
    __dependencies__ = ["apt:ncbi-blast+"]

    def execute(
        self,
        sequence_or_path: str = "",
        sequence_text: str = "",
        database: str = "nt",
        max_target_seqs: int = 10,
        use_remote: bool = True,
        blast_task: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "sequence_or_path": (sequence_or_path or "").strip(),
                "sequence_text": (sequence_text or "").strip(),
                "database": (database or "nt").strip(),
                "use_remote": use_remote,
                "max_target_seqs": max_target_seqs,
                "blast_task": (blast_task or kwargs.get("task") or "").strip(),
            },
        )
        seq, seq_err = resolve_sequence_or_fasta(
            str(filled.get("sequence_or_path") or ""),
            sequence_text=str(filled.get("sequence_text") or ""),
        )
        if seq_err:
            return err(seq_err)

        plan = plan_blastn_execution(
            len(seq),
            str(filled.get("database") or "nt"),
            str(filled.get("blast_task") or ""),
            use_remote=bool(filled.get("use_remote", True)),
        )
        if not plan.get("ok"):
            return err(str(plan.get("message") or "无法启动核酸比对"))

        engine = str(plan["engine"])
        database = str(plan["database"])
        blast_task = str(plan.get("blast_task") or "")
        timeout_s = int(plan["timeout_s"])
        db_arg = str(plan["db_arg"])
        use_remote_flag = bool(plan["use_remote"])
        blastdb_dir = str(plan.get("blastdb_dir") or "")

        blastn = find_cli("blastn")
        if not blastn:
            return err(
                "核酸比对服务暂不可用（未安装 blastn）。请联系管理员在计算环境中安装 ncbi-blast+。"
            )

        query_path = write_temp_fasta(seq, prefix="nucl_query")
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
        cmd: List[str] = [
            blastn,
            "-query",
            query_path,
            "-db",
            db_arg,
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            str(max(1, min(int(filled.get("max_target_seqs") or 10), 50))),
        ]
        if blast_task:
            cmd.extend(["-task", blast_task])
        if use_remote_flag:
            cmd.append("-remote")

        env = {**os.environ, "PYTHONUNBUFFERED": "1"}
        if blastdb_dir and not use_remote_flag:
            env["BLASTDB"] = blastdb_dir

        if use_remote_flag:
            emit_tool_log(
                "📡 正在向 NCBI 全球公共集群提交比对任务，因官方队列波动，通常需要 1~5 分钟，请耐心等待...",
                state="running",
            )

        try:
            if use_remote_flag:
                code, stdout, stderr = run_cli_with_progress(
                    cmd,
                    timeout_s=timeout_s,
                    env=env,
                    progress_interval_s=30.0,
                    progress_message="⏳ NCBI 远程排队中...（官方集群负载波动属正常现象，请勿关闭页面）",
                )
            else:
                code, stdout, stderr = run_cli(cmd, timeout_s=timeout_s, env=env)
        except Exception:
            return err(
                user_facing_blastn_error(
                    code=-1,
                    stderr="",
                    stdout="",
                    engine=engine,
                    timeout_s=timeout_s,
                ),
                engine=engine,
            )
        finally:
            try:
                os.unlink(query_path)
            except OSError:
                pass

        if code != 0:
            return err(
                user_facing_blastn_error(
                    code=code,
                    stderr=stderr,
                    stdout=stdout,
                    engine=engine,
                    timeout_s=timeout_s,
                ),
                engine=engine,
                database=database,
            )

        hits: List[Dict[str, str]] = []
        for line in stdout.splitlines():
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                continue
            hits.append(
                {
                    "qseqid": cols[0],
                    "sseqid": cols[1],
                    "pident": cols[2],
                    "length": cols[3],
                    "evalue": cols[11],
                    "bitscore": cols[12],
                    "title": cols[13] if len(cols) > 13 else "",
                }
            )

        if engine == "local":
            note = f"（本地库 {database}，稳定快速模式）"
        else:
            note = (
                f"（远程 NCBI：db={database}"
                + (f", task={blast_task}" if blast_task else "")
                + "）"
            )
        return ok(
            f"blastn 完成，返回 {len(hits)} 条命中{note}",
            query_length=len(seq),
            database=database,
            blast_task=blast_task or None,
            engine=engine,
            remote=use_remote_flag,
            hits=hits,
            count=len(hits),
        )
