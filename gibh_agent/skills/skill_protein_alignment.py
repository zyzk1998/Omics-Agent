# -*- coding: utf-8 -*-
"""蛋白质序列 BLAST 比对 — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict, List

from gibh_agent.skills._skill_common import (
    err,
    find_cli,
    ok,
    resolve_remote_blastp_options,
    resolve_sequence_or_fasta,
    run_cli,
    write_temp_fasta,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class ProteinSequenceAlignmentSkill(BaseSkill):
    __abstractskill__ = False

    """
    蛋白质序列比对技能（NCBI BLAST blastp）。

    对单条蛋白查询序列执行 blastp，默认使用 NCBI 远程 nr 库。
    输入可为氨基酸序列文本或 FASTA 文件路径。

    参数:
        sequence_or_path: 蛋白序列或 FASTA 文件路径。
        sequence_text: 序列文本（二选一）。
        database: 数据库名，默认 nr。
        max_target_seqs: 最大命中数，默认 10。
        use_remote: 是否远程 BLAST，默认 true。
    """

    skill_id = "protein_sequence_blast"
    display_name = "蛋白质序列比对"
    description = "使用 NCBI BLAST blastp 对蛋白质序列进行同源性搜索（支持 FASTA 或序列文本）。"
    category = "生物医药"
    sub_category = "数据分析"
    aliases = ["blastp", "蛋白BLAST", "蛋白质序列比对", "BLAST蛋白"]
    required_parameters = []
    tool_chain_key = ""
    __dependencies__ = ["apt:ncbi-blast+"]

    def execute(
        self,
        sequence_or_path: str = "",
        sequence_text: str = "",
        database: str = "nr",
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
                "database": (database or "nr").strip(),
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
        use_remote = bool(filled.get("use_remote", True))
        max_target_seqs = int(filled.get("max_target_seqs") or 10)
        database, blast_task, timeout_s = resolve_remote_blastp_options(
            len(seq),
            str(filled.get("database") or "nr"),
            str(filled.get("blast_task") or ""),
            use_remote=use_remote,
        )

        blastp = find_cli("blastp")
        if not blastp:
            return err("未找到 blastp，请安装 ncbi-blast+")

        query_path = write_temp_fasta(seq, prefix="prot_query")
        outfmt = "6 qseqid sseqid pident length evalue bitscore stitle"
        cmd: List[str] = [
            blastp,
            "-query",
            query_path,
            "-db",
            (database or "nr").strip(),
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            str(max(1, min(int(max_target_seqs or 10), 50))),
        ]
        if blast_task:
            cmd.extend(["-task", blast_task])
        if use_remote:
            cmd.append("-remote")

        try:
            code, stdout, stderr = run_cli(cmd, timeout_s=timeout_s)
        finally:
            try:
                import os

                os.unlink(query_path)
            except OSError:
                pass

        if code != 0:
            return err(f"blastp 执行失败: {stderr or stdout}")

        hits: List[Dict[str, str]] = []
        for line in stdout.splitlines():
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            hits.append(
                {
                    "qseqid": cols[0],
                    "sseqid": cols[1],
                    "pident": cols[2],
                    "length": cols[3],
                    "evalue": cols[4],
                    "bitscore": cols[5],
                    "title": cols[6] if len(cols) > 6 else "",
                }
            )

        return ok(
            f"blastp 完成，返回 {len(hits)} 条命中",
            query_length=len(seq),
            hits=hits,
            count=len(hits),
        )
