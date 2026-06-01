# -*- coding: utf-8 -*-
"""序列格式转换 — BioPython SeqIO（FASTA / GenBank / 纯序列）。"""
from __future__ import annotations

from io import StringIO
from typing import Any, Dict

from gibh_agent.skills._skill_common import resolve_sequence_or_fasta
from gibh_agent.skills._skill_payload import error_result, success_payload
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SUPPORTED_IN = frozenset({"fasta", "genbank", "gb", "raw"})
_SUPPORTED_OUT = frozenset({"fasta", "genbank", "gb"})


class SeqFormatConverterSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "seq_format_converter"
    display_name = "序列格式转换工具"
    description = (
        "使用 BioPython 在 FASTA、GenBank 等核酸/蛋白序列表示之间转换，"
        "输出纯文本结构（非化学小分子格式）。"
    )
    category = "生物医药"
    sub_category = "数据处理"
    aliases = ["FASTA转换", "GenBank转换", "序列格式", "SeqIO"]
    required_parameters = ["sequence_text 或 sequence_or_path"]
    tool_chain_key = ""
    output_type = "markdown"
    __dependencies__ = ["pip:biopython"]

    def execute(
        self,
        sequence_text: str = "",
        sequence_or_path: str = "",
        input_format: str = "fasta",
        output_format: str = "fasta",
        record_id: str = "seq1",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "sequence_text": (sequence_text or "").strip(),
                "sequence_or_path": (sequence_or_path or "").strip(),
                "input_format": (input_format or "fasta").strip().lower(),
                "output_format": (output_format or "fasta").strip().lower(),
                "record_id": (record_id or "seq1").strip(),
            },
        )
        inp_fmt = str(filled.get("input_format") or "fasta").lower()
        out_fmt = str(filled.get("output_format") or "fasta").lower()
        if inp_fmt == "gb":
            inp_fmt = "genbank"
        if out_fmt == "gb":
            out_fmt = "genbank"
        if inp_fmt not in _SUPPORTED_IN:
            return error_result(f"不支持的 input_format: {inp_fmt}（可选 fasta/genbank/raw）")
        if out_fmt not in _SUPPORTED_OUT:
            return error_result(f"不支持的 output_format: {out_fmt}（可选 fasta/genbank）")

        from pathlib import Path

        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        path_raw = str(filled.get("sequence_or_path") or "").strip()
        text_raw = str(filled.get("sequence_text") or "").strip()
        raw_in = ""
        if path_raw:
            p = Path(path_raw)
            if p.is_file():
                try:
                    raw_in = p.read_text(encoding="utf-8", errors="replace")
                except OSError as exc:
                    return error_result(f"无法读取序列文件: {exc}")
            else:
                raw_in = path_raw
        else:
            raw_in = text_raw

        if not raw_in.strip():
            return error_result("请提供 sequence_text 或 sequence_or_path")

        rid = str(filled.get("record_id") or "seq1")
        try:
            if inp_fmt == "raw":
                seq, err_msg = resolve_sequence_or_fasta(sequence_text=raw_in)
                if err_msg:
                    return error_result(err_msg)
                record = SeqRecord(Seq(seq or ""), id=rid, description="converted")
            else:
                handle = StringIO(raw_in)
                record = SeqIO.read(handle, inp_fmt)
            out_handle = StringIO()
            SeqIO.write(record, out_handle, out_fmt)
            converted = out_handle.getvalue().strip()
        except Exception as exc:
            return error_result(f"序列格式转换失败：{exc}")

        if not converted:
            return error_result("转换结果为空，请检查输入格式与内容")

        md = (
            f"### 序列格式转换（{inp_fmt} → {out_fmt}）\n\n"
            f"**记录 ID**：`{getattr(record, 'id', rid)}`  \n"
            f"**长度**：{len(record.seq)} bp/aa\n\n"
            f"```text\n{converted[:8000]}\n```\n"
        )
        if len(converted) > 8000:
            md += "\n*（正文已截断至 8000 字符，完整序列见 `converted_text` 字段。）*\n"

        return success_payload(
            f"格式转换完成：{inp_fmt} → {out_fmt}",
            markdown=md,
            converted_text=converted,
            input_format=inp_fmt,
            output_format=out_fmt,
            sequence_length=len(record.seq),
        )
