#!/usr/bin/env python3
"""
doctor_omics_env.py — 组学宿主环境自检与 DevOps 修补指引。

用法（仓库根目录）：
  python3 scripts/doctor_omics_env.py
  python3 scripts/doctor_omics_env.py --json

退出码：0 = 全部核心 CLI 可解析；1 = 存在缺失。
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.tools.omics_pipeline_env import (  # noqa: E402
    _detect_container_runtime,
    resolve_cli_exe,
    resolve_reference_fasta,
)

CORE_CLIS = (
    "fastp",
    "bwa",
    "samtools",
    "bcftools",
    "bowtie2",
    "gatk",
    "macs2",
)


def _dockerfile_patch_snippet() -> str:
    return """\
# --- 追加到 services/api/Dockerfile（在 apt-get 段或单独 RUN）---
RUN apt-get update && apt-get install -y --no-install-recommends \\
    bwa samtools bcftools bowtie2 fastp \\
    default-jdk-headless \\
    && rm -rf /var/lib/apt/lists/*

# 可选：Miniconda + bioconda 环境（体积较大，适合需要 gatk4/macs2 的场景）
# ARG MINIFORGE_INSTALL_DIR=/opt/conda
# RUN wget -qO /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \\
#     && bash /tmp/miniforge.sh -b -p ${MINIFORGE_INSTALL_DIR} \\
#     && ${MINIFORGE_INSTALL_DIR}/bin/conda create -y -n omics-real -c conda-forge -c bioconda \\
#        gatk4 macs2 htslib samtools bwa bowtie2 fastp bcftools \\
#     && rm /tmp/miniforge.sh
# ENV PATH="${MINIFORGE_INSTALL_DIR}/envs/omics-real/bin:${MINIFORGE_INSTALL_DIR}/bin:${PATH}"
"""


def _compose_mount_snippet() -> str:
    home = os.path.expanduser("~")
    return f"""\
# --- docker-compose.override.yml（本地，已 gitignore）示例：挂载宿主机 Conda bin ---
services:
  api:
    environment:
      - PATH={home}/miniforge3/envs/omics-real/bin:/usr/local/bin:/usr/bin:/bin
      - GIBH_REF_HG38=/refs/GRCh38.fa
    volumes:
      - {home}/miniforge3/envs/omics-real/bin:/opt/host-omics-bin:ro
    # 若使用挂载目录，可将 PATH 首项改为 /opt/host-omics-bin
"""


def _bare_metal_fix_snippet() -> str:
    return """\
# --- 宿主机裸机 ---
bash scripts/install_omics_real_env.sh apt    # 系统包：fastp bwa samtools ...
bash scripts/install_omics_real_env.sh conda  # 或 Conda 环境 omics-real

# 启动 API 前激活环境（示例）
source ~/miniforge3/etc/profile.d/conda.sh
conda activate omics-real
export PATH="$CONDA_PREFIX/bin:$PATH"
python server.py

# systemd 示例（Unit 内 Environment）
# Environment=PATH=/home/<user>/miniforge3/envs/omics-real/bin:/usr/local/bin:/usr/bin:/bin
"""


def run_doctor(*, as_json: bool = False) -> int:
    is_container, container_detail = _detect_container_runtime()
    rows: List[Dict[str, Any]] = []
    missing: List[str] = []

    for cli in CORE_CLIS:
        resolved = resolve_cli_exe(cli)
        row = {
            "cli": cli,
            "resolved": resolved,
            "ok": bool(resolved),
        }
        rows.append(row)
        if not resolved:
            missing.append(cli)

    ref_hg38 = resolve_reference_fasta("hg38")
    report: Dict[str, Any] = {
        "is_container": is_container,
        "container_detail": container_detail,
        "user": os.environ.get("USER") or os.environ.get("LOGNAME"),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "path": os.environ.get("PATH", ""),
        "tools": rows,
        "missing": missing,
        "gibh_ref_hg38_resolved": ref_hg38,
    }

    if as_json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 1 if missing else 0

    print("=" * 60)
    print("GIBH 组学环境 Doctor")
    print("=" * 60)
    print(f"运行环境: {'容器内' if is_container else '宿主机裸机'} ({container_detail})")
    print(f"用户: {report['user'] or 'Unknown'}")
    print(f"CONDA_PREFIX: {report['conda_prefix'] or '(unset)'}")
    print(f"PATH: {report['path'][:200]}{'…' if len(report['path']) > 200 else ''}")
    print()
    print("| CLI | resolve_cli_exe |")
    print("|-----|-----------------|")
    for r in rows:
        status = r["resolved"] or "**MISSING**"
        print(f"| {r['cli']} | `{status}` |")
    print()
    print(f"参考 hg38 (resolve_reference_fasta): {ref_hg38 or '**(未配置 GIBH_REF_HG38)**'}")
    print()

    if missing:
        print(f"❌ 缺失 {len(missing)} 项: {', '.join(missing)}")
    else:
        print("✅ 核心 CLI 均已通过 resolve_cli_exe 解析。")

    print()
    print("-" * 60)
    if is_container:
        print("【容器内修复 — 修改 Dockerfile 并重建镜像】")
        print(_dockerfile_patch_snippet())
        print("【或挂载宿主机工具链（docker-compose.override.yml）】")
        print(_compose_mount_snippet())
        print("重建示例：")
        print("  cd /path/to/GIBH-AGENT-V2")
        print("  docker compose build api --no-cache")
        print("  docker compose up -d api")
    else:
        print("【裸机修复 — Conda / apt + 启动时注入 PATH】")
        print(_bare_metal_fix_snippet())

    print("-" * 60)
    print("完整安装说明: scripts/install_omics_real_env.sh")
    return 1 if missing else 0


def main() -> None:
    parser = argparse.ArgumentParser(description="组学宿主环境自检")
    parser.add_argument("--json", action="store_true", help="JSON 输出")
    args = parser.parse_args()
    raise SystemExit(run_doctor(as_json=args.json))


if __name__ == "__main__":
    main()
