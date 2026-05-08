#!/usr/bin/env python3
"""
从 nf-core/test-datasets 各 pipeline 分支拉取极小测试文件至 test_data/。
优先使用 curl 流式下载（大文件不易整包内存超时）；无 curl 时回退 urllib 分块写入。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "test_data")

# nf-core/test-datasets 数据在分支上；下列路径已校验可访问。
TARGETS = [
    (
        "genomics",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz",
        "sample1_R1.fastq.gz",
    ),
    (
        "proteomics",
        "https://raw.githubusercontent.com/nf-core/test-datasets/quantms/testdata/lfq_ci/BSA/BSA1_F1.mzML",
        "BSA1_F1.mzML",
    ),
    (
        "epigenomics",
        "https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/testdata/SRR1822153_1.fastq.gz",
        "SRR1822153_1.fastq.gz",
    ),
]


def _download_curl(url: str, out_path: str) -> None:
    subprocess.run(
        [
            "curl",
            "-fsSL",
            "--connect-timeout",
            "30",
            "--max-time",
            "900",
            "-o",
            out_path,
            url,
        ],
        check=True,
    )


def _download_urllib(url: str, out_path: str) -> None:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GIBH-AGENT-V2/download_omics_test_data"},
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        with open(out_path, "wb") as f:
            shutil.copyfileobj(resp, f, length=256 * 1024)


def main() -> int:
    for sub in ("genomics", "proteomics", "epigenomics"):
        os.makedirs(os.path.join(BASE, sub), mode=0o755, exist_ok=True)

    use_curl = shutil.which("curl") is not None

    for modality, url, fname in TARGETS:
        out_path = os.path.join(BASE, modality, fname)
        try:
            if use_curl:
                _download_curl(url, out_path)
            else:
                _download_urllib(url, out_path)
        except Exception as e:
            print(f"FAILED {modality} <- {url}: {e}", file=sys.stderr)
            return 1
        if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
            print(f"FAILED empty: {out_path}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
