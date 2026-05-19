#!/usr/bin/env python3
"""
在 OMICS_REF_DIR（默认 <repo>/data/references 或 /app/references）制备三大模态 E2E 用微缩参考库。

- 基因组：hg38 微缩 chr21 FASTA + bwa 索引 + samtools fai + 序列字典
- 表观：Bowtie2 索引前缀（与 hg38 同一 FASTA）
- 蛋白组：最小 UniProt 样例 FASTA

用法（仓库根目录）:
  python3 scripts/init_omics_mock_references.py
  OMICS_REF_DIR=/app/references python3 scripts/init_omics_mock_references.py
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

CHR21_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
HOST_PREPARED = Path("/tmp/omics_test_data/ref/chr21.fa")
TESTDATA_GENOME = _ROOT / "test_data/refdata-gex-GRCh38-2024-A/fasta/genome.fa"


def _ref_root() -> Path:
    raw = (os.environ.get("OMICS_REF_DIR") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    return (_ROOT / "data" / "references").resolve()


def _run(cmd: List[str], *, timeout: int = 7200) -> None:
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, timeout=timeout)


def _which(name: str) -> Optional[str]:
    return shutil.which(name)


def _download_chr21(dest_fa: Path) -> None:
    dest_fa.parent.mkdir(parents=True, exist_ok=True)
    gz_path = dest_fa.with_suffix(dest_fa.suffix + ".gz")
    print(f"[download] {CHR21_URL} -> {gz_path}")
    urllib.request.urlretrieve(CHR21_URL, gz_path)  # noqa: S310
    with gzip.open(gz_path, "rb") as gf, open(dest_fa, "wb") as out:
        shutil.copyfileobj(gf, out)
    gz_path.unlink(missing_ok=True)


def _ensure_hg38_fasta(genomics_dir: Path) -> Path:
    hg38 = genomics_dir / "hg38.fa"
    if hg38.is_file() and hg38.stat().st_size > 1000:
        return hg38
    genomics_dir.mkdir(parents=True, exist_ok=True)
    if HOST_PREPARED.is_file():
        print(f"[copy] {HOST_PREPARED} -> {hg38}")
        shutil.copy2(HOST_PREPARED, hg38)
        return hg38
    if TESTDATA_GENOME.is_file() and TESTDATA_GENOME.stat().st_size < 80_000_000:
        print(f"[copy] {TESTDATA_GENOME} -> {hg38} (small genome.fa)")
        shutil.copy2(TESTDATA_GENOME, hg38)
        return hg38
    _download_chr21(hg38)
    return hg38


def _bwa_index(fa: Path) -> None:
    bwa = _which("bwa")
    if not bwa:
        raise RuntimeError("未找到 bwa，无法建立索引")
    if not (fa.parent / f"{fa.name}.bwt").is_file():
        _run([bwa, "index", str(fa)])


def _samtools_faidx(fa: Path) -> None:
    st = _which("samtools")
    if not st:
        return
    fai = Path(f"{fa}.fai")
    if not fai.is_file():
        _run([st, "faidx", str(fa)])


def _sequence_dict(fa: Path) -> None:
    dct = fa.with_suffix(".dict")
    if dct.is_file():
        return
    gatk = _which("gatk")
    if gatk:
        try:
            _run([gatk, "CreateSequenceDictionary", "-R", str(fa), "-O", str(dct)])
            return
        except subprocess.CalledProcessError:
            pass
    # 最小合法 dict（Picard 格式简化）
    dct.write_text(
        "@HD\tVN:1.0\tSO:unsorted\n"
        "@SQ\tSN:chr21\tLN:46709983\n",
        encoding="utf-8",
    )


def _bowtie2_index(fa: Path, prefix: Path) -> None:
    bt2 = _which("bowtie2-build")
    if not bt2:
        raise RuntimeError("未找到 bowtie2-build")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if not (prefix.parent / f"{prefix.name}.1.bt2").is_file():
        _run([bt2, str(fa), str(prefix)])


def _minimal_uniprot(fasta_path: Path) -> None:
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    if fasta_path.is_file():
        return
    fasta_path.write_text(
        ">sp|P02769|ALBU_HUMAN Albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"
        "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNPQLPRPPKRPQPPPLKPTPAPQKPTPAPQKPTPAPQKPTPAPQKPTPAPQKPTPAPQKPTPAPQK\n"
        ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens GN=TP53 PE=1 SV=4\n"
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLKIRGRKRFEMRRELPSGVEEGERRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLKIRGRKRFEMRRELPSGVEEGERRRTEEEN\n",
        encoding="utf-8",
    )


def build_all(*, skip_indexes: bool = False) -> Dict[str, Any]:
    root = _ref_root()
    genomics_dir = root / "genomics"
    epi_bt2_dir = root / "epigenomics" / "bowtie2" / "hg38"
    bt2_prefix = epi_bt2_dir / "hg38"
    prot_fa = root / "proteomics" / "uniprot_mini.fasta"

    print(f"OMICS_REF_DIR = {root}")
    hg38 = _ensure_hg38_fasta(genomics_dir)

    if not skip_indexes:
        _bwa_index(hg38)
        _samtools_faidx(hg38)
        _sequence_dict(hg38)
        _bowtie2_index(hg38, bt2_prefix)

    _minimal_uniprot(prot_fa)

    manifest: Dict[str, Any] = {
        "version": 1,
        "omics_ref_dir": str(root),
        "genomics": {
            "hg38_fasta": str(hg38),
            "bwa_index": str(hg38.parent / f"{hg38.name}.bwt"),
        },
        "epigenomics": {
            "bowtie2_hg38_prefix": str(bt2_prefix),
        },
        "proteomics": {
            "search_fasta": str(prot_fa),
        },
        "env_exports": {
            "OMICS_REF_DIR": str(root),
            "GIBH_REF_HG38": str(hg38),
            "GIBH_BOWTIE2_HG38_INDEX": str(bt2_prefix),
            "GIBH_PROTEOMICS_FASTA": str(prot_fa),
        },
    }
    manifest_path = root / "manifest.json"
    root.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[ok] manifest -> {manifest_path}")
    return manifest


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="制备组学参考库")
    parser.add_argument("--skip-indexes", action="store_true", help="仅复制 FASTA，不跑 bwa/bowtie2（已索引时）")
    args = parser.parse_args()
    try:
        build_all(skip_indexes=args.skip_indexes)
    except Exception as exc:  # noqa: BLE001
        print(f"[error] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
