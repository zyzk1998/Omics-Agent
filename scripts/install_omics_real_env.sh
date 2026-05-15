#!/usr/bin/env bash
#
# install_omics_real_env.sh — 在 Ubuntu/Debian 或 Conda 中安装三大组学 Agent 依赖的真实 CLI。
#
# 人类参考序列（hg38 / GRCh38）体积巨大，**不得**提交到 Git。请下载后放到固定目录并通过环境变量注入：
#
#   export GIBH_REF_HG38=/opt/refs/GRCh38_no_alt_analysis_set.fa
#   export GIBH_BOWTIE2_HG38_INDEX=/opt/refs/bowtie2/hg38/hg38   # 前缀，不含 .1.bt2
#   export GIBH_REF_HG19=/opt/refs/hg19.fa   # 可选
#
# BWA 索引示例（在 FASTA 同目录）：
#   bwa index /opt/refs/GRCh38_no_alt_analysis_set.fa
# Bowtie2 索引示例：
#   bowtie2-build /opt/refs/GRCh38_no_alt_analysis_set.fa /opt/refs/bowtie2/hg38/hg38
#
# 用法：
#   bash scripts/install_omics_real_env.sh apt      # 需要 sudo：系统包
#   bash scripts/install_omics_real_env.sh conda    # 当前用户：bioconda 环境 omics-real
#   bash scripts/install_omics_real_env.sh pip-sci  # 当前用户：第三批科学计算技能栈（numpy/scipy/pandas/matplotlib/biopython）
#   bash scripts/install_omics_real_env.sh print    # 仅打印说明（默认）
#
set -euo pipefail

MODE="${1:-print}"

log() { printf '%s\n' "$*"; }

install_apt() {
  log "=== APT（Ubuntu/Debian）：基础 CLI ==="
  sudo apt-get update
  sudo apt-get install -y \
    build-essential curl wget git ca-certificates \
    default-jdk-headless \
    bwa samtools bcftools bowtie2 \
    fastp \
    python3 python3-pip \
    zlib1g-dev libbz2-dev liblzma-dev
  log "可选：picard-tools trimmomatic bedtools（若仓库可用）"
  sudo apt-get install -y picard-tools trimmomatic bedtools || true
  log "说明：GATK4 / MACS2 在 Debian 默认仓库可能缺失，请用 conda 段安装或使用官方 GATK 压缩包。"
}

install_conda_env() {
  log "=== Conda（bioconda）：推荐隔离安装 gatk / macs2 ==="
  if ! command -v conda >/dev/null 2>&1; then
    log "未找到 conda。可先安装 Miniforge："
    log "  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh"
    log "  bash /tmp/miniforge.sh -b -p \"\$HOME/miniforge3\""
    log "  source \"\$HOME/miniforge3/etc/profile.d/conda.sh\""
    exit 1
  fi
  conda create -n omics-real -y \
    -c conda-forge -c bioconda \
    gatk4 macs2 htslib samtools bwa bowtie2 fastp bcftools picard pysam
  log "激活环境： conda activate omics-real"
  log "验证： gatk --help && macs2 --help && bwa && bowtie2 --help"
}

install_pip_scientific_stack() {
  log "=== pip：第三批生物医药科学计算栈（numpy / scipy / pandas / matplotlib / biopython）==="
  if ! command -v python3 >/dev/null 2>&1; then
    log "未找到 python3，请先执行 apt 子命令或自行安装 Python。"
    exit 1
  fi
  python3 -m pip install --user --upgrade pip
  python3 -m pip install --user "numpy>=1.24,<2" "scipy>=1.11" "pandas>=2" "matplotlib>=3.7" "biopython>=1.83" "promb>=0.1.0"
  log "验证： python3 -c \"import numpy, scipy, pandas, matplotlib; from Bio import SeqIO; import promb; print('ok')\""
  log "可选 ANARCI（CDR 分区）： apt install hmmer && pip install anarci"
}

print_help() {
  log "未执行安装（MODE=$MODE）。子命令： apt | conda | pip-sci | print"
  log "详细占位说明亦见 scripts/install_omics_dependencies.sh"
}

case "$MODE" in
  apt) install_apt ;;
  conda) install_conda_env ;;
  pip-sci) install_pip_scientific_stack ;;
  print|*) print_help ;;
esac
