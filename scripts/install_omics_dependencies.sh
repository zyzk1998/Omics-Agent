#!/usr/bin/env bash
#
# GIBH Omics 重型 CLI 依赖清单（Ubuntu 示例）。
# 用途：当宿主需真实运行 gibh_agent/tools 中的 bwa/samtools/fastp/gatk 等 subprocess 骨架时，
# 可按本节安装；参考基因组与索引体积巨大，请自备路径并通过环境变量注入（勿写入仓库）。
#
# 环境变量（示例）：
#   export GIBH_REF_HG38=/data/refs/hg38.fa
#   export GIBH_BOWTIE2_HG38_INDEX=/data/indexes/hg38/hg38
#
set -euo pipefail

echo "=== 1) APT：基础构建依赖与常用生信包（Ubuntu 22.04+/24.04） ==="
echo "# sudo apt-get update"
echo "# sudo apt-get install -y \\"
echo "#   build-essential git wget curl unzip \\"
echo "#   default-jdk-headless \\"
echo "#   bwa samtools bcftools \\"
echo "#   bowtie2 \\"
echo "#   fastp \\"
echo "#   picard-tools \\"
echo "#   libbz2-dev zlib1g-dev liblzma-dev \\"
echo "#   r-base r-base-dev"

echo ""
echo "=== 2) Conda/Mamba（推荐隔离安装 GATK、MACS2、部分蛋白组工具） ==="
echo "# wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh"
echo "# bash /tmp/miniforge.sh -b -p \"\$HOME/miniforge3\""
echo "# source \"\$HOME/miniforge3/etc/profile.d/conda.sh\""
echo "# conda create -n omics-cli -y bioconda::gatk4 bioconda::macs2 bioconda::htslib bioconda::samtools bioconda::bwa bioconda::fastp"
echo "# conda activate omics-cli"

echo ""
echo "=== 3) 基因组学管线常用组件 ==="
echo "- fastp：接头修剪与 QC"
echo "- BWA-MEM + samtools：比对生成 BAM（对齐 gibh omics_genomics_runner）"
echo "- Picard MarkDuplicates 或 samtools markdup：重复标记"
echo "- GATK4：BQSR、HaplotypeCaller、CNV（需已知位点资源与 Panel-of-Norms）"
echo "- bcftools：过滤、规范化"
echo "- Ensembl VEP / snpEff：注释（需缓存数据库）"

echo ""
echo "=== 4) 表观组学 ==="
echo "- Bowtie2：ATAC/ChIP 比对索引（bowtie2-build；设置 GIBH_BOWTIE2_*_INDEX）"
echo "- MACS2 / Genrich：Peak calling"
echo "- deepTools、HOMER、MEME：可选 downstream"

echo ""
echo "=== 5) 蛋白质组学 ==="
echo "- ProteoWizard msconvert（厂商 RAW→mzML，常在 Wine 容器或 Windows）"
echo "- DIA-NN / MSFragger / Comet / Percolator：搜库与重打分（发行版二进制放入 PATH）"
echo "- MaxQuant：主要为 GUI/.NET，建议在隔离 Worker 或虚拟机运行"

echo ""
echo "=== 6) Java/GATK 提示 ==="
echo "GATK4 官方分发为 gatk 包装脚本 + 预置 JAR；Conda 包 bioconda::gatk4 通常可直接 `gatk --help`。"

echo ""
echo "本脚本仅打印可复制的安装规划，不执行包管理器，以免误改系统；在目标机器上按需取消注释并运行。"
