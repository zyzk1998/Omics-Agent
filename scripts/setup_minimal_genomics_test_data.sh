#!/usr/bin/env bash
# 一键制备：微缩 hg38 chr21 参考 + NCBI 小型真实 FASTQ（≈5 万条 read），供 Web UI 拖拽端到端联调。
set -euo pipefail

ROOT="/tmp/omics_test_data"
REF_DIR="${ROOT}/ref"
FASTQ_DIR="${ROOT}/fastq"
CHR21_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
CHR21_GZ="${REF_DIR}/chr21.fa.gz"
CHR21_FA="${REF_DIR}/chr21.fa"
CHR21_DICT="${REF_DIR}/chr21.dict"
SRA_TARBALL="/tmp/sratoolkit.tar.gz"
SRA_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz"
SRA_DIR="/tmp/sratoolkit.3.1.1-ubuntu64"
PICARD_JAR="${PICARD_JAR:-/tmp/picard-2.27.5-for-dict.jar}"
# Picard 3.x 需 Java 17+；使用 2.27.5 兼容常见 JDK8/11。
PICARD_DL_URL="https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar"

mkdir -p "${REF_DIR}" "${FASTQ_DIR}"

run_create_sequence_dictionary() {
  local fa="$1"
  local dict="$2"
  if command -v gatk >/dev/null 2>&1; then
    gatk CreateSequenceDictionary -R "${fa}" -O "${dict}"
    return 0
  fi
  if command -v picard >/dev/null 2>&1; then
    picard CreateSequenceDictionary "R=${fa}" "O=${dict}"
    return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    for env in omics-real omics gibh; do
      if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${env}"; then
        if conda run -n "${env}" --no-capture-output gatk --help >/dev/null 2>&1; then
          conda run -n "${env}" --no-capture-output gatk CreateSequenceDictionary -R "${fa}" -O "${dict}"
          return 0
        fi
      fi
    done
  fi
  command -v java >/dev/null 2>&1 || {
    echo "ERROR: 需要 java、gatk 或 picard 之一以生成序列字典" >&2
    exit 1
  }
  if [[ ! -f "${PICARD_JAR}" ]]; then
    echo "下载 Picard JAR（用于 CreateSequenceDictionary）..."
    wget -qO "${PICARD_JAR}" "${PICARD_DL_URL}"
  fi
  java -jar "${PICARD_JAR}" CreateSequenceDictionary "R=${fa}" "O=${dict}"
}

echo "[1/4] 下载 chr21 参考序列..."
if [[ ! -s "${CHR21_FA}" ]]; then
  wget -qO "${CHR21_GZ}" "${CHR21_URL}"
  gzip -dc "${CHR21_GZ}" >"${CHR21_FA}"
fi

echo "[2/4] 建立索引（bwa / samtools / 序列字典）..."
command -v bwa >/dev/null 2>&1 || {
  echo "ERROR: 未找到 bwa，请先安装或加载 env（见 scripts/install_omics_real_env.sh）" >&2
  exit 1
}
command -v samtools >/dev/null 2>&1 || {
  echo "ERROR: 未找到 samtools" >&2
  exit 1
}

if [[ ! -f "${CHR21_FA}.bwt" ]]; then
  bwa index "${CHR21_FA}"
fi
if [[ ! -f "${CHR21_FA}.fai" ]]; then
  samtools faidx "${CHR21_FA}"
fi
if [[ ! -s "${CHR21_DICT}" ]]; then
  run_create_sequence_dictionary "${CHR21_FA}" "${CHR21_DICT}"
fi

echo "[3/4] 下载 SRA Toolkit 并截取测试 FASTQ..."
if [[ ! -x "${SRA_DIR}/bin/fastq-dump" ]]; then
  wget -qO "${SRA_TARBALL}" "${SRA_URL}"
  tar -xzf "${SRA_TARBALL}" -C /tmp
fi
"${SRA_DIR}/bin/fastq-dump" -X 50000 --split-files --gzip SRR622461 -O "${FASTQ_DIR}/"

echo "[4/4] 写入 README..."
cat >"${ROOT}/README.txt" <<'EOF'
omics_test_data — 基因组 UI 联调微缩物料
========================================

ref/
  chr21.fa          微缩参考（hg38 chr21，与 UCSC goldenPath 一致）
  chr21.fa.gz       原始下载（可删）
  chr21.fa.*        bwa 索引
  chr21.fa.fai      samtools faidx
  chr21.dict        GATK/Picard 序列字典

fastq/
  SRR622461_*.fastq.gz   NCBI SRA 拉取并截取前 50000 条 read（配对端 split-files）

设置环境变量（可选，未设置时 Runner 会在文件存在时使用默认路径）:
  export GIBH_REF_HG38=/tmp/omics_test_data/ref/chr21.fa

上游对齐 / HaplotypeCaller 请使用与 chr21 一致的样本或 Expect 仅覆盖 chr21 的区域。
EOF

echo "完成。参考: ${CHR21_FA}"
ls -la "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null || ls -la "${FASTQ_DIR}"
