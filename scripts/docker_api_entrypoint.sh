#!/bin/sh
# API 容器入口：修复挂载目录权限；必要时自动创建 BepiPred 独立 venv（/opt/bepipred3-venv）。
set -e

if [ -d /app/uploads ]; then
  chown -R appuser:appuser /app/uploads 2>/dev/null || true
fi
if [ -d /app/results ]; then
  chown -R appuser:appuser /app/results 2>/dev/null || true
fi
if [ -d /app/data ]; then
  chown -R appuser:appuser /app/data 2>/dev/null || true
fi
# 允许 appuser 在挂载的源码树中创建 BepiPred-3.0/.venv（惰性 pip 安装）
if [ -d /app/third_party/BepiPred-3.0 ]; then
  chown -R appuser:appuser /app/third_party/BepiPred-3.0 2>/dev/null || true
fi

# BepiPred：优先使用挂载的 third_party/.venv；否则镜像内 /opt/bepipred3-venv；再不行则运行时 pip（须 /usr/local/bin/python3）
API_PY="/usr/local/bin/python3"
MOUNTED_BP_VENV="/app/third_party/BepiPred-3.0/.venv/bin/python"
BEPIPRED_BOOT="${BEPIPRED3_ENTRYPOINT_BOOTSTRAP:-1}"

_bepipred_py_ok() {
  [ -n "$1" ] && [ -x "$1" ] && "$1" -c "import torch; import esm" >/dev/null 2>&1
}

if _bepipred_py_ok "$MOUNTED_BP_VENV"; then
  echo "[entrypoint] BepiPred：使用挂载 venv $MOUNTED_BP_VENV"
elif _bepipred_py_ok "/opt/bepipred3-venv/bin/python"; then
  echo "[entrypoint] BepiPred：/opt/bepipred3-venv 已就绪"
elif [ "$BEPIPRED_BOOT" != "0" ] && [ "$BEPIPRED_BOOT" != "false" ]; then
  echo "[entrypoint] BepiPred：/opt/bepipred3-venv 缺失或无法 import torch/esm，尝试运行时安装（$API_PY）…"
  TORCH_IDX="${BEPIPRED_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"
  mkdir -p /opt/bepipred3
  if [ ! -f /opt/bepipred3/requirements.txt ] && [ -f /app/third_party/BepiPred-3.0/requirements.txt ]; then
    cp /app/third_party/BepiPred-3.0/requirements.txt /opt/bepipred3/requirements.txt
  fi
  if [ ! -f /opt/bepipred3/requirements.txt ]; then
    echo "[entrypoint] BepiPred：未找到 requirements.txt，跳过。"
  else
    rm -rf /opt/bepipred3-venv
    "$API_PY" -m venv /opt/bepipred3-venv
    /opt/bepipred3-venv/bin/pip install -U pip 'setuptools>=65,<81' wheel
    /opt/bepipred3-venv/bin/pip install torch torchvision --index-url "$TORCH_IDX"
    /opt/bepipred3-venv/bin/pip install --only-binary=:all: 'numpy==1.22.4' 'pandas==1.4.3' || true
    /opt/bepipred3-venv/bin/pip install --no-build-isolation -r /opt/bepipred3/requirements.txt
    chown -R appuser:appuser /opt/bepipred3-venv
    echo "[entrypoint] BepiPred：/opt/bepipred3-venv 已就绪。"
  fi
fi

# chem_*（RDKit）：镜像层应含 /opt/chem-rdkit-venv；须 numpy<2 + 与 Dockerfile 一致的 Draw 自检（缺 libXrender 或 NumPy2 时会失败）
if [ ! -x /opt/chem-rdkit-venv/bin/python ] || \
   ! /opt/chem-rdkit-venv/bin/python -c "import numpy as np; assert np.__version__.startswith('1.'); from rdkit.Chem import AllChem, Descriptors, Draw" >/dev/null 2>&1; then
  echo "[entrypoint] chem-rdkit：venv 缺失或 RDKit/Draw/numpy 自检失败，尝试运行时重建…"
  rm -rf /opt/chem-rdkit-venv
  "$API_PY" -m venv /opt/chem-rdkit-venv
  /opt/chem-rdkit-venv/bin/pip install -U pip
  /opt/chem-rdkit-venv/bin/pip install 'numpy>=1.24.0,<2'
  /opt/chem-rdkit-venv/bin/pip install --no-cache-dir 'rdkit-pypi>=2022.3.1'
  chown -R appuser:appuser /opt/chem-rdkit-venv
  echo "[entrypoint] chem-rdkit：/opt/chem-rdkit-venv 已就绪。"
fi

# 组学参考库：挂载 data/references 后若缺 manifest 则尝试初始化（需 bwa/bowtie2，镜像内已装）
OMICS_AUTO="${OMICS_REF_AUTO_INIT:-1}"
if [ "$OMICS_AUTO" != "0" ] && [ "$OMICS_AUTO" != "false" ]; then
  REF_DIR="${OMICS_REF_DIR:-/app/references}"
  if [ ! -f "$REF_DIR/manifest.json" ]; then
    echo "[entrypoint] 组学参考库未就绪，运行 init_omics_mock_references.py …"
    mkdir -p "$REF_DIR"
    OMICS_REF_DIR="$REF_DIR" python3 /app/scripts/init_omics_mock_references.py \
      || echo "[entrypoint] 警告：参考库初始化失败，请宿主机执行: python3 scripts/init_omics_mock_references.py"
    chown -R appuser:appuser "$REF_DIR" 2>/dev/null || true
  fi
fi

# 第三批生物医药：Biopython / promb（旧镜像未重建时在此补齐）
BIOML_SCI_BOOT="${BIOML_SCI_ENTRYPOINT_BOOTSTRAP:-1}"
if [ "$BIOML_SCI_BOOT" != "0" ] && [ "$BIOML_SCI_BOOT" != "false" ]; then
  if ! python3 -c "from Bio import AlignIO, Phylo, SeqIO" >/dev/null 2>&1; then
    echo "[entrypoint] Biopython 缺失，尝试 pip install biopython…"
    pip install --no-cache-dir 'biopython>=1.83' || \
      echo "[entrypoint] 警告：biopython 安装失败。请重建：docker compose build --no-cache api-server"
  fi
  if ! python3 -c "import promb" >/dev/null 2>&1; then
    echo "[entrypoint] promb 缺失（OASis 人源性），尝试 pip install promb…"
    pip install --no-cache-dir 'promb>=0.1.0' || \
      echo "[entrypoint] 警告：promb 安装失败。请重建 api-server 镜像。"
  fi
fi

exec gosu appuser "$@"
