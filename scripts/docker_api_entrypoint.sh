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

# BepiPred：镜像构建层应已带 /opt/bepipred3-venv；若使用旧镜像或未 build 该层，则首次启动时尝试安装（需联网）
BEPIPRED_BOOT="${BEPIPRED3_ENTRYPOINT_BOOTSTRAP:-1}"
needs_bepipred_bootstrap=0
if [ "$BEPIPRED_BOOT" != "0" ] && [ "$BEPIPRED_BOOT" != "false" ]; then
  if [ ! -x /opt/bepipred3-venv/bin/python ]; then
    needs_bepipred_bootstrap=1
  elif ! /opt/bepipred3-venv/bin/python -c "import torch; import esm" >/dev/null 2>&1; then
    needs_bepipred_bootstrap=1
  fi
fi

if [ "$needs_bepipred_bootstrap" = "1" ]; then
  echo "[entrypoint] BepiPred：/opt/bepipred3-venv 缺失或无法 import torch/esm，尝试运行时安装…"
  TORCH_IDX="${BEPIPRED_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"
  mkdir -p /opt/bepipred3
  if [ ! -f /opt/bepipred3/requirements.txt ] && [ -f /app/third_party/BepiPred-3.0/requirements.txt ]; then
    cp /app/third_party/BepiPred-3.0/requirements.txt /opt/bepipred3/requirements.txt
  fi
  if [ ! -f /opt/bepipred3/requirements.txt ]; then
    echo "[entrypoint] BepiPred：未找到 requirements.txt（需镜像内 /opt/bepipred3/ 或挂载 /app/third_party/BepiPred-3.0/），跳过。"
  else
    rm -rf /opt/bepipred3-venv
    python3 -m venv /opt/bepipred3-venv
    /opt/bepipred3-venv/bin/pip install -U pip
    /opt/bepipred3-venv/bin/pip install torch torchvision --index-url "$TORCH_IDX"
    /opt/bepipred3-venv/bin/pip install -r /opt/bepipred3/requirements.txt
    chown -R appuser:appuser /opt/bepipred3-venv
    echo "[entrypoint] BepiPred：/opt/bepipred3-venv 已就绪。"
  fi
fi

# chem_*（RDKit）：镜像层应含 /opt/chem-rdkit-venv；须 numpy<2 + 与 Dockerfile 一致的 Draw 自检（缺 libXrender 或 NumPy2 时会失败）
if [ ! -x /opt/chem-rdkit-venv/bin/python ] || \
   ! /opt/chem-rdkit-venv/bin/python -c "import numpy as np; assert np.__version__.startswith('1.'); from rdkit.Chem import AllChem, Descriptors, Draw" >/dev/null 2>&1; then
  echo "[entrypoint] chem-rdkit：venv 缺失或 RDKit/Draw/numpy 自检失败，尝试运行时重建…"
  rm -rf /opt/chem-rdkit-venv
  python3 -m venv /opt/chem-rdkit-venv
  /opt/chem-rdkit-venv/bin/pip install -U pip
  /opt/chem-rdkit-venv/bin/pip install 'numpy>=1.24.0,<2'
  /opt/chem-rdkit-venv/bin/pip install --no-cache-dir 'rdkit-pypi>=2022.3.1'
  chown -R appuser:appuser /opt/chem-rdkit-venv
  echo "[entrypoint] chem-rdkit：/opt/chem-rdkit-venv 已就绪。"
fi

exec gosu appuser "$@"
