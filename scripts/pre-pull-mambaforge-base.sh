#!/usr/bin/env bash
# worker-pyskills 已改为基于 ubuntu:22.04（不再拉 condaforge/mambaforge）。
# 本脚本预拉小体积基底，减轻 compose build 首步等待。
# 用法：bash scripts/pre-pull-mambaforge-base.sh
set -euo pipefail

IMG="${WORKER_PYSKILLS_UBUNTU_BASE:-ubuntu:22.04}"
MAX_ATTEMPTS="${PULL_MAX_ATTEMPTS:-5}"

for i in $(seq 1 "$MAX_ATTEMPTS"); do
  echo "[pre-pull] 尝试 $i/$MAX_ATTEMPTS: docker pull $IMG"
  if docker pull "$IMG"; then
    echo "[pre-pull] 完成: $IMG"
    exit 0
  fi
  echo "[pre-pull] 失败，${i}0 秒后重试…" >&2
  sleep $((i * 10))
done

echo "[pre-pull] 已达最大重试次数。" >&2
exit 1
