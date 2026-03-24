#!/usr/bin/env bash
# 清理 STED-EC 白盒脚本（test_sted_ec_pipeline_*.py）产生的可再生产物，释放磁盘。
# 默认保留：test_data/sted-ec.h5ad（全量输入）、Cell Ranger 参考基因组 refdata-*（若存在）。
#
# 用法：
#   ./scripts/cleanup_sted_ec_test_artifacts.sh              # 仅当前用户可删文件
#   sudo ./scripts/cleanup_sted_ec_test_artifacts.sh         # 一并删掉容器 root 写入的目录
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TD="${ROOT}/test_data"
RF="${ROOT}/results_full"

echo "==> 清空全量测试输出: ${RF}"
mkdir -p "${RF}"
rm -rf "${RF:?}"/*

echo "==> 删除样本测试中间件（可由 sted-ec.h5ad 再生）"
rm -f \
  "${TD}/mini_sted_ec.h5ad" \
  "${TD}/mini_sted_ec_formatted.h5ad" \
  "${TD}/adata_with_tmaps.h5ad" \
  "${TD}/mini_sted_ec.diagnosis.json" \
  "${TD}/sted_ec_results.zip" \
  "${TD}/sted_ec_trajectory_umap.png" \
  "${TD}/umap_time.png" \
  "${TD}/sted_ec_validation_cells_per_time.png" \
  "${TD}/sted_ec_format_cells_per_time.png" \
  2>/dev/null || true

rm -rf "${TD}/cache" 2>/dev/null || true

# 以下目录常由 Docker 以 root 创建；无 sudo 时可能报 Permission denied
for d in "${TD}/sted_ec_report_images" "${TD}/tmaps"; do
  if [ -d "$d" ]; then
    rm -rf "$d" || echo "WARN: 无法删除 $d（请用 sudo 重跑本脚本）"
  fi
done

echo "==> 可选：删除 10x PBMC 1k 演示数据（约 10GB+，需手动取消注释）"
# rm -rf "${TD}/pbmc_1k_v3_fastqs" "${TD}/pbmc_1k_v3_fastqs.tar" \
#   "${TD}/pbmc_1k_v3_output" "${TD}/pbmc_1k_v3_filtered.h5ad"

echo "完成。剩余大块一般是 sted-ec.h5ad 与 refdata-gex-*。"
echo "提示：工作流测试输出在 data/results/run_*，可用 ./scripts/cleanup_data_results.sh 按次数清理。"
du -sh "${ROOT}/results_full" "${TD}" 2>/dev/null || true
