#!/usr/bin/env python3
"""
STED-EC 全量分析白盒测试：使用未拆分的原始 sted-ec.h5ad，输出到 results_full/，带内存监控与 OOM 防御。

用法：
  docker compose run --rm api-server python test_sted_ec_pipeline_full.py
  或：python test_sted_ec_pipeline_full.py（宿主机需已安装依赖）
"""
import gc
import json
import os
import sys
import traceback

_ROOT = os.path.abspath(os.path.dirname(__file__) or ".")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from loguru import logger
except ImportError:
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    logger = logging.getLogger(__name__)

from gibh_agent.tools.sted_ec_tools import (
    sted_ec_data_validation,
    sted_ec_time_series_formatting,
    sted_ec_moscot_trajectory,
    sted_ec_plot_trajectory,
)

_TEST_DATA_DIR = os.environ.get("TEST_DATA_DIR", os.path.join(_ROOT, "test_data"))
# 全量脚本：使用原始全量 h5ad，不切 mini
FULL_H5AD = os.path.join(_TEST_DATA_DIR, "sted-ec.h5ad")
# 输出隔离：保存到 results_full/，避免与 sample 版覆盖
OUTPUT_DIR_FULL = os.path.join(_ROOT, "results_full")


def _memory_mb() -> float:
    """当前进程内存占用（MB），用于终端观察水位线。无 psutil 时返回 0。"""
    try:
        import psutil
        return psutil.Process().memory_info().rss / (1024 * 1024)
    except Exception:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        except Exception:
            return 0.0


def _log_memory(label: str) -> None:
    mb = _memory_mb()
    if mb > 0:
        logger.info("[内存监控] {}: 当前进程 ≈ {:.1f} MB".format(label, mb))


def main() -> None:
    try:
        os.makedirs(OUTPUT_DIR_FULL, exist_ok=True)
        logger.info("🚀 STED-EC 全量分析测试（输出目录: {}）", OUTPUT_DIR_FULL)
        logger.info("=" * 60)
        if not os.path.isfile(FULL_H5AD):
            logger.error("全量文件不存在: {}", FULL_H5AD)
            sys.exit(1)
        _log_memory("启动前")
        gc.collect()

        # ---------- Step 1: 数据与元数据校验 ----------
        logger.info("▶️ [Step 1/4] 数据与元数据校验...")
        _log_memory("Step1 前")
        res1 = sted_ec_data_validation(
            h5ad_path=FULL_H5AD,
            time_key="day",
            cell_type_key="cell_type",
            output_dir=OUTPUT_DIR_FULL,
        )
        _log_memory("Step1 后")
        gc.collect()
        if res1.get("status") != "success":
            logger.error("❌ Step 1 失败: {}", json.dumps(res1, indent=2, ensure_ascii=False))
            sys.exit(1)
        time_key = res1.get("time_key", "day")
        cell_type_key = res1.get("cell_type_key")
        logger.info("✅ Step 1 完成。time_key={}, cell_type_key={}", time_key, cell_type_key)
        logger.info("-" * 60)

        # ---------- Step 2: 时间序列标准化 ----------
        logger.info("▶️ [Step 2/4] 时间序列标准化...")
        _log_memory("Step2 前")
        res2 = sted_ec_time_series_formatting(
            h5ad_path=res1["h5ad_path"],
            time_key=time_key,
            cell_type_key=cell_type_key,
            output_dir=OUTPUT_DIR_FULL,
        )
        _log_memory("Step2 后")
        gc.collect()
        if res2.get("status") != "success":
            logger.error("❌ Step 2 失败: {}", json.dumps(res2, indent=2, ensure_ascii=False))
            sys.exit(1)
        formatted_path = res2.get("h5ad_path")
        logger.info("✅ Step 2 完成。formatted h5ad: {}", formatted_path)
        logger.info("-" * 60)

        # ---------- Step 3: Moscot 轨迹推断（全量易 OOM，工具内已有 del + gc） ----------
        logger.info("▶️ [Step 3/4] Moscot 轨迹推断（全量）...")
        _log_memory("Step3 前")
        res3 = sted_ec_moscot_trajectory(
            h5ad_path=formatted_path,
            time_key=time_key,
            cell_type_key=cell_type_key,
            output_dir=OUTPUT_DIR_FULL,
        )
        _log_memory("Step3 后")
        gc.collect()
        if res3.get("status") != "success":
            logger.error("❌ Step 3 失败: {}", json.dumps(res3, indent=2, ensure_ascii=False))
            sys.exit(1)
        trajectory_path = res3.get("h5ad_path")
        logger.info("✅ Step 3 完成。轨迹 h5ad: {}", trajectory_path)
        logger.info("-" * 60)

        # ---------- Step 4: 轨迹可视化 ----------
        logger.info("▶️ [Step 4/4] 轨迹可视化...")
        _log_memory("Step4 前")
        res4 = sted_ec_plot_trajectory(
            trajectory_data_path=trajectory_path,
            output_dir=OUTPUT_DIR_FULL,
            time_key=time_key,
            cell_type_key=cell_type_key,
        )
        _log_memory("Step4 后")
        gc.collect()
        if res4.get("status") != "success":
            logger.error("❌ Step 4 失败: {}", json.dumps(res4, indent=2, ensure_ascii=False))
            sys.exit(1)
        report_data = res4.get("report_data") or {}
        logger.info("✅ Step 4 完成。生成图片: {}", report_data.get("images", []))
        _log_memory("全部完成")
        logger.info("=" * 60)
        logger.info("✅ 全量测试完成。输出目录: {}", OUTPUT_DIR_FULL)
    except Exception as e:
        logger.error("未捕获异常: {}", e)
        logger.error("Traceback:\n{}", traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
