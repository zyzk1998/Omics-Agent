#!/usr/bin/env python3
"""
STED-EC 白盒集成测试：脱离 Agent 前端与 LLM，直接串联三个核心 Tool 验证数据流转与算法执行。

用法（在项目根目录执行）：
    python tests/test_sted_ec_pipeline.py
或：
    python -m tests.test_sted_ec_pipeline

请先将下方 TEST_H5AD_PATH 等替换为真实测试数据路径与列名。
"""
import json
import logging
import os
import sys

# 保证从项目根运行时能导入 gibh_agent
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from gibh_agent.tools.sted_ec_tools import (
    sted_ec_plot_trajectory,
    sted_ec_moscot_trajectory,
    sted_ec_preprocess,
)

# ---------------------------------------------------------------------------
# 硬编码测试参数（请替换为真实路径与列名）
# ---------------------------------------------------------------------------
TEST_H5AD_PATH = "你的测试文件路径.h5ad"  # 替换为真实 h5ad 路径
TIME_KEY = "development_stage"  # 替换为测试数据的时间列
CELL_TYPE_KEY = "cell_type"  # 替换为测试数据的细胞类型列
OUTPUT_DIR = "./test_sted_ec_output"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    logger.info("🚀 开始 STED-EC 独立链路测试...")
    logger.info("=" * 60)

    # -----------------------------------------------------------------------
    # Step 1: 预处理
    # -----------------------------------------------------------------------
    logger.info("▶️ [Step 1/3] 执行预处理...")
    res1 = sted_ec_preprocess(
        h5ad_path=TEST_H5AD_PATH,
        time_key=TIME_KEY,
        cell_type_key=CELL_TYPE_KEY,
    )
    if res1.get("status") != "success":
        logger.error(
            "❌ Step 1 失败: %s",
            json.dumps(res1, indent=2, ensure_ascii=False),
        )
        if res1.get("traceback"):
            logger.error("Traceback:\n%s", res1["traceback"])
        sys.exit(1)
    processed_h5ad = res1.get("h5ad_path")
    logger.info("✅ Step 1 完成，h5ad_path: %s", processed_h5ad)
    logger.info("-" * 60)

    # -----------------------------------------------------------------------
    # Step 2: 最优传输轨迹推断（核心，易触发依赖错误）
    # -----------------------------------------------------------------------
    logger.info("▶️ [Step 2/3] 执行 Moscot 轨迹推断...")
    res2 = sted_ec_moscot_trajectory(
        h5ad_path=processed_h5ad,
        time_key=TIME_KEY,
        cell_type_key=CELL_TYPE_KEY,
        output_dir=OUTPUT_DIR,
    )
    if res2.get("status") != "success":
        logger.error(
            "❌ Step 2 失败: %s",
            json.dumps(res2, indent=2, ensure_ascii=False),
        )
        if res2.get("traceback"):
            logger.error("Traceback:\n%s", res2["traceback"])
        sys.exit(1)
    trajectory_h5ad = res2.get("h5ad_path")
    logger.info("✅ Step 2 完成，轨迹 h5ad_path: %s", trajectory_h5ad)
    logger.info("-" * 60)

    # -----------------------------------------------------------------------
    # Step 3: 轨迹可视化
    # -----------------------------------------------------------------------
    logger.info("▶️ [Step 3/3] 执行轨迹可视化...")
    res3 = sted_ec_plot_trajectory(
        trajectory_data_path=trajectory_h5ad,
        time_key=TIME_KEY,
        cell_type_key=CELL_TYPE_KEY,
        output_dir=OUTPUT_DIR,
    )
    if res3.get("status") != "success":
        logger.error(
            "❌ Step 3 失败: %s",
            json.dumps(res3, indent=2, ensure_ascii=False),
        )
        if res3.get("traceback"):
            logger.error("Traceback:\n%s", res3["traceback"])
        sys.exit(1)
    report_data = res3.get("report_data") or {}
    images = report_data.get("images") or []
    logger.info("✅ Step 3 完成，生成图片: %s", images)
    logger.info("=" * 60)
    logger.info(
        "✅ 测试圆满成功！报告数据: %s",
        json.dumps(report_data, indent=2, ensure_ascii=False),
    )


if __name__ == "__main__":
    main()
