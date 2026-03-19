#!/usr/bin/env python3
"""
STED-EC 白盒集成测试：自动切片大文件 + 串联 4 步核心 Tool。

- 大文件：/app/test_data/sted-ec.h5ad（容器内映射）
- 微缩版：/app/test_data/mini_sted_ec.h5ad（3000 细胞，用于快速测试）

用法：
  1) 推荐（一次性容器，不拖垮常驻 api-server，避免 OOM/重启）：
     docker compose run --rm api-server python test_sted_ec_pipeline.py

  2) 在已运行的 api-server 内执行（可能与 Gunicorn 争资源导致容器重启）：
     docker compose exec api-server python test_sted_ec_pipeline.py

  3) 宿主机（需已安装依赖）：python test_sted_ec_pipeline.py
"""
import json
import os
import sys
import traceback

# 保证从项目根运行时能导入 gibh_agent
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

# ---------------------------------------------------------------------------
# 路径定义：优先用环境变量 TEST_DATA_DIR，否则用项目根下的 test_data（宿主机与容器通用）
# 若 test_data 不可写则自动使用 test_data_run 并在结束后将图片/zip 复制回 test_data
# ---------------------------------------------------------------------------
_DEFAULT_TEST_DATA = os.path.join(_ROOT, "test_data")
_TEST_DATA_DIR = os.environ.get("TEST_DATA_DIR", _DEFAULT_TEST_DATA)
_RUN_DIR = _TEST_DATA_DIR  # 实际运行写入的目录（可能与 _TEST_DATA_DIR 不同）
_FINAL_OUTPUT_DIR = _DEFAULT_TEST_DATA  # 最终希望输出到的目录


def _can_write_to(dirpath: str) -> bool:
    """检查是否能在目录下创建/写入文件。"""
    if not dirpath or not os.path.isdir(dirpath):
        return False
    probe = os.path.join(dirpath, ".write_probe")
    try:
        with open(probe, "w") as f:
            f.write("")
        os.remove(probe)
        return True
    except Exception:
        return False


def _ensure_writable_output_dir() -> None:
    """若默认 test_data 不可写，则使用 test_data_run 并准备输入数据。"""
    global _RUN_DIR, _TEST_DATA_DIR, BIG_H5AD, MINI_H5AD, OUTPUT_DIR
    tmaps_dir = os.path.join(_TEST_DATA_DIR, "tmaps")
    os.makedirs(tmaps_dir, exist_ok=True)
    if _can_write_to(tmaps_dir):
        _RUN_DIR = _TEST_DATA_DIR
        OUTPUT_DIR = _RUN_DIR
        return
    run_dir = os.path.join(_ROOT, "test_data_run")
    os.makedirs(run_dir, exist_ok=True)
    for name in ("sted-ec.h5ad", "mini_sted_ec.h5ad"):
        src = os.path.join(_DEFAULT_TEST_DATA, name)
        dst = os.path.join(run_dir, name)
        if os.path.isfile(src) and (not os.path.isfile(dst) or os.path.getmtime(src) > os.path.getmtime(dst)):
            import shutil
            logger.info("复制 %s -> %s（因 test_data 不可写）", name, run_dir)
            shutil.copy2(src, dst)
    _RUN_DIR = run_dir
    _TEST_DATA_DIR = run_dir
    BIG_H5AD = os.path.join(_RUN_DIR, "sted-ec.h5ad")
    MINI_H5AD = os.path.join(_RUN_DIR, "mini_sted_ec.h5ad")
    OUTPUT_DIR = _RUN_DIR


def _copy_results_to_test_data() -> None:
    """将本次运行生成的图片与 zip 复制到 test_data（若运行目录与 test_data 不同）。"""
    if _RUN_DIR == _FINAL_OUTPUT_DIR:
        return
    import shutil
    for name in ("sted_ec_results.zip",):
        src = os.path.join(_RUN_DIR, name)
        if os.path.isfile(src):
            dst = os.path.join(_FINAL_OUTPUT_DIR, name)
            try:
                shutil.copy2(src, dst)
                logger.info("已复制到 test_data: %s", name)
            except Exception as e:
                logger.warning("复制到 test_data 失败 %s: %s", name, e)
    img_dir = os.path.join(_RUN_DIR, "sted_ec_report_images")
    if os.path.isdir(img_dir):
        out_img = os.path.join(_FINAL_OUTPUT_DIR, "sted_ec_report_images")
        for f in os.listdir(img_dir):
            if f.endswith(".png"):
                src = os.path.join(img_dir, f)
                if os.path.isfile(src) and os.path.getsize(src) > 0:
                    try:
                        os.makedirs(out_img, exist_ok=True)
                        shutil.copy2(src, os.path.join(out_img, f))
                        logger.info("已复制到 test_data: sted_ec_report_images/%s", f)
                    except Exception as e:
                        dst_flat = os.path.join(_FINAL_OUTPUT_DIR, f)
                        try:
                            shutil.copy2(src, dst_flat)
                            logger.info("已复制到 test_data: %s", f)
                        except Exception as e2:
                            logger.warning("复制图片失败 %s: %s", f, e2)


BIG_H5AD = os.path.join(_TEST_DATA_DIR, "sted-ec.h5ad")
MINI_H5AD = os.path.join(_TEST_DATA_DIR, "mini_sted_ec.h5ad")
OUTPUT_DIR = _TEST_DATA_DIR
N_SUBSAMPLE = 3000


def ensure_mini_h5ad() -> str:
    """若微缩版不存在则从大文件切片并保存；返回用于流程的 h5ad 路径。"""
    logger.info("数据目录: %s，大文件路径: %s" % (_TEST_DATA_DIR, BIG_H5AD))
    if os.path.isfile(MINI_H5AD):
        logger.info("微缩版已存在，直接使用: %s", MINI_H5AD)
        return MINI_H5AD
    if not os.path.isfile(BIG_H5AD):
        logger.error("大文件不存在: %s（请确认文件已放入 test_data 或设置 TEST_DATA_DIR）" % (BIG_H5AD,))
        sys.exit(1)
    import scanpy as sc
    logger.info("正在加载大文件: %s", BIG_H5AD)
    adata = sc.read_h5ad(BIG_H5AD)
    logger.info("adata.obs.columns 完整列表: %s", list(adata.obs.columns))
    if adata.n_obs <= N_SUBSAMPLE:
        logger.info("细胞数 %d <= %d，不切片，直接另存为微缩版", adata.n_obs, N_SUBSAMPLE)
        adata.write(MINI_H5AD)
        return MINI_H5AD
    sc.pp.subsample(adata, n_obs=N_SUBSAMPLE, copy=False)
    os.makedirs(os.path.dirname(MINI_H5AD), exist_ok=True)
    adata.write(MINI_H5AD)
    logger.info("已保存微缩版 (n_obs=%d): %s", adata.n_obs, MINI_H5AD)
    return MINI_H5AD


def main() -> None:
    try:
        logger.info("🚀 开始 STED-EC 四步白盒集成测试...")
        logger.info("=" * 60)
        _ensure_writable_output_dir()
        logger.info("运行/输出目录: %s", _RUN_DIR)

        # ---------- 自动切片与 Metadata 嗅探 ----------
        h5ad_path = ensure_mini_h5ad()
        logger.info("-" * 60)

        # ---------- Step 1: 数据与元数据校验 ----------
        logger.info("▶️ [Step 1/4] 正在执行数据与元数据校验...")
        res1 = sted_ec_data_validation(
            h5ad_path=h5ad_path,
            time_key="day",
            cell_type_key="cell_type",
        )
        if res1.get("status") != "success":
            logger.error("❌ Step 1 失败: %s", json.dumps(res1, indent=2, ensure_ascii=False))
            if res1.get("traceback"):
                logger.error("Traceback:\n%s", res1["traceback"])
            sys.exit(1)
        time_key = res1.get("time_key", "day")
        cell_type_key = res1.get("cell_type_key")
        logger.info("✅ Step 1 完成。time_key=%s, cell_type_key=%s", time_key, cell_type_key)
        logger.info("-" * 60)

        # ---------- Step 2: 时间序列标准化 ----------
        logger.info("▶️ [Step 2/4] 正在执行时间序列标准化...")
        res2 = sted_ec_time_series_formatting(
            h5ad_path=res1["h5ad_path"],
            time_key=time_key,
            cell_type_key=cell_type_key,
            output_dir=OUTPUT_DIR,
        )
        if res2.get("status") != "success":
            logger.error("❌ Step 2 失败: %s", json.dumps(res2, indent=2, ensure_ascii=False))
            if res2.get("traceback"):
                logger.error("Traceback:\n%s", res2["traceback"])
            sys.exit(1)
        formatted_path = res2.get("h5ad_path")
        logger.info("✅ Step 2 完成。formatted h5ad: %s", formatted_path)
        logger.info("-" * 60)

        # ---------- Step 3: Moscot 轨迹推断 ----------
        logger.info("▶️ [Step 3/4] 正在执行 Moscot 轨迹推断...")
        res3 = sted_ec_moscot_trajectory(
            h5ad_path=formatted_path,
            time_key=time_key,
            cell_type_key=cell_type_key,
            output_dir=OUTPUT_DIR,
        )
        if res3.get("status") != "success":
            logger.error("❌ Step 3 失败: %s", json.dumps(res3, indent=2, ensure_ascii=False))
            if res3.get("traceback"):
                logger.error("Traceback:\n%s", res3["traceback"])
            sys.exit(1)
        trajectory_path = res3.get("h5ad_path")
        logger.info("✅ Step 3 完成。轨迹 h5ad: %s", trajectory_path)
        logger.info("-" * 60)

        # ---------- Step 4: 轨迹可视化 ----------
        logger.info("▶️ [Step 4/4] 正在执行轨迹可视化...")
        res4 = sted_ec_plot_trajectory(
            trajectory_data_path=trajectory_path,
            output_dir=OUTPUT_DIR,
            time_key=time_key,
            cell_type_key=cell_type_key,
        )
        if res4.get("status") != "success":
            logger.error("❌ Step 4 失败: %s", json.dumps(res4, indent=2, ensure_ascii=False))
            if res4.get("traceback"):
                logger.error("Traceback:\n%s", res4["traceback"])
            sys.exit(1)
        report_data = res4.get("report_data") or {}
        logger.info("✅ Step 4 完成。生成图片: %s", report_data.get("images", []))
        _copy_results_to_test_data()
        logger.info("=" * 60)
        logger.info("✅ 测试圆满成功！report_data: %s", json.dumps(report_data, indent=2, ensure_ascii=False))

    except Exception as e:
        logger.error("未捕获异常: %s", e)
        logger.error("Traceback:\n%s", traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
