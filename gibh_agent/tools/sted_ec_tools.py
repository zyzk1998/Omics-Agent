"""
STED-EC 单细胞时空轨迹推断工具（JiekaiLab / moscot）。

逻辑提取自 https://github.com/JiekaiLab/STED-EC（6-trajetory_analysis.ipynb 与 Treesing/treesing.py）。
使用 moscot.TemporalProblem 对连续时间点做最优传输，输出 transport 矩阵；轨迹图基于 UMAP + 时间/细胞类型着色。

全量分析：本流程对传入的 h5ad 做全量计算，无降采样；大文件场景依赖 sparse 矩阵与循环内显式 gc 防 OOM。
"""
import gc
import logging
import os
# 根治 Docker 环境下 Numba 编译缓存无权限/找不到路径的问题
os.environ["NUMBA_CACHE_DIR"] = "/tmp/numba_cache"
import shutil
import traceback
from pathlib import Path
from typing import Dict, Any, Optional, List

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)

# STED-EC 依赖说明：请使用 moscot>=0.3.5 + ott-jax>=0.4.6（与 pandas>=2、jax 0.4.x 兼容）；旧版 ott-jax 0.3.x 会触发 register_dataclass 错误


def _log_memory_mb(label: str) -> None:
    """可选：打印当前进程内存占用（MB），便于监控全量 6GB h5ad 表现。无依赖时跳过。"""
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_SELF)
        # ru_maxrss 在 Linux 上为 KB
        mb = usage.ru_maxrss / 1024.0
        logger.info("[OOM 监控] %s: 当前进程 maxrss ≈ %.1f MB", label, mb)
    except Exception:
        try:
            import psutil
            proc = psutil.Process()
            mb = proc.memory_info().rss / (1024 * 1024)
            logger.info("[OOM 监控] %s: 当前进程 RSS ≈ %.1f MB", label, mb)
        except Exception:
            pass


def _df_head_markdown(df: Any, title: str, n: int = 5) -> str:
    """生成 DataFrame 前 n 行的 Markdown 摘要（tabulate 缺失时自动降级）。"""
    try:
        if df is None:
            return ""
        head = df.head(n)
        try:
            md = head.to_markdown(index=False)
        except Exception:
            # tabulate 缺失或不兼容时兜底
            md = "```\n" + head.to_string(index=False) + "\n```"
        return f"\n\n**{title}（前{n}行）**\n\n{md}\n"
    except Exception:
        return ""


@registry.register(
    name="sted_ec_data_validation",
    description="STED-EC 数据与元数据校验：加载 h5ad，校验 time_key 与 cell_type_key 是否存在并解析实际列名。",
    category="STED_EC",
    output_type="json",
)
def sted_ec_data_validation(
    h5ad_path: str,
    time_key: str = "day",
    cell_type_key: str = "cell_type",
    output_dir: Optional[str] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """仅负责加载 h5ad，校验 time_key 和 cell_type_key 是否存在；返回实际使用的列名与路径；可选生成中间态「时间点细胞数」图。"""
    try:
        import scanpy as sc
        logger.info("正在加载并校验单细胞数据: %s", h5ad_path)
        final_path = str(Path(h5ad_path).resolve())
        if not os.path.isfile(final_path):
            raise FileNotFoundError(f"文件不存在: {h5ad_path}")
        adata = sc.read_h5ad(final_path)
        actual_time_key = time_key
        actual_cell_type_key = cell_type_key
        valid_time_keys = ["development_stage", "day", "time", "age", "timepoint"]
        valid_cell_type_keys = ["cell_type", "celltype", "annotation", "clusters", "louvain", "leiden", "label"]
        if actual_time_key not in adata.obs.columns:
            found = False
            for k in valid_time_keys:
                if k in adata.obs.columns:
                    actual_time_key = k
                    logger.info("智能嗅探到时间序列列: %s", k)
                    found = True
                    break
            if not found:
                raise ValueError(
                    "数据校验失败：找不到时间序列列。请确保数据包含 "
                    f"{valid_time_keys} 中的任意一列，或通过参数明确指定。"
                )
        if actual_cell_type_key not in adata.obs.columns:
            for k in valid_cell_type_keys:
                if k in adata.obs.columns:
                    actual_cell_type_key = k
                    logger.info("智能嗅探到细胞类型列: %s", k)
                    break
            else:
                actual_cell_type_key = None
                logger.warning("未找到细胞类型列，轨迹图将仅按时间着色。")
        n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
        logger.info("数据校验通过。包含 %d 个细胞。时间列: %s, 细胞类型列: %s", n_obs, actual_time_key, actual_cell_type_key or "(无)")
        summary = f"数据校验通过，共 {n_obs} 细胞、{n_vars} 基因；时间列: {actual_time_key}，细胞类型列: {actual_cell_type_key or '(无)'}。"
        summary += _df_head_markdown(adata.obs, "Obs 数据透视摘要")
        images_list: List[Dict[str, str]] = []
        out_dir = output_dir or str(Path(final_path).parent)
        if actual_time_key in adata.obs.columns:
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                counts = adata.obs[actual_time_key].value_counts().sort_index()
                fig, ax = plt.subplots(figsize=(8, 4))
                counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="navy", alpha=0.8)
                ax.set_title("QC: Cells per time point (Data Validation)")
                ax.set_xlabel(actual_time_key)
                ax.set_ylabel("N cells")
                plt.xticks(rotation=45, ha="right")
                plt.tight_layout()
                os.makedirs(out_dir, exist_ok=True)
                fig_path = os.path.join(out_dir, "sted_ec_validation_cells_per_time.png")
                plt.savefig(fig_path, bbox_inches="tight", dpi=150)
                plt.close()
                images_list.append({"title": "各时间点细胞数", "path": fig_path})
                logger.info("已生成中间态图: %s", fig_path)
            except Exception as e:
                logger.warning("生成校验阶段中间图失败（已跳过）: %s", e)
        return {
            "status": "success",
            "message": summary,
            "h5ad_path": final_path,
            "time_key": actual_time_key,
            "cell_type_key": actual_cell_type_key,
            "n_obs": n_obs,
            "n_vars": n_vars,
            "summary": summary,
            "report_data": {"images": images_list, "download_links": [], "summary": summary} if images_list else None,
        }
    except Exception as e:
        tb_str = traceback.format_exc()
        logger.error("sted_ec_data_validation 失败: %s\n%s", e, tb_str)
        return {"status": "error", "error": str(e), "traceback": tb_str}


@registry.register(
    name="sted_ec_time_series_formatting",
    description="STED-EC 时间序列标准化：对上一步校验后的 h5ad 执行时间序列清洗与格式转换，输出标准化 h5ad。",
    category="STED_EC",
    output_type="json",
)
def sted_ec_time_series_formatting(
    h5ad_path: str,
    time_key: str = "day",
    cell_type_key: Optional[str] = None,
    output_dir: Optional[str] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """接收上一步的 h5ad_path，执行时间序列清洗和转换（三级容错：纯数字 / 正则提取 / 自然排序编码），写出标准化 h5ad；返回新路径。"""
    try:
        import re
        import scanpy as sc
        import numpy as np
        import pandas as pd

        logger.info("执行时间序列标准化: %s", h5ad_path)
        path = str(Path(h5ad_path).resolve())
        if not os.path.isfile(path):
            return {"status": "error", "error": f"文件不存在: {path}"}
        adata = sc.read_h5ad(path)
        if time_key not in adata.obs.columns:
            return {"status": "error", "error": f"obs 中缺少时间列: {time_key}"}

        col = adata.obs[time_key].copy()
        n_total = len(col)

        def _extract_first_number(s: str) -> float:
            """从字符串中提取第一个浮点数或整数。"""
            s = str(s).strip()
            m = re.search(r"[-+]?\d*\.?\d+|\d+", s)
            if m:
                try:
                    return float(m.group())
                except ValueError:
                    return np.nan
            return np.nan

        # 防线 1：已是数值型 -> 转 float 放行
        if pd.api.types.is_numeric_dtype(col):
            time_numeric = col.astype(float)
            logger.info("时间列已是数值型，直接转为 float。")
        else:
            # 防线 2：正则提取
            time_numeric = col.astype(str).apply(_extract_first_number)
            nan_count = time_numeric.isna().sum()
            # 防线 3：若 NaN 超过一半，视为纯文本，改用 Categorical codes
            if nan_count > n_total / 2:
                logger.info("正则提取 NaN 过多(>50%%)，改用自然排序编码。")
                raw_str = col.astype(str).str.strip()
                uniq = np.sort(raw_str.unique())
                cat = pd.Categorical(raw_str, categories=uniq)
                time_numeric = pd.Series(cat.codes.astype(float), index=col.index)
            else:
                # 善后：极少数 NaN 用最大值 + 1.0 填充
                if nan_count > 0:
                    fill_val = float(time_numeric.max()) + 1.0
                    time_numeric = time_numeric.fillna(fill_val)
                    logger.info("已用 max+1 填充 %d 个正则提取失败的 NaN。", nan_count)

        adata.obs[time_key] = time_numeric
        # 去掉空行（若仍有异常值可在此过滤，当前 time_numeric 已无 NaN）
        adata = adata[~adata.obs[time_key].isna()].copy()
        n_obs = int(adata.n_obs)

        out_dir = output_dir or str(Path(path).parent)
        os.makedirs(out_dir, exist_ok=True)
        base_name = Path(path).stem
        formatted_path = str(Path(out_dir) / f"{base_name}_formatted.h5ad")
        adata.write(formatted_path)
        logger.info("已写出标准化 h5ad: %s，保留 %d 细胞", formatted_path, n_obs)
        summary = f"时间序列标准化完成，保留 {n_obs} 细胞，已写出 {formatted_path}。"
        summary += _df_head_markdown(adata.obs, "标准化后 Obs 数据透视摘要")
        images_list: List[Dict[str, str]] = []
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            counts = adata.obs[time_key].value_counts().sort_index()
            fig, ax = plt.subplots(figsize=(8, 4))
            counts.plot(kind="bar", ax=ax, color="teal", edgecolor="darkgreen", alpha=0.8)
            ax.set_title("Cells per time point after formatting (Time Series)")
            ax.set_xlabel(time_key)
            ax.set_ylabel("N cells")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            fig_path = os.path.join(out_dir, "sted_ec_format_cells_per_time.png")
            plt.savefig(fig_path, bbox_inches="tight", dpi=150)
            plt.close()
            images_list.append({"title": "标准化后各时间点细胞数", "path": fig_path})
            logger.info("已生成中间态图: %s", fig_path)
        except Exception as e:
            logger.warning("生成时序标准化中间图失败（已跳过）: %s", e)
        del adata
        gc.collect()
        return {
            "status": "success",
            "message": summary,
            "h5ad_path": formatted_path,
            "n_obs": n_obs,
            "summary": summary,
            "report_data": {"images": images_list, "download_links": [], "summary": summary} if images_list else None,
        }
    except Exception as e:
        tb_str = traceback.format_exc()
        logger.error("sted_ec_time_series_formatting 失败: %s\n%s", e, tb_str)
        return {"status": "error", "error": str(e), "traceback": tb_str}


@registry.register(
    name="sted_ec_preprocess",
    description="STED-EC 预处理：对 h5ad 执行基础过滤与时间序列 Metadata 校验，为轨迹推断做准备。",
    category="STED_EC",
    output_type="json",
)
def sted_ec_preprocess(
    h5ad_path: str,
    time_key: str = "day",
    cell_type_key: str = "cell_type",
    **kwargs: Any,
) -> Dict[str, Any]:
    """执行基础过滤与时间序列 Metadata 校验；仅在严格生物学白名单内智能嗅探，绝不以 batch/condition 作时间。"""
    try:
        import scanpy as sc

        logger.info("正在加载并校验单细胞数据: %s", h5ad_path)
        final_path = str(Path(h5ad_path).resolve())
        if not os.path.isfile(final_path):
            raise FileNotFoundError(f"文件不存在: {h5ad_path}")
        adata = sc.read_h5ad(final_path)

        # 1. 优先使用用户传入的参数
        actual_time_key = time_key
        actual_cell_type_key = cell_type_key

        # 2. 严格生物学白名单：绝不可将 batch/condition/sample 当作时间序列
        valid_time_keys = ["development_stage", "day", "time", "age", "timepoint"]
        valid_cell_type_keys = ["cell_type", "celltype", "annotation", "clusters", "louvain", "leiden"]

        if actual_time_key not in adata.obs.columns:
            found = False
            for k in valid_time_keys:
                if k in adata.obs.columns:
                    actual_time_key = k
                    logger.info("智能嗅探到时间序列列: %s", k)
                    found = True
                    break
            if not found:
                raise ValueError(
                    "数据校验失败：找不到时间序列列。请确保数据包含 "
                    f"{valid_time_keys} 中的任意一列，或通过参数明确指定。绝不可使用 batch 或 condition 作为时间！"
                )

        if actual_cell_type_key not in adata.obs.columns:
            for k in valid_cell_type_keys:
                if k in adata.obs.columns:
                    actual_cell_type_key = k
                    logger.info("智能嗅探到细胞类型列: %s", k)
                    break
            else:
                actual_cell_type_key = None
                logger.warning("未找到细胞类型列，轨迹图将仅按时间着色。")

        logger.info("数据校验通过。包含 %d 个细胞。时间列: %s, 细胞类型列: %s", adata.n_obs, actual_time_key, actual_cell_type_key or "(无)")
        return {
            "status": "success",
            "h5ad_path": final_path,
            "time_key": actual_time_key,
            "cell_type_key": actual_cell_type_key,
        }
    except Exception as e:
        tb_str = traceback.format_exc()
        logger.error("🔥[sted_ec_preprocess] 崩溃:\n%s", tb_str)
        return {"status": "error", "error": str(e), "traceback": tb_str}


@registry.register(
    name="sted_ec_moscot_trajectory",
    description="基于 moscot 的最优传输轨迹推断：输入 h5ad 与时间列名，对连续时间点求解 OT，输出 transport 矩阵与带 tmaps 路径的 h5ad。",
    category="STED_EC",
    output_type="json",
)
def sted_ec_moscot_trajectory(
    h5ad_path: str,
    time_key: str = "day",
    cell_type_key: Optional[str] = None,
    n_pcs: int = 40,
    epsilon: float = 1e-3,
    tau_a: float = 0.99,
    tau_b: float = 0.999,
    scale_cost: str = "mean",
    max_iterations: int = 20,
    tag: str = "sted_ec_ott",
    output_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    对单细胞时间序列数据执行 moscot 最优传输轨迹推断（复刻 JiekaiLab/STED-EC 的 TemporalProblem 逻辑）。

    对每个连续时间对 (t1, t2)，在联合 PCA 空间上构建 TemporalProblem，求解后保存 transport 矩阵到 tmaps/。
    参数与仓库一致：epsilon, tau_a, tau_b, scale_cost, max_iterations。

    Args:
        h5ad_path: 预处理后的 h5ad 路径（需含 obs[time_key]）。
        time_key: 时间节点列名（如 'day', 'time_point'），值需可排序。
        cell_type_key: 细胞类型列名，可选（仅用于后续可视化）。
        n_pcs: 联合 PCA 维度，用于 OT 代价空间。
        epsilon: OT 熵正则（默认 1e-3）。
        tau_a: 源边际正则（默认 0.99）。
        tau_b: 目标边际正则（默认 0.999）。
        scale_cost: 代价缩放方式（"mean" / "max" 等）。
        max_iterations: 求解器最大迭代次数。
        tag: tmap 文件名前缀。
        output_dir: 输出目录；默认使用 h5ad 同目录下的 sted_ec_out。

    Returns:
        包含 status, h5ad_path（写入 uns['sted_ec_tmaps_dir']）, tmaps_dir, time_key, error 的字典。
    """
    try:
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import anndata
        from scipy.sparse import csr_matrix
        from moscot.problems.time import TemporalProblem
    except ImportError as e:
        logger.error("sted_ec_moscot_trajectory 依赖缺失: %s", e, exc_info=True)
        return {"status": "error", "error": f"缺少依赖: {e}"}
    except TypeError as e:
        if "register_dataclass" in str(e):
            msg = (
                "ott-jax 与当前 jax/flax 版本不兼容（register_dataclass API 变更）。"
                "请使用 moscot==0.3.5、ott-jax==0.4.6、jax==0.4.20、jaxlib==0.4.20 并重新构建 Docker 镜像。"
            )
            logger.error("🚨 %s 原始错误: %s", msg, e)
            return {"status": "error", "error": msg}
        raise

    try:
        h5ad_path = str(Path(h5ad_path).resolve())
        if not os.path.isfile(h5ad_path):
            return {"status": "error", "error": f"文件不存在: {h5ad_path}"}

        adata = sc.read_h5ad(h5ad_path)
        if time_key not in adata.obs.columns:
            return {"status": "error", "error": f"obs 中缺少时间列: {time_key}"}

        # 按数值排序时间点，避免字符串排序导致 19 排在 7 前（7.0 -> 7.5 -> ... -> 19.0）
        try:
            times = np.sort(adata.obs[time_key].astype(float).unique())
        except (ValueError, TypeError):
            times = np.sort(adata.obs[time_key].astype(str).unique())
        if len(times) < 2:
            return {"status": "error", "error": "至少需要两个时间点才能做轨迹推断"}

        out_base = output_dir or os.path.join(os.path.dirname(h5ad_path), "sted_ec_out")
        tmaps_dir = os.path.join(out_base, "tmaps")
        os.makedirs(tmaps_dir, exist_ok=True)

        col = adata.obs[time_key]
        use_numeric = np.issubdtype(times.dtype, np.floating) or times.dtype.kind == "f"
        pairs_done = 0
        for i in range(len(times) - 1):
            t1, t2 = times[i], times[i + 1]
            logger.info("sted_ec_moscot_trajectory: coupling %s -> %s", t1, t2)
            if use_numeric:
                col_f = col.astype(float)
                mask_t1 = np.abs(col_f - float(t1)) < 1e-9
                mask_t2 = np.abs(col_f - float(t2)) < 1e-9
            else:
                mask_t1 = (col.astype(str) == str(t1))
                mask_t2 = (col.astype(str) == str(t2))
            # .copy() 避免视图导致隐性引用，便于循环内释放
            adata1 = adata[mask_t1].copy()
            adata2 = adata[mask_t2].copy()
            if adata1.n_obs == 0 or adata2.n_obs == 0:
                logger.warning("时间点 %s 或 %s 无细胞，跳过", t1, t2)
                continue

            # 联合 PCA 空间（复刻 Treesing runPCA 思路：拼接后 PCA）
            merged = anndata.concat([adata1, adata2], join="inner")
            if "highly_variable" not in merged.var.columns:
                sc.pp.highly_variable_genes(merged, min_mean=0.0125, max_mean=3, min_disp=0.7, n_top_genes=min(2000, merged.n_vars))
            feats = merged.var_names[merged.var.get("highly_variable", np.ones(merged.n_vars, dtype=bool))].tolist()
            if len(feats) < 2:
                feats = merged.var_names[: min(500, merged.n_vars)].tolist()
            sub = merged[:, feats]
            sc.pp.scale(sub, max_value=10)
            sc.tl.pca(sub, n_comps=min(n_pcs, sub.n_obs - 1, sub.n_vars - 1), svd_solver="arpack")
            merged.obsm["X_pca"] = sub.obsm["X_pca"]

            day_vals = merged.obs[time_key].values
            try:
                merged.obs[time_key] = pd.Categorical(day_vals, categories=np.sort(np.unique(day_vals)))
            except Exception:
                merged.obs[time_key] = np.array(day_vals)

            tp = TemporalProblem(merged).prepare(
                time_key=time_key,
                joint_attr="X_pca",
            ).solve(
                epsilon=epsilon,
                tau_a=tau_a,
                tau_b=tau_b,
                scale_cost=scale_cost,
                min_iterations=1,
                max_iterations=max_iterations,
            )
            t1_val, t2_val = str(t1), str(t2)
            if (t1_val, t2_val) not in tp.solutions:
                keys = [k for k in tp.solutions.keys() if str(k[0]) == t1_val and str(k[1]) == t2_val]
                if not keys:
                    logger.warning("未找到解 (%s, %s)，跳过", t1_val, t2_val)
                    del merged, sub, tp
                    gc.collect()
                    continue
                key = keys[0]
            else:
                key = (t1_val, t2_val)
            transport = tp.solutions[key].transport_matrix
            if hasattr(transport, "toarray"):
                transport = transport.toarray()
            cps_adata = anndata.AnnData(
                X=csr_matrix(transport),
                obs=pd.DataFrame(index=adata1.obs_names),
                var=pd.DataFrame(index=adata2.obs_names),
            )
            out_path = os.path.join(tmaps_dir, f"{tag}_{t1_val}_{t2_val}.h5ad")
            cps_adata.write(out_path)
            logger.info("已写入 %s", out_path)
            pairs_done += 1
            # 🔥 OOM 防御：释放本时间对的大对象，避免 AnnData 视图/引用滞留
            del merged, sub, tp, cps_adata, transport
            gc.collect()
            _log_memory_mb(f"moscot 时间对 {t1_val}->{t2_val} 完成后")

        adata.uns["sted_ec_tmaps_dir"] = tmaps_dir
        adata.uns["sted_ec_time_key"] = time_key
        if cell_type_key and cell_type_key in adata.obs.columns:
            adata.uns["sted_ec_cell_type_key"] = cell_type_key
        result_path = os.path.join(out_base, "adata_with_tmaps.h5ad")
        adata.write(result_path)
        summary = f"已完成 {pairs_done} 个时间对 OT 求解，transport 矩阵已写入 tmaps；结果 h5ad 已保存。"
        # 🔥 404 终极修复：StaticFiles 不能下载目录，必须打包为 zip 再返回路径
        tmaps_zip_path = os.path.join(out_base, "tmaps.zip")
        try:
            if os.path.isdir(tmaps_dir):
                if os.path.isfile(tmaps_zip_path):
                    os.remove(tmaps_zip_path)
                shutil.make_archive(os.path.join(out_base, "tmaps"), "zip", out_base, "tmaps")
                download_links = [{"title": "tmaps (transport 矩阵) ZIP", "path": tmaps_zip_path}]
            else:
                download_links = [{"title": "tmaps (transport 矩阵)", "path": tmaps_dir}]
        except Exception as e:
            logger.warning("tmaps 打包 zip 失败，仍返回目录路径: %s", e)
            download_links = [{"title": "tmaps 目录 (transport 矩阵)", "path": tmaps_dir}]
        return {
            "status": "success",
            "message": "moscot 轨迹推断完成，transport 矩阵已写入 tmaps",
            "h5ad_path": result_path,
            "tmaps_dir": tmaps_dir,
            "time_key": time_key,
            "summary": summary,
            "report_data": {"images": [], "download_links": download_links, "summary": summary},
        }
    except Exception as e:
        logger.error("sted_ec_moscot_trajectory 执行失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="sted_ec_plot_trajectory",
    description="STED-EC 轨迹可视化：根据 h5ad 与时间/细胞类型列绘制多维 UMAP 与细胞类型演化图，保存到 sted_ec_report_images 并打包 ZIP 供下载。",
    category="STED_EC",
    output_type="mixed",
)
def sted_ec_plot_trajectory(
    trajectory_data_path: str,
    output_plot_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    plot_type: str = "umap",
    time_key: Optional[str] = None,
    cell_type_key: Optional[str] = None,
) -> Dict[str, Any]:
    """
    加载轨迹推断结果 h5ad，生成多维顶刊级图表：时空 UMAP、细胞类型 UMAP、细胞类型演化堆叠柱状图；
    全部保存到 output_dir/sted_ec_report_images，并打包为 sted_ec_results.zip 供下载。

    Args:
        trajectory_data_path: 轨迹结果 h5ad 路径（可含 uns['sted_ec_time_key'] / sted_ec_cell_type_key）。
        output_plot_path: 单张图路径（兼容旧参数，当前所有图统一写入 sted_ec_report_images）。
        output_dir: 输出目录；图表子目录与 ZIP 均在此目录下。
        plot_type: 图类型，当前支持 'umap'。
        time_key: 时间列名；若不提供则从 adata.uns['sted_ec_time_key'] 读取。
        cell_type_key: 细胞类型列名；若不提供则从 adata.uns['sted_ec_cell_type_key'] 读取。

    Returns:
        包含 status, report_data: { images: [{title, path}], download_links: [{title, path}], summary }, error 的字典。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        import scanpy as sc
        import seaborn as sns
        import numpy as np
    except ImportError as e:
        logger.error("sted_ec_plot_trajectory 依赖缺失: %s", e, exc_info=True)
        return {"status": "error", "error": str(e), "report_data": {"images": [], "download_links": [], "summary": ""}}

    try:
        # 未解析的占位符（上一步失败时）给出明确提示
        if isinstance(trajectory_data_path, str) and trajectory_data_path.strip().startswith("<") and trajectory_data_path.strip().endswith(">"):
            return {
                "status": "error",
                "error": "轨迹数据路径未解析（上一步「最优传输轨迹推断」可能未成功）。请先完成 moscot 轨迹推断步骤。",
                "report_data": {"images": [], "download_links": [], "summary": ""},
            }
        trajectory_data_path = str(Path(trajectory_data_path).resolve())
        if not os.path.isfile(trajectory_data_path):
            return {"status": "error", "error": f"文件不存在: {trajectory_data_path}", "report_data": {"images": [], "download_links": [], "summary": ""}}

        adata = sc.read_h5ad(trajectory_data_path)
        time_key = time_key or adata.uns.get("sted_ec_time_key", "day")
        # 🔥 画图步骤独立智能嗅探：前端传空字符串时不再降级，在函数内自决 cell_type_key
        if not cell_type_key or (isinstance(cell_type_key, str) and cell_type_key.strip() == ""):
            cell_type_key = adata.uns.get("sted_ec_cell_type_key")  # 上一步可能已存
        if not cell_type_key or (isinstance(cell_type_key, str) and cell_type_key.strip() == ""):
            valid_cell_type_keys = ["cell_type", "celltype", "annotation", "clusters", "louvain", "leiden", "label"]
            for k in valid_cell_type_keys:
                if k in adata.obs.columns:
                    cell_type_key = k
                    logger.info("sted_ec_plot_trajectory 智能嗅探到细胞类型列: %s", k)
                    break
        if time_key not in adata.obs.columns:
            return {"status": "error", "error": f"obs 中缺少时间列: {time_key}", "report_data": {"images": [], "download_links": [], "summary": ""}}

        if "X_umap" not in adata.obsm:
            logger.info("未检测到 X_umap，正在计算 PCA + 邻域 + UMAP")
            if "X_pca" not in adata.obsm:
                n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
                if n_comps < 2:
                    return {"status": "error", "error": "细胞或基因数过少，无法计算 PCA/UMAP", "report_data": {"images": [], "download_links": [], "summary": ""}}
                sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
            sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=min(50, adata.obsm["X_pca"].shape[1]))
            sc.tl.umap(adata)

        out_dir = output_dir or os.path.dirname(trajectory_data_path)
        if output_plot_path:
            out_dir = os.path.dirname(output_plot_path)
        img_dir = os.path.join(out_dir, "sted_ec_report_images")
        os.makedirs(img_dir, exist_ok=True)

        images_list: List[Dict[str, str]] = []

        # 图表 1：时空分布 UMAP (必定执行) —— 严格使用独立变量 path1
        sc.pl.umap(adata, color=time_key, show=False, title="Spatiotemporal UMAP (by Time)")
        path1 = os.path.join(img_dir, "umap_time.png")
        plt.savefig(path1, bbox_inches="tight", dpi=300)
        plt.close()
        images_list.append({"title": "Spatiotemporal UMAP", "path": path1})
        logger.info("Saved Spatiotemporal UMAP: %s", path1)

        # 图表 2 & 3：细胞类型相关 (动态执行，有则画，无则跳过) —— 独立使用 path2/path3，避免变量污染
        if cell_type_key and cell_type_key in adata.obs.columns:
            try:
                # 图表 2：细胞类型 UMAP
                sc.pl.umap(adata, color=cell_type_key, show=False, title="Cell Type UMAP")
                path2 = os.path.join(img_dir, "umap_celltype.png")
                plt.savefig(path2, bbox_inches="tight", dpi=300)
                plt.close()
                images_list.append({"title": "Cell Type UMAP", "path": path2})
                logger.info("Saved Cell Type UMAP: %s", path2)

                # 图表 3：演化堆叠柱状图
                cross_tab = pd.crosstab(adata.obs[time_key], adata.obs[cell_type_key], normalize="index")
                cross_tab.plot(kind="bar", stacked=True, figsize=(10, 6))
                plt.title("Cell Type Dynamics Over Time")
                plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
                plt.tight_layout()
                path3 = os.path.join(img_dir, "cell_type_dynamics.png")
                plt.savefig(path3, dpi=300)
                plt.close()
                images_list.append({"title": "Cell Type Dynamics", "path": path3})
                logger.info("Saved Cell Type Dynamics: %s", path3)
            except Exception as e:
                logger.warning("绘制细胞类型相关图表失败，已跳过: %s", e, exc_info=True)
        else:
            obs_cols = list(adata.obs.columns)
            logger.warning(
                "期望的细胞类型列名为 [%s]，但数据中只包含以下列: [%s]。已触发降级，跳过细胞类型图表绘制。",
                cell_type_key or "(未指定)",
                ", ".join(obs_cols) if obs_cols else "(空)",
            )

        # ZIP: pack sted_ec_report_images into sted_ec_results.zip under output_dir
        zip_base = os.path.join(out_dir, "sted_ec_results")
        zip_path = zip_base + ".zip"
        if os.path.isfile(zip_path):
            os.remove(zip_path)
        shutil.make_archive(zip_base, "zip", out_dir, "sted_ec_report_images")
        logger.info("Packed report archive: %s", zip_path)

        download_links = [{"title": "📦 Download Full Results (ZIP)", "path": zip_path}]
        summary = (
            "已生成时空 UMAP、细胞类型 UMAP、细胞类型演化图；所有图表已打包至 ZIP 供下载。"
            if images_list
            else "轨迹图表已生成并打包至 ZIP。"
        )
        del adata
        gc.collect()
        return {
            "status": "success",
            "message": "轨迹多维图表已保存并打包",
            "output_plot_path": path1,
            "report_data": {
                "images": images_list,
                "download_links": download_links,
                "summary": summary,
            },
        }
    except Exception as e:
        logger.error("sted_ec_plot_trajectory 执行失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e), "report_data": {"images": [], "download_links": [], "summary": ""}}
