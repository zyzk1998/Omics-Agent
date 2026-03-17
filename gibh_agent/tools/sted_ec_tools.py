"""
STED-EC 单细胞时空轨迹推断工具（JiekaiLab / moscot）。

逻辑提取自 https://github.com/JiekaiLab/STED-EC（6-trajetory_analysis.ipynb 与 Treesing/treesing.py）。
使用 moscot.TemporalProblem 对连续时间点做最优传输，输出 transport 矩阵；轨迹图基于 UMAP + 时间/细胞类型着色。
"""
import logging
import os
import traceback
from pathlib import Path
from typing import Dict, Any, Optional, List

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)

# STED-EC 依赖说明：请使用 moscot>=0.3.5 + ott-jax>=0.4.6（与 pandas>=2、jax 0.4.x 兼容）；旧版 ott-jax 0.3.x 会触发 register_dataclass 错误


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
    **kwargs: Any,
) -> Dict[str, Any]:
    """仅负责加载 h5ad，校验 time_key 和 cell_type_key 是否存在；返回实际使用的列名与路径。"""
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
        logger.info("数据校验通过。包含 %d 个细胞。时间列: %s, 细胞类型列: %s", adata.n_obs, actual_time_key, actual_cell_type_key or "(无)")
        return {
            "status": "success",
            "h5ad_path": final_path,
            "time_key": actual_time_key,
            "cell_type_key": actual_cell_type_key,
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

        out_dir = output_dir or str(Path(path).parent)
        os.makedirs(out_dir, exist_ok=True)
        base_name = Path(path).stem
        formatted_path = str(Path(out_dir) / f"{base_name}_formatted.h5ad")
        adata.write(formatted_path)
        logger.info("已写出标准化 h5ad: %s", formatted_path)
        return {"status": "success", "h5ad_path": formatted_path}
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

        times = np.sort(adata.obs[time_key].astype(str).unique())
        if len(times) < 2:
            return {"status": "error", "error": "至少需要两个时间点才能做轨迹推断"}

        out_base = output_dir or os.path.join(os.path.dirname(h5ad_path), "sted_ec_out")
        tmaps_dir = os.path.join(out_base, "tmaps")
        os.makedirs(tmaps_dir, exist_ok=True)

        for i in range(len(times) - 1):
            t1, t2 = times[i], times[i + 1]
            logger.info("sted_ec_moscot_trajectory: coupling %s -> %s", t1, t2)
            mask_t1 = (adata.obs[time_key].astype(str) == t1)
            mask_t2 = (adata.obs[time_key].astype(str) == t2)
            adata1 = adata[mask_t1].copy()
            adata2 = adata[mask_t2].copy()
            if adata1.n_obs == 0 or adata2.n_obs == 0:
                logger.warning("时间点 %s 或 %s 无细胞，跳过", t1, t2)
                continue

            # 联合 PCA 空间（复刻 Treesing runPCA 思路：拼接后 PCA）
            merged = adata1.concatenate(adata2)
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
                # 兼容 key 为数值类型
                keys = [k for k in tp.solutions.keys() if str(k[0]) == t1_val and str(k[1]) == t2_val]
                if not keys:
                    logger.warning("未找到解 (%s, %s)，跳过", t1_val, t2_val)
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

        adata.uns["sted_ec_tmaps_dir"] = tmaps_dir
        adata.uns["sted_ec_time_key"] = time_key
        if cell_type_key and cell_type_key in adata.obs.columns:
            adata.uns["sted_ec_cell_type_key"] = cell_type_key
        result_path = os.path.join(out_base, "adata_with_tmaps.h5ad")
        adata.write(result_path)

        return {
            "status": "success",
            "message": "moscot 轨迹推断完成，transport 矩阵已写入 tmaps",
            "h5ad_path": result_path,
            "tmaps_dir": tmaps_dir,
            "time_key": time_key,
        }
    except Exception as e:
        logger.error("sted_ec_moscot_trajectory 执行失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="sted_ec_plot_trajectory",
    description="STED-EC 轨迹可视化：根据 h5ad 与时间/细胞类型列绘制 UMAP 轨迹图，保存到 output_dir 并返回 report_data.images。",
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
    加载轨迹推断结果 h5ad（或 sted_ec_moscot_trajectory 输出的 adata_with_tmaps.h5ad），
    在 UMAP 上按时间/细胞类型着色，复现顶刊轨迹图风格；图片保存到 output_dir 或 output_plot_path。

    Args:
        trajectory_data_path: 轨迹结果 h5ad 路径（可含 uns['sted_ec_time_key'] / sted_ec_cell_type_key）。
        output_plot_path: 单张图路径；若提供则优先使用。
        output_dir: 输出目录；若不提供 output_plot_path 则在此目录下生成 sted_ec_trajectory_umap.png 等。
        plot_type: 图类型，当前支持 'umap'。
        time_key: 时间列名；若不提供则从 adata.uns['sted_ec_time_key'] 读取。
        cell_type_key: 细胞类型列名；若不提供则从 adata.uns['sted_ec_cell_type_key'] 读取。

    Returns:
        包含 status, output_plot_path, report_data: { "images": [path, ...] }, error 的字典。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import scanpy as sc
        import numpy as np
    except ImportError as e:
        logger.error("sted_ec_plot_trajectory 依赖缺失: %s", e, exc_info=True)
        return {"status": "error", "error": str(e), "report_data": {"images": []}}

    try:
        # 未解析的占位符（上一步失败时）给出明确提示
        if isinstance(trajectory_data_path, str) and trajectory_data_path.strip().startswith("<") and trajectory_data_path.strip().endswith(">"):
            return {
                "status": "error",
                "error": f"轨迹数据路径未解析（上一步「最优传输轨迹推断」可能未成功）。请先完成 moscot 轨迹推断步骤。",
                "report_data": {"images": []},
            }
        trajectory_data_path = str(Path(trajectory_data_path).resolve())
        if not os.path.isfile(trajectory_data_path):
            return {"status": "error", "error": f"文件不存在: {trajectory_data_path}", "report_data": {"images": []}}

        adata = sc.read_h5ad(trajectory_data_path)
        time_key = time_key or adata.uns.get("sted_ec_time_key", "day")
        cell_type_key = cell_type_key or adata.uns.get("sted_ec_cell_type_key")
        if time_key not in adata.obs.columns:
            return {"status": "error", "error": f"obs 中缺少时间列: {time_key}", "report_data": {"images": []}}

        if "X_umap" not in adata.obsm:
            logger.info("未检测到 X_umap，正在计算 PCA + 邻域 + UMAP")
            if "X_pca" not in adata.obsm:
                n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
                if n_comps < 2:
                    return {"status": "error", "error": "细胞或基因数过少，无法计算 PCA/UMAP", "report_data": {"images": []}}
                sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
            sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=min(50, adata.obsm["X_pca"].shape[1]))
            sc.tl.umap(adata)

        out_dir = output_dir or os.path.dirname(trajectory_data_path)
        if output_plot_path:
            out_dir = os.path.dirname(output_plot_path)
        os.makedirs(out_dir, exist_ok=True)
        images: List[str] = []

        # 1) UMAP 按时间着色（轨迹主图）
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        sc.pl.umap(adata, color=time_key, ax=ax, show=False, title=f"Trajectory by {time_key}")
        path_time = output_plot_path or os.path.join(out_dir, "sted_ec_trajectory_umap.png")
        plt.tight_layout()
        plt.savefig(path_time, dpi=150, bbox_inches="tight")
        plt.close()
        images.append(path_time)
        logger.info("已保存轨迹 UMAP 图: %s", path_time)

        # 2) 若有细胞类型列，再画一张按细胞类型着色的 UMAP
        if cell_type_key and cell_type_key in adata.obs.columns:
            fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
            sc.pl.umap(adata, color=cell_type_key, ax=ax2, show=False, title=f"Trajectory by {cell_type_key}")
            path_cell = os.path.join(out_dir, "sted_ec_trajectory_umap_celltype.png")
            plt.tight_layout()
            plt.savefig(path_cell, dpi=150, bbox_inches="tight")
            plt.close()
            images.append(path_cell)
            logger.info("已保存细胞类型 UMAP 图: %s", path_cell)

        return {
            "status": "success",
            "message": "轨迹 UMAP 图已保存",
            "output_plot_path": path_time,
            "report_data": {"images": images},
        }
    except Exception as e:
        logger.error("sted_ec_plot_trajectory 执行失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e), "report_data": {"images": []}}
