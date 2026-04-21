"""
STED_EC 域工具集：单细胞时空动力学分析（JiekaiLab / moscot 实现，底层 tool_id 不变）。

逻辑提取自 https://github.com/JiekaiLab/STED-EC（6-trajetory_analysis.ipynb 与 Treesing/treesing.py）。
使用 moscot.TemporalProblem 对连续时间点做最优传输，输出 transport 矩阵；轨迹图基于 UMAP + 时间/细胞类型着色。

全量分析：本流程对传入的 h5ad 做全量计算，无降采样；大文件场景依赖 sparse 矩阵与循环内显式 gc 防 OOM。
"""
import gc
import inspect
import json
import logging
import os
# 根治 Docker 环境下 Numba 编译缓存无权限/找不到路径的问题
os.environ["NUMBA_CACHE_DIR"] = "/tmp/numba_cache"
import shutil
import traceback
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

from ..core.tool_registry import registry
from ..core.utils import sanitize_for_json

logger = logging.getLogger(__name__)

# 时间列嗅探白名单（顺序有意义）：先 day（官方 sted-ec），再 time（肾/肺演示集等）。
# 列名不在此列表时，必须在数据校验 / 工作流表单中显式传入 time_key。
# 逗号分隔额外列名：环境变量 STED_EC_EXTRA_TIME_KEYS，例如 "Phase,dpt_pseudotime"
STED_EC_VALID_TIME_KEYS: Tuple[str, ...] = (
    "day",
    "time",
    "timepoint",
    "time_point",
    "age",
    "age_days",
    "development_stage",
    "developmental_stage",
    "embryonic_period",
    "library_time",
    "capture_time",
    "sample_time",
)

STED_EC_VALID_CELL_TYPE_KEYS: Tuple[str, ...] = (
    "cell_type",
    "celltype",
    "annotation",
    "clusters",
    "louvain",
    "leiden",
    "label",
    "level1_group",
    "level3_celltypes",
    "level2_organORlineages",
    "new_level3_celltypes",
)


def _sted_ec_time_keys_ordered() -> Tuple[str, ...]:
    extra = os.environ.get("STED_EC_EXTRA_TIME_KEYS", "").strip()
    if not extra:
        return STED_EC_VALID_TIME_KEYS
    more = tuple(x.strip() for x in extra.split(",") if x.strip())
    return STED_EC_VALID_TIME_KEYS + more


def _normalize_enabled_mcps_for_tools(raw: Any) -> List[str]:
    if raw is None:
        return []
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except Exception:
            return [s.strip() for s in raw.split(",") if s.strip()]
    if not isinstance(raw, list):
        return []
    return [str(x).strip() for x in raw if x and str(x).strip()]


def _assistant_msg_to_openai_dict(msg: Any) -> Dict[str, Any]:
    out: Dict[str, Any] = {"role": "assistant", "content": getattr(msg, "content", None) or ""}
    tcs = getattr(msg, "tool_calls", None)
    if not tcs:
        return out
    tool_calls = []
    for tc in tcs:
        fn = getattr(tc, "function", None)
        tool_calls.append(
            {
                "id": getattr(tc, "id", ""),
                "type": getattr(tc, "type", "function") or "function",
                "function": {
                    "name": getattr(fn, "name", "") if fn is not None else "",
                    "arguments": (getattr(fn, "arguments", None) or "{}") if fn is not None else "{}",
                },
            }
        )
    out["tool_calls"] = tool_calls
    return out


def _extract_references_from_content(content: str) -> Tuple[str, List[Dict[str, str]]]:
    """从定界 JSON 或文末 ```json 数组解析参考文献；返回 (正文, references)。"""
    text = (content or "").strip()
    refs: List[Dict[str, str]] = []
    start_tag = "<<<REFERENCES_JSON>>>"
    end_tag = "<<<END_REFERENCES>>>"
    if start_tag in text and end_tag in text:
        a = text.find(start_tag)
        b = text.find(end_tag, a)
        chunk = text[a + len(start_tag) : b].strip()
        body = (text[:a].rstrip() + "\n" + text[b + len(end_tag) :]).strip()
        try:
            data = json.loads(chunk)
            if isinstance(data, list):
                for item in data:
                    if not isinstance(item, dict):
                        continue
                    u = item.get("url") or item.get("href") or item.get("link")
                    if not u:
                        continue
                    refs.append(
                        {
                            "title": str(item.get("title") or item.get("name") or "Reference")[:500],
                            "url": str(u)[:2048],
                        }
                    )
        except json.JSONDecodeError:
            pass
        return body, refs

    # 回退：最后一个 ```json ... ``` 代码块且解析为数组
    fence = "```json"
    if fence in text.lower():
        idx = text.lower().rfind(fence.lower())
        rest = text[idx + len(fence) :]
        end_fence = rest.find("```")
        if end_fence != -1:
            chunk = rest[:end_fence].strip()
            try:
                data = json.loads(chunk)
                if isinstance(data, list) and data and isinstance(data[0], dict):
                    for item in data:
                        u = item.get("url") or item.get("href")
                        if u:
                            refs.append(
                                {
                                    "title": str(item.get("title") or "Reference")[:500],
                                    "url": str(u)[:2048],
                                }
                            )
                    if refs:
                        body = (text[:idx].rstrip() + rest[end_fence + 3 :].lstrip()).strip()
                        return body, refs
            except json.JSONDecodeError:
                pass
    return text, []


def _sted_safe_return(payload: Dict[str, Any]) -> Dict[str, Any]:
    """保证工具返回可被 JSON / SSE 序列化，避免 NaN、numpy 标量等导致流中断。"""
    try:
        return sanitize_for_json(payload)
    except Exception as e:
        logger.error("_sted_safe_return sanitize 失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e)}


def _robust_sanitize_and_preprocess(adata: Any, time_key: str, *, logger_label: str = "sted_ec") -> None:
    """
    防弹清洗：X 中 NaN/Inf → 0（稀疏/稠密均处理 .data 或整块矩阵）；
    若 X.max()>50 视为 Raw Counts，则 normalize_total + log1p 后再 nan_to_num；
    时间列 time_key：pd.to_numeric(coerce) + Inf 清掉 + NaN 用列中位数填充；
    其余 obs 数值列做有限化（与官方 sted-ec.h5ad 兼容：log 矩阵不触发归一化）。
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.sparse import issparse

    X = adata.X
    if issparse(X):
        X = X.copy()
        if X.nnz > 0:
            fixed = np.nan_to_num(np.asarray(X.data, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
            X.data = fixed.astype(X.data.dtype, copy=False)
        adata.X = X
    else:
        adata.X = np.nan_to_num(np.asarray(X, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    try:
        X2 = adata.X
        if issparse(X2):
            x_max = float(np.max(X2.data)) if X2.nnz else 0.0
        else:
            x_max = float(np.max(np.asarray(X2))) if np.size(X2) else 0.0
        if x_max > 50:
            logger.info("[%s] X_max=%.3f > 50，执行 normalize_total + log1p", logger_label, x_max)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            X3 = adata.X
            if issparse(X3) and X3.nnz > 0:
                X3 = X3.copy()
                fixed = np.nan_to_num(np.asarray(X3.data, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
                X3.data = fixed.astype(X3.data.dtype, copy=False)
                adata.X = X3
            elif not issparse(X3):
                adata.X = np.nan_to_num(np.asarray(X3, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    except Exception as ex:
        logger.warning("[%s] 智能归一化分支异常（已跳过）: %s", logger_label, ex)

    if time_key and time_key in adata.obs.columns:
        try:
            tk = adata.obs[time_key]
            vals = pd.to_numeric(tk, errors="coerce").replace([np.inf, -np.inf], np.nan)
            med = float(np.nanmedian(vals.to_numpy(dtype=np.float64)))
            if not np.isfinite(med):
                med = 0.0
            adata.obs[time_key] = vals.fillna(med)
        except Exception as tk_ex:
            logger.warning("[%s] 时间列 %s 清洗失败（已跳过该步）: %s", logger_label, time_key, tk_ex)

    for col in list(adata.obs.columns):
        if col == time_key:
            continue
        s = adata.obs[col]
        try:
            if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_categorical_dtype(s):
                continue
            if pd.api.types.is_bool_dtype(s):
                continue
            if pd.api.types.is_numeric_dtype(s):
                vals = pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan)
                med = float(np.nanmedian(vals.to_numpy(dtype=np.float64)))
                if not np.isfinite(med):
                    med = 0.0
                adata.obs[col] = vals.fillna(med)
                continue
            if s.dtype == object:
                num = pd.to_numeric(s, errors="coerce")
                if num.notna().sum() < max(3, len(s) // 20):
                    continue
                vals = num.replace([np.inf, -np.inf], np.nan)
                med = float(np.nanmedian(vals.to_numpy(dtype=np.float64)))
                if not np.isfinite(med):
                    med = 0.0
                adata.obs[col] = vals.fillna(med)
        except Exception as col_ex:
            logger.debug("[%s] obs 列 %s 清洗跳过: %s", logger_label, col, col_ex)


def _sted_obsm_all_finite(adata: Any, key: str) -> bool:
    import numpy as np

    if key not in adata.obsm:
        return False
    v = np.asarray(adata.obsm[key])
    if v.size == 0:
        return False
    return bool(np.all(np.isfinite(v)))


def _nan_to_num_obsm(adata: Any, key: str) -> None:
    import numpy as np

    if key not in adata.obsm:
        return
    v = np.asarray(adata.obsm[key], dtype=np.float64)
    adata.obsm[key] = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)


def _sted_ec_safe_recompute_embedding(adata: Any, logger_label: str) -> None:
    """在已清洗 X 的前提下重算 PCA → neighbors → UMAP，并对 obsm 再 nan_to_num。"""
    import numpy as np
    import scanpy as sc

    n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
    if n_comps < 2:
        raise ValueError("细胞或基因数过少，无法计算 PCA/UMAP")
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    _nan_to_num_obsm(adata, "X_pca")
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=min(50, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata)
    _nan_to_num_obsm(adata, "X_umap")
    u = np.asarray(adata.obsm["X_umap"], dtype=np.float64)
    if not np.all(np.isfinite(u)):
        adata.obsm["X_umap"] = np.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0)


_STED_PLOT_PLACEHOLDER_MSG = (
    "Visualization Failed due to data distribution issues, but pipeline continues."
)


def _write_sted_plot_placeholder(path: str, message: str = _STED_PLOT_PLACEHOLDER_MSG) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11, wrap=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def _sanitize_obs_color_column(adata: Any, col: Optional[str]) -> None:
    """绘图 color 用数值列时，保证无 Inf/NaN（category 列跳过）。"""
    if not col or col not in adata.obs.columns:
        return
    import numpy as np
    import pandas as pd

    s = adata.obs[col]
    if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_categorical_dtype(s):
        return
    if not pd.api.types.is_numeric_dtype(s) and s.dtype != object:
        return
    vals = pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan)
    med = float(np.nanmedian(vals.to_numpy(dtype=np.float64)))
    if not np.isfinite(med):
        med = 0.0
    adata.obs[col] = vals.fillna(med)


def _safe_sted_pl_umap(
    adata: Any,
    color: Optional[str],
    title: str,
    path_out: str,
    images_list: List[Dict[str, str]],
    list_title: str,
    logger_label: str,
) -> bool:
    import matplotlib.pyplot as plt
    import scanpy as sc

    try:
        sc.pl.umap(adata, color=color, show=False, title=title)
        plt.savefig(path_out, bbox_inches="tight", dpi=300)
        plt.close()
        images_list.append({"title": list_title, "path": path_out})
        return True
    except Exception as e:
        logger.error("%s: sc.pl.umap 失败 color=%s err=%s", logger_label, color, e, exc_info=True)
        try:
            plt.close()
        except Exception:
            pass
        try:
            _write_sted_plot_placeholder(path_out, _STED_PLOT_PLACEHOLDER_MSG)
            images_list.append({"title": list_title + " (placeholder)", "path": path_out})
        except Exception as pe:
            logger.error("%s: 占位图写入失败: %s", logger_label, pe, exc_info=True)
        return False


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
        valid_time_keys = _sted_ec_time_keys_ordered()
        valid_cell_type_keys = STED_EC_VALID_CELL_TYPE_KEYS
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
                    "数据校验失败：找不到时间序列列。请确保 obs 含下列任一时间相关列，或在参数中显式指定 time_key："
                    f"{list(STED_EC_VALID_TIME_KEYS)}；也可设置环境变量 STED_EC_EXTRA_TIME_KEYS 追加列名（逗号分隔）。"
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

        # 2. 严格生物学白名单：绝不可将 batch/condition 当作时间序列（sample 亦不在白名单）
        valid_time_keys = _sted_ec_time_keys_ordered()
        valid_cell_type_keys = STED_EC_VALID_CELL_TYPE_KEYS

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
                    "数据校验失败：找不到时间序列列。请确保 obs 含白名单列或在参数中显式指定 time_key："
                    f"{list(STED_EC_VALID_TIME_KEYS)}；可设置 STED_EC_EXTRA_TIME_KEYS 追加。绝不可使用 batch 或 condition 作为时间！"
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

        # 在 moscot / 任何 OT 计算前强制清洗 X 与时间列，避免 NaN/Inf 导致 Sinkhorn 发散
        _robust_sanitize_and_preprocess(adata, time_key, logger_label="Moscot")

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
        result_path = os.path.join(out_base, "adata_with_tmaps.h5ad")

        col = adata.obs[time_key]
        use_numeric = np.issubdtype(times.dtype, np.floating) or times.dtype.kind == "f"
        pairs_done = 0
        moscot_fail_msg = (
            "Moscot trajectory inference failed due to numerical instability or convergence issues. "
            "Downstream trajectory-driven analysis may be affected."
        )

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

            merged = None
            sub = None
            tp = None
            cps_adata = None
            transport = None
            try:
                # 联合 PCA 空间（复刻 Treesing runPCA 思路：拼接后 PCA）
                merged = anndata.concat([adata1, adata2], join="inner")
                if "highly_variable" not in merged.var.columns:
                    sc.pp.highly_variable_genes(
                        merged, min_mean=0.0125, max_mean=3, min_disp=0.7, n_top_genes=min(2000, merged.n_vars)
                    )
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
            except Exception as pair_e:
                logger.error(
                    "sted_ec_moscot_trajectory: Moscot 时间对 %s -> %s 失败: %s",
                    t1,
                    t2,
                    pair_e,
                    exc_info=True,
                )
            finally:
                del merged, sub, tp, cps_adata, transport
                gc.collect()
                _log_memory_mb(f"moscot 时间对 {t1}->{t2} 完成后")

        adata.uns["sted_ec_tmaps_dir"] = tmaps_dir
        adata.uns["sted_ec_time_key"] = time_key
        if cell_type_key and cell_type_key in adata.obs.columns:
            adata.uns["sted_ec_cell_type_key"] = cell_type_key

        try:
            adata.write(result_path)
            logger.info("sted_ec_moscot_trajectory: 已写出清洗后 AnnData %s", result_path)
        except Exception as w_err:
            logger.error("sted_ec_moscot_trajectory: 写出 h5ad 失败: %s", w_err, exc_info=True)
            return {
                "status": "error",
                "error": str(w_err),
                "message": moscot_fail_msg if pairs_done == 0 else f"{moscot_fail_msg} (additionally failed to save h5ad)",
            }

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

        if pairs_done == 0:
            summary = (
                f"{moscot_fail_msg} 已写出清洗后的全量 h5ad 供下游可视化：{result_path}。"
            )
            return {
                "status": "error",
                "error": moscot_fail_msg,
                "message": moscot_fail_msg,
                "h5ad_path": result_path,
                "tmaps_dir": tmaps_dir,
                "time_key": time_key,
                "moscot_pairs_solved": 0,
                "summary": summary,
                "report_data": {"images": [], "download_links": download_links, "summary": summary},
            }

        summary = f"已完成 {pairs_done} 个时间对 OT 求解，transport 矩阵已写入 tmaps；清洗后结果 h5ad 已保存。"
        return {
            "status": "success",
            "message": "moscot 轨迹推断完成，transport 矩阵已写入 tmaps",
            "h5ad_path": result_path,
            "tmaps_dir": tmaps_dir,
            "time_key": time_key,
            "moscot_pairs_solved": pairs_done,
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
            for k in STED_EC_VALID_CELL_TYPE_KEYS:
                if k in adata.obs.columns:
                    cell_type_key = k
                    logger.info("sted_ec_plot_trajectory 智能嗅探到细胞类型列: %s", k)
                    break
        if time_key not in adata.obs.columns:
            return {"status": "error", "error": f"obs 中缺少时间列: {time_key}", "report_data": {"images": [], "download_links": [], "summary": ""}}

        _robust_sanitize_and_preprocess(adata, time_key, logger_label="sted_ec_plot_trajectory")

        out_dir = output_dir or os.path.dirname(trajectory_data_path)
        if output_plot_path:
            out_dir = os.path.dirname(output_plot_path)
        img_dir = os.path.join(out_dir, "sted_ec_report_images")
        os.makedirs(img_dir, exist_ok=True)
        h5ad_out = os.path.join(out_dir, "sted_ec_after_umap.h5ad")

        need_reembed = ("X_umap" not in adata.obsm) or (not _sted_obsm_all_finite(adata, "X_umap"))
        if need_reembed:
            logger.info("sted_ec_plot_trajectory: X_umap 缺失或含非有限值，重算 PCA + neighbors + UMAP")
            try:
                _sted_ec_safe_recompute_embedding(adata, "sted_ec_plot_trajectory")
            except Exception as emb_e:
                logger.error("sted_ec_plot_trajectory: 重算嵌入失败（仍将尝试绘图占位并写出 h5ad）: %s", emb_e, exc_info=True)
                _nan_to_num_obsm(adata, "X_pca")
                _nan_to_num_obsm(adata, "X_umap")
        else:
            _nan_to_num_obsm(adata, "X_pca")
            _nan_to_num_obsm(adata, "X_umap")

        images_list: List[Dict[str, str]] = []

        path1 = os.path.join(img_dir, "umap_time.png")
        _sanitize_obs_color_column(adata, time_key)
        _safe_sted_pl_umap(
            adata,
            time_key,
            "Spatiotemporal UMAP (by Time)",
            path1,
            images_list,
            "Spatiotemporal UMAP",
            "sted_ec_plot_trajectory",
        )
        logger.info("sted_ec_plot_trajectory: 时空 UMAP 输出路径 %s", path1)

        # 图表 2 & 3：细胞类型相关 (动态执行，有则画，无则跳过) —— 独立使用 path2/path3，避免变量污染
        if cell_type_key and cell_type_key in adata.obs.columns:
            try:
                path2 = os.path.join(img_dir, "umap_celltype.png")
                _sanitize_obs_color_column(adata, cell_type_key)
                _safe_sted_pl_umap(
                    adata,
                    cell_type_key,
                    "Cell Type UMAP",
                    path2,
                    images_list,
                    "Cell Type UMAP",
                    "sted_ec_plot_trajectory",
                )
                logger.info("sted_ec_plot_trajectory: Cell Type UMAP 输出路径 %s", path2)

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

        # DAG 铁律：无论 UMAP 绘图是否成功，优先写出清洗后的 h5ad，供基因挖掘等下游使用
        try:
            adata.write(h5ad_out)
            logger.info("sted_ec_plot_trajectory: 已写出下游用 AnnData %s", h5ad_out)
        except Exception as w_err:
            logger.error("sted_ec_plot_trajectory: 写出 h5ad 失败: %s", w_err, exc_info=True)
            del adata
            gc.collect()
            return {
                "status": "error",
                "error": f"无法写出下游 h5ad: {w_err}",
                "report_data": {"images": images_list, "download_links": [], "summary": ""},
            }

        zip_base = os.path.join(out_dir, "sted_ec_results")
        zip_path = zip_base + ".zip"
        download_links: List[Dict[str, str]] = []
        try:
            if os.path.isfile(zip_path):
                os.remove(zip_path)
            shutil.make_archive(zip_base, "zip", out_dir, "sted_ec_report_images")
            logger.info("Packed report archive: %s", zip_path)
            download_links = [{"title": "📦 Download Full Results (ZIP)", "path": zip_path}]
        except Exception as zip_e:
            logger.warning("sted_ec_plot_trajectory: ZIP 打包失败（h5ad 已写出，下游可继续）: %s", zip_e, exc_info=True)
            download_links = []

        plot_degraded = any("(placeholder)" in str(it.get("title", "")) for it in images_list)
        summary = (
            "已生成时空 UMAP、细胞类型 UMAP、细胞类型演化图；所有图表已打包至 ZIP 供下载。"
            if images_list
            else "轨迹图表已生成并打包至 ZIP。"
        )
        if plot_degraded:
            summary += " 部分 UMAP 因数据分布异常已降级为占位图，下游 h5ad 仍可用于基因挖掘。"
        if not download_links:
            summary += " 报告 ZIP 打包失败，请直接使用 sted_ec_after_umap.h5ad 与单张 PNG。"
        success_msg = "轨迹多维图表已保存；已输出含 UMAP 的 h5ad 供下游分析"
        if download_links:
            success_msg += "并打包 ZIP"
        if plot_degraded:
            success_msg += "（部分图表为占位图，详见 report_data.summary）"
        del adata
        gc.collect()
        return {
            "status": "success",
            "message": success_msg,
            "output_plot_path": path1,
            "h5ad_path": h5ad_out,
            "output_h5ad": h5ad_out,
            "time_key_used": time_key,
            "cell_type_key_used": cell_type_key,
            "report_data": {
                "images": images_list,
                "download_links": download_links,
                "summary": summary,
            },
        }
    except Exception as e:
        logger.error("sted_ec_plot_trajectory 执行失败: %s", e, exc_info=True)
        return {"status": "error", "error": str(e), "report_data": {"images": [], "download_links": [], "summary": ""}}


@registry.register(
    name="sted_ec_driver_gene_extraction",
    description="轨迹驱动基因挖掘：对第 4 步输出的含 UMAP/时序 obs 的 h5ad 使用 scanpy.tl.rank_genes_groups (wilcoxon) 提取候选基因表 driver_genes.csv，供 GSEA 等下游使用。",
    category="STED_EC",
    output_type="mixed",
)
def sted_ec_driver_gene_extraction(
    h5ad_path: str,
    output_dir: Optional[str] = None,
    top_n: int = 100,
    groupby_mode: str = "auto",
    time_key: Optional[str] = None,
    cell_type_key: Optional[str] = None,
    pval_adj_max: float = 0.05,
    min_abs_logfoldchange: float = 0.25,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    承上启下：仅使用 Scanpy 原生 API，不做自定义假设检验实现。

    Args:
        h5ad_path: 通常为第 4 步 sted_ec_plot_trajectory 写出的 sted_ec_after_umap.h5ad。
        output_dir: CSV 输出目录，默认与 h5ad 同目录。
        top_n: 导出表格保留的基因行数上限（去重后）。
        groupby_mode: 'auto' | 'time' | 'cell_type' — rank_genes_groups 的分组列优先策略。
        time_key / cell_type_key: 可选；缺省时从 adata.uns 与前序列名推断。
    """
    try:
        import numpy as np
        import pandas as pd
        import scanpy as sc
    except ImportError as e:
        trace_str = traceback.format_exc()
        logger.error("Tool execution failed: sted_ec_driver_gene_extraction (ImportError)\n%s", trace_str)
        return _sted_safe_return({"status": "error", "error": str(e), "traceback": trace_str})

    if isinstance(h5ad_path, str) and h5ad_path.strip().startswith("<") and h5ad_path.strip().endswith(">"):
        return _sted_safe_return({
            "status": "error",
            "error": "h5ad 路径未解析。请先成功完成「时空动力学图可视化」步骤。",
        })

    path_resolved = str(Path(h5ad_path).resolve())
    if not os.path.isfile(path_resolved):
        return _sted_safe_return({"status": "error", "error": f"文件不存在: {path_resolved}"})

    out_base = output_dir or os.path.dirname(path_resolved)
    os.makedirs(out_base, exist_ok=True)
    csv_path = os.path.join(out_base, "driver_genes.csv")

    adata = None
    try:
        adata = sc.read_h5ad(path_resolved)
        tk = time_key or adata.uns.get("sted_ec_time_key", "day")
        ck = cell_type_key if (cell_type_key and str(cell_type_key).strip()) else adata.uns.get("sted_ec_cell_type_key")

        gb_source: Optional[str] = None
        mode = (groupby_mode or "auto").lower().strip()
        if mode == "time":
            gb_source = tk if tk and tk in adata.obs.columns else None
        elif mode == "cell_type":
            gb_source = ck if ck and ck in adata.obs.columns else None
        else:
            if tk and tk in adata.obs.columns:
                nu_t = adata.obs[tk].astype(str).nunique()
                if 2 <= nu_t <= 80:
                    gb_source = tk
            if gb_source is None and ck and ck in adata.obs.columns:
                nu_c = adata.obs[ck].astype(str).nunique()
                if nu_c >= 2:
                    gb_source = ck
            if gb_source is None and tk and tk in adata.obs.columns:
                gb_source = tk

        if not gb_source or gb_source not in adata.obs.columns:
            return _sted_safe_return({
                "status": "error",
                "error": "无法在 obs 中解析用于 rank_genes_groups 的分组列（time_key/cell_type_key）。",
            })

        target_col = gb_source
        # 1. 提取目标列，并强制、立刻转换为普通的 object/string 类型！绝不保留 Categorical！
        raw_series = adata.obs[target_col].astype(str)

        # 2. 将所有可能的空值形态（"nan", "NaN", "None", "", 纯空格）统一替换为 "Unknown"
        import numpy as np

        raw_series = raw_series.replace(["nan", "NaN", "None", "", " "], "Unknown")
        raw_series = raw_series.replace(r"^\s*$", "Unknown", regex=True)
        raw_series = raw_series.fillna("Unknown")

        # 3. 赋值给新的临时列
        adata.obs["_sted_driver_groupby"] = raw_series

        # 4. 检查唯一值数量
        unique_vals = adata.obs["_sted_driver_groupby"].nunique()
        if unique_vals <= 1:
            del adata.obs["_sted_driver_groupby"]
            return _sted_safe_return({
                "status": "error",
                "error": f"分组变量 {target_col} 只有一个有效类别，无法进行差异基因分析。",
            })

        # 5. 最后一步：安全地转换为 category 类型供 Scanpy 使用
        adata.obs["_sted_driver_groupby"] = adata.obs["_sted_driver_groupby"].astype("category")

        use_raw = adata.raw is not None
        sc.tl.rank_genes_groups(
            adata,
            groupby="_sted_driver_groupby",
            method="wilcoxon",
            use_raw=use_raw,
            key_added="rank_genes_driver",
        )

        df = None
        try:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_driver", group=None)
        except Exception:
            df = None
        if df is None or df.empty:
            del adata.obs["_sted_driver_groupby"]
            return _sted_safe_return({"status": "error", "error": "rank_genes_groups 未产生可用结果（get.rank_genes_groups_df 为空）。"})

        df_f = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["pvals_adj"])
        df_f = df_f[df_f["pvals_adj"] < float(pval_adj_max)]
        if "logfoldchanges" in df_f.columns:
            df_f = df_f[np.isfinite(df_f["logfoldchanges"])]
            df_f = df_f[df_f["logfoldchanges"].abs() >= float(min_abs_logfoldchange)]

        note_suffix = ""
        if df_f.empty:
            df_f = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["scores"])
            df_f = df_f.sort_values("scores", ascending=False)
            note_suffix = "（统计阈值下过检基因为空，已按 score 取 Top 去重）"
        else:
            df_f = df_f.sort_values("scores", ascending=False)

        df_f = df_f.drop_duplicates(subset=["names"]).head(int(top_n))
        df_f.to_csv(csv_path, index=False)

        # Top 5 Markdown（展示用）；清洗 NaN/Inf，避免 markdown / JSON 序列化异常
        # 注意：rank_genes_groups_df 常含 Categorical 列（如 group），对 Categorical 直接 fillna("")
        # 会触发「Cannot setitem on a Categorical with a new category ()」，须先转为 object/str
        preview = df_f.head(5).copy()
        preview = preview.replace([np.inf, -np.inf], np.nan)
        for _pcol in preview.columns:
            if pd.api.types.is_categorical_dtype(preview[_pcol]):
                preview[_pcol] = preview[_pcol].astype(object)
        preview = preview.fillna("")
        try:
            md_table = preview.to_markdown(index=False)
        except Exception:
            md_table = "```\n" + preview.to_string(index=False) + "\n```"
        driver_genes_markdown = (
            f"**轨迹驱动基因挖掘**（groupby=`{gb_source}`，method=wilcoxon，导出 Top {len(df_f)} → `driver_genes.csv`）{note_suffix}\n\n"
            f"{md_table}\n"
        )

        # UI：折叠面板 summary + Markdown 表格；report_data.download_links 供前端下载按钮（title/path 与 index.html 一致）
        summary_str = (
            f"已导出驱动基因表: {csv_path}\n\n**驱动基因 Top 5 预览：**\n{driver_genes_markdown}"
        )
        report_data = {
            "images": [],
            "download_links": [{"title": "下载驱动基因表 (CSV)", "path": csv_path}],
            "summary": summary_str,
        }

        return _sted_safe_return({
            "status": "success",
            # message 须含完整 Markdown：executor 用 message 优先作为 step_detail.summary
            "message": summary_str,
            "summary": summary_str,
            "driver_genes_csv": csv_path,
            "csv_path": csv_path,
            "output_file": csv_path,
            "groupby_used": str(gb_source),
            "n_genes_written": int(len(df_f)),
            "driver_genes_markdown": driver_genes_markdown,
            "report_data": report_data,
        })
    except Exception as e:
        trace_str = traceback.format_exc()
        logger.error("Tool execution failed: sted_ec_driver_gene_extraction\n%s", trace_str)
        return _sted_safe_return({"status": "error", "error": str(e), "traceback": trace_str})
    finally:
        if adata is not None:
            try:
                if "_sted_driver_groupby" in adata.obs.columns:
                    del adata.obs["_sted_driver_groupby"]
            except Exception:
                pass
            del adata
        gc.collect()


@registry.register(
    name="sted_ec_pathway_enrichment",
    description="STED_EC 通路富集：读取轨迹驱动基因 CSV，调用 gseapy.enrichr（KEGG/GO 等 Enrichr 库），输出富集表 CSV 与条形图/气泡图 PNG，供专家报告引用。",
    category="STED_EC",
    output_type="mixed",
)
def sted_ec_pathway_enrichment(
    driver_genes_csv: str,
    output_dir: Optional[str] = None,
    organism: str = "Human",
    gene_sets: Optional[str] = None,
    top_pathways_plot: int = 15,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    与 spatial_pathway_enrichment 能力对齐：复用相同 gseapy.enrichr + 本地绑图模式（no_plot=True）。

    Args:
        driver_genes_csv: 第 5 步 driver_genes.csv 的绝对路径（含 scanpy 的 names 列）。
        output_dir: 结果输出目录；默认同 CSV 所在目录。
        organism: Enrichr organism（Human / Mouse 等）。
        gene_sets: 逗号分隔的 Enrichr 库名；默认 ``KEGG_2021_Human,GO_Biological_Process_2021``。
        top_pathways_plot: 绑图中的通路条数上限。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except ImportError as e:
        return _sted_safe_return({"status": "error", "error": str(e)})

    try:
        # Executor 占位符若未命中路径字段，可能把整段 step_result dict 传入，需兜底解析
        if isinstance(driver_genes_csv, dict):
            _d = driver_genes_csv
            driver_genes_csv = (
                _d.get("driver_genes_csv")
                or _d.get("csv_path")
                or _d.get("output_file")
                or ""
            )
            if not isinstance(driver_genes_csv, str) or not str(driver_genes_csv).strip():
                return _sted_safe_return({
                    "status": "error",
                    "error": "driver_genes_csv 应为 CSV 文件路径，当前收到 dict 且无法提取路径。",
                })

        if isinstance(driver_genes_csv, str) and driver_genes_csv.strip().startswith("<") and driver_genes_csv.strip().endswith(">"):
            return _sted_safe_return({
                "status": "error",
                "error": "driver_genes_csv 未解析，请先成功完成「轨迹驱动基因挖掘」步骤。",
            })

        csv_in = str(Path(str(driver_genes_csv)).resolve())
        if not os.path.isfile(csv_in):
            return _sted_safe_return({"status": "error", "error": f"驱动基因文件不存在: {csv_in}"})
    except Exception as e:
        logger.error("sted_ec_pathway_enrichment 路径校验异常: %s", e, exc_info=True)
        return _sted_safe_return({"status": "error", "error": str(e)})

    out_base = output_dir or os.path.dirname(csv_in)
    os.makedirs(out_base, exist_ok=True)
    base_tag = "sted_ec_pathway_enrichment"
    out_csv = os.path.join(out_base, f"{base_tag}_enrichment_results.csv")
    bar_png = os.path.join(out_base, f"{base_tag}_bar.png")
    bubble_png = os.path.join(out_base, f"{base_tag}_bubble.png")

    gs_default = "KEGG_2021_Human,GO_Biological_Process_2021"
    gs_raw = (gene_sets or gs_default).strip()
    gene_sets_list = [s.strip() for s in gs_raw.split(",") if s.strip()]
    if not gene_sets_list:
        gene_sets_list = ["KEGG_2021_Human", "GO_Biological_Process_2021"]

    try:
        df_in = pd.read_csv(csv_in)
        col_gene = None
        for c in ("names", "gene", "Gene", "symbol", "gene_symbol"):
            if c in df_in.columns:
                col_gene = c
                break
        if col_gene is None and len(df_in.columns) > 0:
            col_gene = df_in.columns[0]
        if col_gene is None:
            return _sted_safe_return({"status": "error", "error": "driver_genes.csv 中未找到基因名列（期望 names）。"})

        genes_ser = (
            df_in[col_gene]
            .dropna()
            .astype(str)
            .str.strip()
        )
        genes = genes_ser[genes_ser.str.len() > 0].unique().tolist()
        if len(genes) < 3:
            return _sted_safe_return({"status": "error", "error": f"有效基因数过少（{len(genes)}），无法进行富集。"})

        enrichment_table: Optional[pd.DataFrame] = None
        fallback_message = None

        try:
            import gseapy as gp

            # gseapy.enrichr 仅接受小写 organism（如 human），拒绝 "Human"；空串回退 human
            _org = str(organism or "").strip().lower()
            safe_organism = _org if _org else "human"

            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_sets_list,
                organism=safe_organism,
                outdir=None,
                no_plot=True,
            )
            enrichment_table = getattr(enr, "results", enr) if enr is not None and hasattr(enr, "results") else enr
            if not isinstance(enrichment_table, pd.DataFrame) or enrichment_table.empty:
                enrichment_table = None
                fallback_message = "gseapy.enrichr 返回空结果。"
        except Exception as ex:
            logger.warning("sted_ec_pathway_enrichment: gseapy 失败: %s", ex, exc_info=True)
            fallback_message = f"gseapy 未执行成功: {ex}"
            enrichment_table = None

        if enrichment_table is None or enrichment_table.empty:
            pd.DataFrame({"gene": genes[:500]}).to_csv(os.path.join(out_base, f"{base_tag}_input_genes_fallback.csv"), index=False)
            return _sted_safe_return({
                "status": "error",
                "error": str(fallback_message or "富集无结果"),
                "message": str(fallback_message or "富集无结果"),
                "n_genes_input": int(len(genes)),
                "report_data": {
                    "images": [],
                    "download_links": [],
                    "summary": str(fallback_message or "富集失败"),
                },
            })

        enrichment_table = enrichment_table.copy()
        enrichment_table = enrichment_table.replace([np.inf, -np.inf], np.nan)
        enrichment_table.to_csv(out_csv, index=False)

        pcol = "Adjusted P-value" if "Adjusted P-value" in enrichment_table.columns else "P-value"
        term_col = "Term" if "Term" in enrichment_table.columns else enrichment_table.columns[0]
        plot_df = enrichment_table.sort_values(pcol, ascending=True).head(int(top_pathways_plot)).copy()
        pvals = plot_df[pcol].astype(float)
        pvals = np.clip(pvals, 1e-300, 1.0)
        logp = -np.log10(pvals)
        y_pos = np.arange(len(plot_df))

        fig1, ax1 = plt.subplots(figsize=(9, max(4.5, len(plot_df) * 0.38)))
        ax1.barh(y_pos, logp, color="steelblue", alpha=0.85)
        ax1.set_yticks(y_pos)
        labels = plot_df[term_col].astype(str).str.slice(0, 60)
        ax1.set_yticklabels(labels, fontsize=9)
        ax1.set_xlabel("-log10(Adjusted P-value)" if pcol == "Adjusted P-value" else "-log10(P-value)")
        ax1.set_title("STED-EC pathway enrichment (Enrichr)")
        ax1.invert_yaxis()
        fig1.tight_layout()
        fig1.savefig(bar_png, dpi=150)
        plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(9, max(4.5, len(plot_df) * 0.38)))
        area = 80.0 + 420.0 * (logp / (float(np.max(logp)) or 1.0))
        ax2.scatter(logp, y_pos, s=area, c=logp, cmap="viridis", alpha=0.75, edgecolors="k", linewidths=0.3)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(labels, fontsize=9)
        ax2.set_xlabel("-log10 adj.P / P")
        ax2.set_title("Pathway enrichment (bubble view)")
        ax2.invert_yaxis()
        fig2.tight_layout()
        fig2.savefig(bubble_png, dpi=150)
        plt.close(fig2)

        preview = enrichment_table.head(5).copy()
        preview = preview.replace([np.inf, -np.inf], np.nan)
        for _pcol in preview.columns:
            if pd.api.types.is_categorical_dtype(preview[_pcol]):
                preview[_pcol] = preview[_pcol].astype(object)
        preview = preview.fillna("")
        try:
            md_table = preview.to_markdown(index=False)
        except Exception:
            md_table = "```\n" + preview.to_string(index=False) + "\n```"
        pathway_enrichment_markdown = (
            f"**通路富集（Enrichr: {', '.join(gene_sets_list)}）**\n前 5 条：\n\n{md_table}\n"
        )

        summary_text = (
            f"已对 {len(genes)} 个输入基因做 Enrichr 富集；"
            f"显著条目数 {len(enrichment_table)}；结果表见 {os.path.basename(out_csv)}。"
        )

        gc.collect()

        return _sted_safe_return({
            "status": "success",
            "message": summary_text,
            "enrichment_csv": str(out_csv),
            "csv_path": str(out_csv),
            "output_file": str(out_csv),
            "pathway_enrichment_markdown": pathway_enrichment_markdown,
            "n_genes_input": int(len(genes)),
            "n_terms": int(len(enrichment_table)),
            "gene_sets_used": [str(x) for x in gene_sets_list],
            "summary": {
                "text": summary_text,
                "pathway_enrichment_markdown": pathway_enrichment_markdown,
                "n_terms": int(len(enrichment_table)),
            },
            "report_data": {
                "images": [
                    {"title": "通路富集 条形图", "path": str(bar_png)},
                    {"title": "通路富集 气泡图", "path": str(bubble_png)},
                ],
                "download_links": [
                    {"title": "富集结果 CSV", "path": str(out_csv)},
                ],
                "summary": pathway_enrichment_markdown,
            },
        })
    except Exception as e:
        logger.error("sted_ec_pathway_enrichment 失败: %s", e, exc_info=True)
        gc.collect()
        return _sted_safe_return({"status": "error", "error": str(e)})


@registry.register(
    name="sted_ec_expert_report",
    description="[遗留/手动] STED-EC 专家解读：仅建议在非工作流场景单独调用。标准 SPATIOTEMPORAL/STED 流程的报告由 Orchestrator 在 execute 之后经 BaseAgent._generate_analysis_summary 生成。",
    category="STED_EC",
    output_type="json",
)
def sted_ec_expert_report(
    read_only_context: str = "{}",
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    专家解读报告（原子 DAG 节点）。可选启用 MCP：与闲聊一致注入 mcp_web_search / mcp_ncbi_search 工具并多轮调用。

    Args:
        read_only_context: Executor 从前序 sted_ec_* 步骤组装的 JSON 字符串（只读）。
        enabled_mcps: 由 Executor 从 step_context 注入，如 ["web_search","authority_db"]。

    Returns:
        含 expert_report_markdown / report、report_data.references，供 SSE steps_details 与前端「参考文献」面板。
    """
    from gibh_agent.agents.reporting.sted_ec_expert_report_prompts import (
        STED_EC_CHILD_AGENT_SYSTEM_PROMPT,
        STED_EC_EXPERT_SYSTEM_WITH_MCP,
        build_sted_ec_child_user_message,
    )
    from gibh_agent.core.llm_client import LLMClientFactory
    from gibh_agent.core.openai_tools import (
        apply_hpc_mcp_tool_policy,
        hpc_chat_system_suffix,
        tool_names_to_openai_tools,
    )

    try:
        ctx_raw = (read_only_context or "").strip() or "{}"
        try:
            parsed = json.loads(ctx_raw)
        except json.JSONDecodeError:
            return {
                "status": "error",
                "error": "read_only_context 不是合法 JSON",
                "message": "专家报告步骤失败：上下文解析错误",
            }
        steps = parsed.get("steps") if isinstance(parsed, dict) else None
        if not steps:
            return {
                "status": "error",
                "error": "缺少前序 STED-EC 步骤上下文，无法生成解读",
                "message": "专家报告步骤失败：前序步骤结果为空",
            }

        enabled_mcps = _normalize_enabled_mcps_for_tools(kwargs.get("enabled_mcps"))
        tool_names: List[str] = []
        if "web_search" in enabled_mcps:
            tool_names.append("mcp_web_search")
        if "authority_db" in enabled_mcps:
            tool_names.append("mcp_ncbi_search")
        tool_names = apply_hpc_mcp_tool_policy(tool_names, enabled_mcps)

        system_prompt = (
            STED_EC_EXPERT_SYSTEM_WITH_MCP if tool_names else STED_EC_CHILD_AGENT_SYSTEM_PROMPT
        )
        system_prompt += hpc_chat_system_suffix(enabled_mcps)
        user_content = build_sted_ec_child_user_message(ctx_raw)

        from gibh_agent.core.llm_cloud_providers import validate_and_resolve_model_name

        _mn = validate_and_resolve_model_name((kwargs.get("model_name") or "").strip() or None)
        client = LLMClientFactory.create_for_model(_mn)
        messages: List[Dict[str, Any]] = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content},
        ]

        tools = tool_names_to_openai_tools(tool_names) if tool_names else []
        max_rounds = 12
        content = ""

        if not tools:
            completion = client.chat(
                messages=messages,
                temperature=0.35,
                max_tokens=8192,
            )
            ch = getattr(completion, "choices", None) or []
            if not ch:
                return {
                    "status": "error",
                    "error": "LLM 返回空 choices",
                    "message": "专家报告生成失败：模型无输出",
                }
            content = (getattr(ch[0].message, "content", None) or "").strip()
        else:
            for _round in range(max_rounds):
                completion = client.chat(
                    messages=messages,
                    tools=tools,
                    tool_choice="auto",
                    temperature=0.35,
                    max_tokens=8192,
                )
                ch = getattr(completion, "choices", None) or []
                if not ch:
                    return {
                        "status": "error",
                        "error": "LLM 返回空 choices",
                        "message": "专家报告生成失败：模型无输出",
                    }
                msg = ch[0].message
                tcalls = getattr(msg, "tool_calls", None)
                if tcalls:
                    messages.append(_assistant_msg_to_openai_dict(msg))
                    for tc in tcalls:
                        fn = getattr(tc, "function", None)
                        name = getattr(fn, "name", None) if fn is not None else None
                        args_raw = (getattr(fn, "arguments", None) if fn is not None else None) or "{}"
                        if not name:
                            continue
                        try:
                            args = json.loads(args_raw)
                        except json.JSONDecodeError:
                            args = {}
                        tool_fn = registry.get_tool(name)
                        if name.startswith("hpc_mcp_") and "compute_scheduler" not in enabled_mcps:
                            payload = {
                                "status": "error",
                                "message": "HPC Error: 未授权使用超算工具，请提示用户开启「计算资源智能调度」。",
                            }
                        elif not tool_fn:
                            payload = {"status": "error", "error": f"未知工具: {name}"}
                        else:
                            try:
                                if inspect.iscoroutinefunction(tool_fn):
                                    import asyncio

                                    payload = asyncio.run(tool_fn(**args))
                                else:
                                    payload = tool_fn(**args)
                            except Exception as tool_exc:  # noqa: BLE001
                                logger.exception("sted_ec_expert_report 工具执行异常: %s", name)
                                if name.startswith("hpc_mcp_"):
                                    payload = {
                                        "status": "error",
                                        "message": f"HPC Error: execution failed, please inform the user. ({tool_exc})",
                                    }
                                else:
                                    payload = {"status": "error", "error": str(tool_exc), "message": str(tool_exc)}
                        messages.append(
                            {
                                "role": "tool",
                                "tool_call_id": getattr(tc, "id", ""),
                                "content": json.dumps(payload, ensure_ascii=False),
                            }
                        )
                    continue
                content = (getattr(msg, "content", None) or "").strip()
                break
            else:
                return {
                    "status": "error",
                    "error": "工具轮次超限",
                    "message": "专家报告生成失败：工具调用轮次过多",
                }

        if not content:
            return {
                "status": "error",
                "error": "LLM 内容为空",
                "message": "专家报告生成失败：正文为空",
            }
        if content.startswith("```"):
            lines = content.split("\n")
            if len(lines) >= 2 and lines[0].strip().startswith("```"):
                lines = lines[1:]
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            content = "\n".join(lines).strip()

        report_md, references = _extract_references_from_content(content)

        payload = {
            "status": "success",
            "message": "STED-EC 专家解读报告已生成",
            "summary": {
                "expert_report_generated": True,
                "source": "sted_ec_expert_report",
                "mcp_tools_used": tool_names,
            },
            "report": report_md,
            "expert_report_markdown": report_md,
            "report_node": "sted_ec_expert_report",
            "report_data": {
                "references": references,
                "report": report_md,
            },
        }
        return _sted_safe_return(payload)
    except Exception as e:
        logger.error("sted_ec_expert_report 失败: %s", e, exc_info=True)
        return {
            "status": "error",
            "error": str(e),
            "message": f"专家报告步骤异常: {e}",
        }
