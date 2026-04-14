"""
PySkills Worker：传统生信算法综合微服务（ViennaRNA / PyMOL / RDKit / gseapy）。
上传目录与主 API 通过 docker-compose 挂载同一宿主机路径：/app/uploads。
富媒体产物写入 /app/uploads/results/pyskills/<run_id>/，对外仅返回 /uploads/results/... URL（禁止泄露容器内绝对路径）。
"""

from __future__ import annotations

import csv
import io
import logging
import os
import re
import subprocess
import sys
import threading
import traceback
import uuid
from pathlib import Path

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

import RNA
import pymol
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski as RDKitLipinski
import gseapy

from sa_score_compat import load_calculate_score

_calculate_score_fn = None


def _get_calculate_score():
    global _calculate_score_fn
    if _calculate_score_fn is None:
        _calculate_score_fn = load_calculate_score()
    return _calculate_score_fn

app = FastAPI(
    title="PySkills Worker API",
    description="传统生信算法综合微服务（ViennaRNA, PyMOL, RDKit, gseapy）",
)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# 与 nginx `location /uploads/` -> /app/uploads/ 对齐；响应中只出现 /uploads/... 相对路径
PYSKILLS_PUBLIC_BASE = "/uploads/results/pyskills"
_pymol_lock = threading.Lock()
_pymol_launcher_started = False

# 与 docker-compose 挂载一致：宿主机 ./third_party/drug-similarity -> 容器 /app/third_party/drug-similarity
DRUG_SIM_BUNDLE = Path(os.environ.get("DRUG_SIMILARITY_BUNDLE", "/app/third_party/drug-similarity")).resolve()


def _pyskills_out_dir() -> Path:
    d = Path("/app/uploads/results/pyskills")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _public_url(run_id: str, filename: str) -> str:
    return f"{PYSKILLS_PUBLIC_BASE}/{run_id}/{filename}"


class ToolRequest(BaseModel):
    file_path: str


class RnafoldRequest(BaseModel):
    """与智能体填参对齐：FASTA 路径 + 可选温度与配对概率（RNAfold -p 风格）。"""

    file_path: str
    temperature_celsius: float = Field(default=37.0, ge=0.0, le=100.0, description="折叠温度（℃），启用温度控制时生效")
    use_temperature_control: bool = Field(default=True, description="False 时使用 ViennaRNA 默认模型细节（不覆写 md.temperature）")
    include_pairing_probabilities: bool = Field(
        default=False,
        description="True 时计算分区函数并返回 pair_prob_structure（点括号+概率信息，用于 dot plot / 矩阵可视化）",
    )


def _read_text_safe(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        raise HTTPException(status_code=400, detail=f"无法读取文件: {e}") from e


def _first_fasta_sequence(text: str) -> str:
    """从 FASTA 文本中提取第一条序列（忽略 header 行）。"""
    lines: list[str] = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith(">"):
            if lines:
                break
            continue
        lines.append(s)
    seq = "".join(lines).strip().upper()
    seq = re.sub(r"\s+", "", seq)
    if not seq:
        raise ValueError("FASTA 中未找到有效序列")
    allowed = set("AUCGTUN")
    bad = [c for c in seq if c not in allowed]
    if bad:
        raise ValueError(f"序列含非法碱基字符: {set(bad)!r}（仅允许 RNA: AUCGTU 或 N）")
    return seq


def _first_smiles_from_text(text: str) -> str:
    """从文本中取第一个非注释、非空行作为 SMILES（或行内第一个 token）。"""
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith(">"):
            continue
        token = s.split()[0].strip()
        if token:
            return token
    raise ValueError("未找到 SMILES（请每行一个 SMILES 或首列为 SMILES）")


def _rna_fold_advanced(
    sequence: str,
    *,
    use_temperature_control: bool,
    temperature_celsius: float,
    include_pairing_probabilities: bool,
) -> dict:
    """
    MFE + 可选分区函数（碱基配对概率，RNAfold -p 类输出）。
    返回 dict：structure, mfe, 以及可选 pair_prob_structure, partition_free_energy。
    """
    try:
        if not use_temperature_control and not include_pairing_probabilities:
            structure, mfe = RNA.fold(sequence)
            mfe_f = float(mfe)
            return {
                "structure": str(structure),
                "mfe": mfe_f,
                "temperature_celsius": None,
                "use_temperature_control": False,
            }

        md = RNA.md()
        if use_temperature_control:
            md.temperature = float(temperature_celsius)
        fc = RNA.fold_compound(sequence, md)
        structure, mfe = fc.mfe()
        mfe_f = float(mfe)
        out: dict = {
            "structure": str(structure),
            "mfe": mfe_f,
            "temperature_celsius": float(temperature_celsius) if use_temperature_control else None,
            "use_temperature_control": use_temperature_control,
        }
        if include_pairing_probabilities:
            fc.exp_params_rescale(mfe)
            pf_ret = fc.pf()
            if isinstance(pf_ret, tuple) and len(pf_ret) >= 2:
                pp_str, pf_energy = pf_ret[0], pf_ret[1]
            else:
                pp_str, pf_energy = pf_ret, 0.0
            out["pair_prob_structure"] = str(pp_str)
            try:
                out["partition_free_energy"] = float(pf_energy)
            except (TypeError, ValueError):
                out["partition_free_energy"] = float(getattr(pf_energy, "real", pf_energy))
        return out
    except Exception as exc:
        logger.exception("ViennaRNA fold / pf")
        raise ValueError(f"ViennaRNA 计算失败: {exc}") from exc


@app.get("/health")
async def health():
    return {"status": "ok", "service": "worker-pyskills"}


@app.post("/api/rnafold")
async def run_rnafold(req: RnafoldRequest):
    try:
        p = Path(req.file_path)
        if not p.is_file():
            raise HTTPException(status_code=400, detail=f"文件不存在: {req.file_path}")
        raw = _read_text_safe(p)
        try:
            sequence = _first_fasta_sequence(raw).replace("T", "U")
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e)) from e
        try:
            data = _rna_fold_advanced(
                sequence,
                use_temperature_control=req.use_temperature_control,
                temperature_celsius=req.temperature_celsius,
                include_pairing_probabilities=req.include_pairing_probabilities,
            )
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e)) from e
        return {"status": "success", "data": data}
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("run_rnafold")
        raise HTTPException(status_code=500, detail=traceback.format_exc()) from e


def _pymol_ensure_launcher() -> None:
    global _pymol_launcher_started
    with _pymol_lock:
        if _pymol_launcher_started:
            return
        try:
            pymol.finish_launching(["pymol", "-c"])
        except Exception as e:
            logger.warning("PyMOL finish_launching: %s", e)
        _pymol_launcher_started = True


def _pymol_color_cartoon_robust(selection: str) -> None:
    """
    为 cartoon 着色。说明：
    - 禁止 cmd.color(\"ss\", …)：\"ss\" 不是颜色名 → Unknown color。
    - 部分精简构建下 util.cbag 使用的 \"atomic\" 等着色关键字也会报 Unknown color，
      故优先 spectrum / 基础色名，再尝试 util。
    """
    from pymol import cmd

    for expr, palette in (
        ("count", "slate_cyan"),
        ("count", "blue_red"),
        ("b", "blue_white_red"),
        ("count", "rainbow"),
    ):
        try:
            cmd.spectrum(expr, palette, selection)
            return
        except Exception as e:
            logger.debug("PyMOL spectrum(%s,%s) 失败: %s", expr, palette, e)

    for cname in (
        "cyan",
        "marine",
        "yellow",
        "red",
        "blue",
        "green",
        "magenta",
        "orange",
        "violet",
        "salmon",
        "tv_blue",
        "grey70",
        "gray70",
        "gray",
    ):
        try:
            cmd.color(cname, selection)
            return
        except Exception:
            continue

    try:
        from pymol import util  # type: ignore

        for fn_name, args in (("cbag", (selection,)), ("cbac", (selection,))):
            fn = getattr(util, fn_name, None)
            if callable(fn):
                try:
                    fn(*args)
                    return
                except Exception as e:
                    logger.debug("PyMOL util.%s 失败: %s", fn_name, e)

        if callable(getattr(util, "cbss", None)):
            try:
                cmd.dss(selection)
            except Exception:
                pass
            try:
                util.cbss(selection)
                return
            except Exception as e:
                logger.debug("PyMOL util.cbss 失败: %s", e)
    except Exception as e:
        logger.debug("PyMOL util 导入或调用异常: %s", e)

    logger.warning("PyMOL 未能为 %s 设置颜色，使用默认 cartoon 颜色", selection)


def _pymol_color_sticks_simple(selection: str) -> None:
    """单残基/寡原子：仅用基础色名，避免 spectrum / util 在极简结构上异常。"""
    from pymol import cmd

    for cname in ("marine", "cyan", "red", "blue", "yellow", "magenta", "green", "orange", "violet"):
        try:
            cmd.color(cname, selection)
            return
        except Exception:
            continue
    logger.warning("PyMOL sticks 未能为 %s 着色", selection)


def _png_file_valid(path: Path) -> bool:
    try:
        if not path.is_file():
            return False
        if path.stat().st_size < 64:
            return False
        with path.open("rb") as f:
            return f.read(8) == b"\x89PNG\r\n\x1a\n"
    except OSError:
        return False


def _pymol_render_png(pdb_path: Path, out_png: Path) -> None:
    from pymol import cmd

    _pymol_ensure_launcher()
    with _pymol_lock:
        try:
            cmd.reinitialize()
            cmd.delete("all")
            cmd.load(str(pdb_path), "struct")
            # 浅灰底：避免 cartoon/sticks 与纯白背景融在一起（用户见「PNG 有效但画面全白」）
            try:
                cmd.set("ray_opaque_background", 1)
            except Exception:
                pass
            cmd.bg_color("grey90")

            try:
                n_atom = int(cmd.count_atoms("struct"))
            except Exception:
                n_atom = 999
            try:
                n_ca = int(cmd.count_atoms("struct and name CA"))
            except Exception:
                n_ca = 0
            # 模板用单残基片段时 cartoon 往往极细或不可见，改 lines+sticks
            use_sticks = n_atom <= 50 or n_ca < 2

            cmd.hide("everything", "struct")
            if use_sticks:
                cmd.show("lines", "struct")
                cmd.show("sticks", "struct")
                try:
                    cmd.set("stick_radius", 0.14)
                except Exception:
                    pass
            else:
                try:
                    cmd.show_as("cartoon", "struct")
                except Exception as e:
                    logger.warning("PyMOL cartoon 失败，回退 lines+sticks: %s", e)
                    cmd.show("lines", "struct")
                    cmd.show("sticks", "struct")
                    use_sticks = True

            try:
                if use_sticks:
                    _pymol_color_sticks_simple("struct")
                else:
                    _pymol_color_cartoon_robust("struct")
            except Exception as e:
                logger.warning("PyMOL 着色失败，仍尝试导出 PNG: %s", e)
            cmd.orient("struct")
            cmd.zoom("struct", buffer=4.0)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            try:
                if out_png.is_file():
                    out_png.unlink()
            except OSError:
                pass
            cmd.png(str(out_png), width=1400, height=1000, dpi=120, ray=0)
            if not _png_file_valid(out_png):
                raise RuntimeError("PyMOL 未写出有效 PNG（文件缺失、过小或非 PNG 头）")
        finally:
            try:
                cmd.reinitialize()
            except Exception:
                pass


@app.post("/api/pymol")
async def run_pymol(req: ToolRequest):
    try:
        p = Path(req.file_path).expanduser()
        if not p.is_file():
            raise HTTPException(status_code=400, detail="结构文件不存在或路径无效，请确认已通过上传获得有效路径。")
        run_id = uuid.uuid4().hex[:16]
        out_dir = _pyskills_out_dir() / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_png = out_dir / "structure_cartoon.png"
        try:
            _pymol_render_png(p, out_png)
        except Exception as e:
            logger.exception("PyMOL render")
            raise HTTPException(
                status_code=500,
                detail=f"PyMOL 渲染失败: {e}",
            ) from e
        if not out_png.is_file():
            raise HTTPException(status_code=500, detail="PyMOL 未生成图像文件。")
        return {
            "status": "success",
            "message": "蛋白质三维结构渲染完成。",
            "data": {
                "image_url": _public_url(run_id, "structure_cartoon.png"),
            },
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("run_pymol")
        raise HTTPException(status_code=500, detail=f"PyMOL 服务异常: {e}") from e


@app.post("/api/sascore")
async def run_sascore(req: ToolRequest):
    try:
        p = Path(req.file_path)
        if not p.is_file():
            raise HTTPException(status_code=400, detail=f"文件不存在: {req.file_path}")
        raw = _read_text_safe(p)
        try:
            smiles = _first_smiles_from_text(raw)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e)) from e
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail=f"无法解析 SMILES: {smiles!r}")
        try:
            score = float(_get_calculate_score()(mol))
        except Exception as e:
            logger.exception("SA_Score calculateScore")
            raise HTTPException(status_code=500, detail=f"SA_Score 计算失败: {e}") from e
        return {"status": "success", "data": {"sa_score": score}}
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("run_sascore")
        raise HTTPException(status_code=500, detail=traceback.format_exc()) from e


class LipinskiRequest(BaseModel):
    """与 sascore 一致：file_path 指向含 SMILES 的文本文件（共享 uploads 卷）。"""

    file_path: str


@app.post("/api/lipinski")
async def run_lipinski(req: LipinskiRequest):
    """Lipinski Rule of Five：MW、LogP、HBD、HBA 与 is_druglike（违反 ≤1 条为 true）。"""
    try:
        p = Path(req.file_path).expanduser()
        if not p.is_file():
            raise HTTPException(status_code=400, detail=f"文件不存在: {req.file_path}")
        raw = _read_text_safe(p)
        smiles = _first_smiles_from_text(raw)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail=f"无法解析 SMILES: {smiles!r}")
        mw = float(round(Descriptors.MolWt(mol), 4))
        logp = float(round(Descriptors.MolLogP(mol), 4))
        hbd = int(RDKitLipinski.NumHDonors(mol))
        hba = int(RDKitLipinski.NumHAcceptors(mol))
        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        is_druglike = violations <= 1
        return {
            "status": "success",
            "data": {
                "MW": mw,
                "logP": logp,
                "HBD": hbd,
                "HBA": hba,
                "lipinski_violations": violations,
                "is_druglike": is_druglike,
            },
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("run_lipinski")
        raise HTTPException(status_code=500, detail=traceback.format_exc()) from e


class DrugSimilarityRequest(BaseModel):
    file_path: str
    top_n: int = Field(default=10, ge=1, le=50)
    fingerprint: str = Field(default="morgan_ecfp")
    metric: str = Field(default="tanimoto")
    databases: list[str] = Field(default_factory=lambda: ["pubchem", "chembl"])
    # 与 third_party/drug-similarity/find_similar.py → get_similar_from_all 默认对齐
    similarity_threshold: float = Field(default=0.5, ge=0.05, le=0.99)
    max_results_per_db: int = Field(default=100, ge=10, le=500)


@app.post("/api/drug_similarity")
async def run_drug_similarity(req: DrugSimilarityRequest):
    """
    调用技能包 scripts/drug_similarity.py（需联网访问公共化学库）。
    产物写入 /app/uploads/results/drug_similarity/<run_id>/similarity_report.html 与 similarity_results.csv。
    """
    script = DRUG_SIM_BUNDLE / "scripts" / "drug_similarity.py"
    if not script.is_file():
        raise HTTPException(
            status_code=503,
            detail="药物相似性技能包未挂载：缺少 /app/third_party/drug-similarity/scripts/drug_similarity.py（请检查 compose 卷）。",
        )
    try:
        p = Path(req.file_path).expanduser()
        if not p.is_file():
            raise HTTPException(status_code=400, detail=f"文件不存在: {req.file_path}")
        raw = _read_text_safe(p)
        smiles = _first_smiles_from_text(raw)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail=f"无法解析 SMILES: {smiles!r}")
        canon = Chem.MolToSmiles(mol)
        dbs = [str(d).lower().strip() for d in (req.databases or []) if str(d).strip()]
        if not dbs:
            dbs = ["pubchem", "chembl"]
        allowed_db = {"pubchem", "chembl", "drugbank", "bindingdb", "zinc", "chemspider"}
        dbs = [d for d in dbs if d in allowed_db] or ["pubchem", "chembl"]
        run_id = uuid.uuid4().hex[:16]
        out_dir = Path("/app/uploads/results/drug_similarity") / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_html = out_dir / "similarity_report.html"
        cmd: list[str] = [
            sys.executable,
            str(script),
            "--input",
            canon,
            "--input-type",
            "smiles",
            "--top-n",
            str(int(req.top_n)),
            "--similarity-threshold",
            str(float(req.similarity_threshold)),
            "--max-results-per-db",
            str(int(req.max_results_per_db)),
            "--fingerprint",
            str(req.fingerprint or "morgan_ecfp").lower(),
            "--metric",
            str(req.metric or "tanimoto").lower(),
            "--output-format",
            "html",
            "--output",
            str(out_html),
        ]
        for db in dbs:
            cmd.extend(["--databases", db])
        logger.info("drug_similarity CLI: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd,
            cwd=str(DRUG_SIM_BUNDLE),
            capture_output=True,
            text=True,
            timeout=900,
            env={**os.environ, "PYTHONNOUSERSITE": "1"},
        )
        if proc.returncode != 0:
            err = (proc.stderr or proc.stdout or "").strip() or f"exit {proc.returncode}"
            logger.error(
                "drug_similarity CLI 失败 rc=%s stderr=%s stdout=%s",
                proc.returncode,
                (proc.stderr or "")[:8000],
                (proc.stdout or "")[:4000],
            )
            logger.warning("drug_similarity CLI 失败: %s", err[:2000])
            raise HTTPException(status_code=500, detail=f"drug_similarity 执行失败: {err[:4000]}")
        if not out_html.is_file():
            raise HTTPException(status_code=500, detail="未生成 similarity_report.html")
        public_html = f"/uploads/results/drug_similarity/{run_id}/similarity_report.html"
        out_csv = out_dir / "similarity_results.csv"
        csv_public = f"/uploads/results/drug_similarity/{run_id}/similarity_results.csv"
        data_payload: dict = {
            "html_url": public_html,
            "query_smiles": canon,
            "stdout_tail": (proc.stdout or "")[-4000:],
        }
        # 与 drug_similarity.py 同目录落盘的 Top-N 命中表（pandas/csv）；供前端「下载 CSV」
        if out_csv.is_file():
            data_payload["csv_url"] = csv_public
        return {
            "status": "success",
            "message": "结构相似性检索完成，已生成交互式 HTML 报告（基于公共库检索，结果依赖外网数据可用性）。",
            "data": data_payload,
        }
    except HTTPException:
        raise
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="drug_similarity 执行超时（>900s）") from None
    except Exception as e:
        logger.exception("run_drug_similarity")
        raise HTTPException(status_code=500, detail=traceback.format_exc()) from e


def _read_gene_list_from_table(path: Path) -> list[str]:
    raw = _read_text_safe(path)
    genes: list[str] = []
    suf = path.suffix.lower()
    if suf == ".csv":
        reader = csv.reader(io.StringIO(raw))
        rows = [r for r in reader if any(c.strip() for c in r)]
        if not rows:
            raise ValueError("CSV 为空")
        header = [c.strip().lower() for c in rows[0]]
        if "gene" in header:
            gi = header.index("gene")
            for row in rows[1:]:
                if len(row) > gi and row[gi].strip():
                    genes.append(row[gi].strip())
        else:
            for row in rows:
                if row and row[0].strip() and not str(row[0]).strip().startswith("#"):
                    genes.append(str(row[0]).strip())
    else:
        for line in raw.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s.split()[0].strip())
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    if len(out) < 3:
        raise ValueError("有效基因数量过少（至少 3 个），无法进行富集分析。")
    return out


def _gseapy_plot_enrichment(df, png_path: Path, title: str) -> None:
    """优先 gseapy 内置 barplot，失败则用 Matplotlib 简易条形图。"""
    try:
        if hasattr(gseapy, "barplot"):
            gseapy.barplot(
                df,
                column="Adjusted P-value",
                title=title,
                cutoff=0.05,
                ofname=str(png_path),
            )
        else:
            from gseapy.plot import barplot as plot_barplot

            plot_barplot(
                df,
                column="Adjusted P-value",
                title=title,
                cutoff=0.05,
                ofname=str(png_path),
            )
        if png_path.is_file():
            return
    except Exception as e:
        logger.warning("gseapy barplot 失败，回退 Matplotlib: %s", e)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        col_p = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        col_term = "Term" if "Term" in df.columns else df.columns[0]
        sub = df.copy().head(20)
        if col_p not in sub.columns:
            raise ValueError("结果表缺少 P-value 列")
        sub = sub.assign(_neglog=-np.log10(sub[col_p].clip(lower=1e-300)))
        sub = sub.sort_values("_neglog", ascending=True)
        fig, ax = plt.subplots(figsize=(9, max(4, len(sub) * 0.25)))
        ax.barh(sub[col_term].astype(str), sub["_neglog"], color="#3b6ea5")
        ax.set_xlabel("-log10(adjusted P)" if col_p == "Adjusted P-value" else "-log10(P)")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(str(png_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        logger.exception("Matplotlib 富集图失败")
        raise RuntimeError(f"无法生成富集图: {e}") from e


# Enrichr 不可用时用于 ORA 的极小人类相关基因集（离线、无网络）
_LOCAL_ORA_GENE_SETS_BASE: dict[str, list[str]] = {
    "Synthetic_oncogenic_signaling": [
        "TP53",
        "MYC",
        "KRAS",
        "EGFR",
        "BRAF",
        "MAPK1",
        "MAPK3",
        "STAT3",
        "PIK3CA",
        "AKT1",
        "PTEN",
        "MTOR",
    ],
    "Synthetic_cell_cycle_DNA_repair": [
        "TP53",
        "BRCA1",
        "BRCA2",
        "ATM",
        "CHEK2",
        "E2F1",
        "CCND1",
        "CDK4",
        "RB1",
        "CDKN1A",
        "MDM2",
    ],
    "Synthetic_transcription_apoptosis": [
        "MYC",
        "STAT3",
        "STAT1",
        "JAK2",
        "NFKB1",
        "RELA",
        "BCL2",
        "BAX",
        "CASP3",
        "CYCS",
    ],
}
# 合并集：输入基因只要命中任一子集即可产生 ORA 行，减少「有输入却无重叠」的假阴性
_LOCAL_ORA_GENE_SETS: dict[str, list[str]] = {
    **_LOCAL_ORA_GENE_SETS_BASE,
    "Synthetic_combined_reference": sorted(
        {g for genes in _LOCAL_ORA_GENE_SETS_BASE.values() for g in genes}
    ),
}


def _gseapy_results_nonempty(enr) -> bool:
    try:
        r = getattr(enr, "results", None)
        return r is not None and len(r) > 0
    except Exception:
        return False


def _gseapy_local_ora_enrich(gene_list: list[str]):
    """离线 ORA，不访问 Enrichr；返回与 enrichr 接口类似的 Enrichr 对象（含 .results）。"""
    try:
        return gseapy.enrich(
            gene_list=gene_list,
            gene_sets=_LOCAL_ORA_GENE_SETS,
            organism="Human",
            outdir=None,
            no_plot=True,
            verbose=False,
        )
    except TypeError:
        return gseapy.enrich(
            gene_list=gene_list,
            gene_sets=_LOCAL_ORA_GENE_SETS,
            outdir=None,
            no_plot=True,
        )


@app.post("/api/gseapy")
async def run_gseapy(req: ToolRequest):
    try:
        p = Path(req.file_path).expanduser()
        if not p.is_file():
            raise HTTPException(status_code=400, detail="基因列表文件不存在或路径无效。")
        try:
            gene_list = _read_gene_list_from_table(p)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e)) from e

        gene_set = "KEGG_2021_Human"
        run_id = uuid.uuid4().hex[:16]
        out_dir = _pyskills_out_dir() / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "enrichr_results.csv"
        png_path = out_dir / "enrichr_barplot.png"
        used_backend = "enrichr"
        enr = None

        try:
            enr = gseapy.enrichr(
                gene_list=gene_list,
                description="gseapy_worker",
                gene_sets=gene_set,
                organism="Human",
                outdir=None,
                no_plot=True,
            )
        except Exception as e:
            logger.warning("gseapy.enrichr 异常（将尝试离线 ORA）: %s", e)

        # Enrichr 可能不抛错但返回空表；或网络超时已被库吞掉——统一在无结果时走离线 ORA
        if not _gseapy_results_nonempty(enr):
            if enr is not None:
                logger.warning("Enrichr 无有效富集行，改用离线内置基因集 ORA")
            try:
                enr_local = _gseapy_local_ora_enrich(gene_list)
                if _gseapy_results_nonempty(enr_local):
                    enr = enr_local
                    gene_set = "LOCAL_SYNTHETIC_PATHWAYS_V1"
                    used_backend = "local_ora"
            except Exception as e2:
                logger.exception("离线 ORA 回退失败: %s", e2)

        if not _gseapy_results_nonempty(enr):
            raise HTTPException(
                status_code=502,
                detail=(
                    "基因集富集失败：Enrichr 不可用或未返回结果（常见于内网/防火墙），"
                    "且离线内置基因集也未产生有效富集。请确认基因为人类 HGNC 符号并重试；"
                    "部署侧请更新 worker-pyskills 并保证可访问 maayanlab.cloud 或接受离线结果。"
                ),
            )

        try:
            enr.results.to_csv(csv_path, index=False, encoding="utf-8")
        except Exception as e:
            logger.exception("写入富集 CSV")
            raise HTTPException(status_code=500, detail=f"保存富集结果表失败: {e}") from e

        try:
            title_prefix = "Enrichr" if used_backend == "enrichr" else "本地 ORA"
            _gseapy_plot_enrichment(
                enr.results.head(50), png_path, title=f"{title_prefix} · {gene_set}"
            )
        except Exception as e:
            logger.exception("富集可视化")
            raise HTTPException(status_code=500, detail=str(e)) from e

        try:
            ap_col = enr.results["Adjusted P-value"]
            n_sig = int((ap_col.astype(float) < 0.05).sum())
        except Exception:
            n_sig = 0

        msg = "基因富集分析完成。"
        if used_backend == "local_ora":
            msg += "（Enrichr 不可用时已使用内置小型基因集做离线 ORA，结果仅供演示。）"

        return {
            "status": "success",
            "message": msg,
            "data": {
                "csv_url": _public_url(run_id, "enrichr_results.csv"),
                "image_url": _public_url(run_id, "enrichr_barplot.png"),
                "gene_set_library": gene_set,
                "enrichment_backend": used_backend,
                "input_gene_count": len(gene_list),
                "enriched_term_count": len(enr.results),
                "significant_term_count_approx": n_sig,
            },
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("run_gseapy")
        raise HTTPException(status_code=500, detail=f"GSEApy 服务异常: {e}") from e


from radiomics_api import router as radiomics_router
from dynamic_runner import router as dynamic_runner_router

app.include_router(radiomics_router)
app.include_router(dynamic_runner_router)
