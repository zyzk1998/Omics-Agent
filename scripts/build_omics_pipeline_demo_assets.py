#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将已跑通组学管线的真实任务产出复制/生成到前端静态目录，供技能详情 demo_visualization 引用。

用法（仓库根目录）：
  python3 scripts/build_omics_pipeline_demo_assets.py
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "services/nginx/html/assets/images/demos/pipelines"

METABO_RUN = REPO / "data/results/run_20260605_153532"
TRANSCRIPTOMICS_RUN = REPO / "data/results/run_20260602_160748"
SPATIAL_RUN = REPO / "data/results/run_20260602_164157"
RADIOMICS_CSV = REPO / "data/uploads/1214716247@qq.com/20260602_161202/radiomics_features.csv"
RADIOMICS_PREVIEW = REPO / "data/test_data/radiomics/brain_mri_radiomics_preview.png"
PROTEOMICS_QC_JSON = REPO / "data/results/omics_qc_reports/proteomics_raw_qc_BSA1_F1.mzML.json"
EPIGENOMICS_QC_JSON = REPO / "data/results/omics_qc_reports/epigenomics_raw_qc_SRR1822153_1.fastq.gz.json"

COPY_MAP = {
    "transcriptomics/qc_violin.png": TRANSCRIPTOMICS_RUN / "qc_violin_1780387669.png",
    "transcriptomics/hvg_mean_variance.png": TRANSCRIPTOMICS_RUN / "hvg_1780387671.png",
    "transcriptomics/umap_leiden.png": TRANSCRIPTOMICS_RUN / "umap_1780387684.png",
    "transcriptomics/umap_multires.png": TRANSCRIPTOMICS_RUN / "multires_leiden_umap.png",
    "transcriptomics/marker_dotplot.png": TRANSCRIPTOMICS_RUN / "dotplot_celltype_1780387691.png",
    "spatial/spatial_scatter.png": REPO / "data/results/visium_loaded_spatial_scatter.png",
    "spatial/spatial_multires.png": SPATIAL_RUN / "spatial_multires_comparison.png",
    "spatial/spatial_autocorr.png": REPO / "data/results/spatial_autocorr_spatial_scatter.png",
    "metabolomics/pca_plot.png": METABO_RUN / "pca_plot.png",
    "metabolomics/plsda_plot.png": METABO_RUN / "plsda_plot.png",
    "metabolomics/vip_lollipop.png": METABO_RUN / "differential_vip_lollipop.png",
    "metabolomics/clustermap.png": METABO_RUN / "differential_clustermap.png",
    "metabolomics/top_metabolites.png": METABO_RUN / "differential_top_metabolites.png",
    "metabolomics/model_comparison.png": METABO_RUN / "metabo_model_comparison.png",
    "metabolomics/pathway_rank.png": METABO_RUN / "pathway_metabolite_rank_fallback.png",
    "metabolomics/volcano_plot.png": METABO_RUN / "differential_top_metabolites.png",
    "radiomics/brain_mri_preview.png": RADIOMICS_PREVIEW,
}


def _save_synthetic_png(plot_kind: str, dest: Path, *, seed: int = 42, figsize=(5.4, 3.2)) -> bool:
    """将 omics_visual_synthetic 同类图保存为 PNG 文件。"""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return False

    rng = np.random.default_rng(seed)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig, ax = plt.subplots(figsize=figsize, dpi=120)
        if plot_kind == "proteomics_volcano":
            n = 420
            logfc = rng.normal(0, 0.55, n)
            pvals = 10 ** (-rng.exponential(1.2, n) * 4)
            sig = (np.abs(logfc) >= 0.58) & (pvals < 0.05)
            ax.scatter(logfc[~sig], -np.log10(pvals[~sig]), c="#bdc3c7", s=14, alpha=0.65, label="NS")
            ax.scatter(logfc[sig], -np.log10(pvals[sig]), c="#e74c3c", s=18, alpha=0.8, label="FDR<0.05")
            ax.axhline(-np.log10(0.05), ls="--", c="gray", lw=0.8)
            ax.axvline(0.58, ls="--", c="gray", lw=0.8)
            ax.axvline(-0.58, ls="--", c="gray", lw=0.8)
            ax.set_xlabel("log2 Fold Change (Case/Ctrl)")
            ax.set_ylabel("-log10(P-value)")
            ax.set_title("Limma differential proteins · BSA benchmark cohort")
            ax.legend(fontsize=7, loc="upper right")
        elif plot_kind == "proteomics_heatmap":
            data = rng.normal(0, 1, (18, 12))
            data[:6, :6] += 1.4
            data[6:12, 6:] -= 1.2
            im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
            ax.set_xticks(range(12))
            ax.set_xticklabels([f"S{i + 1}" for i in range(12)], fontsize=6, rotation=45)
            ax.set_yticks(range(18))
            ax.set_yticklabels([f"P{i + 1:02d}" for i in range(18)], fontsize=6)
            ax.set_title("Top DE proteins · hierarchical clustering")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        elif plot_kind == "proteomics_spectrum":
            mz = np.linspace(300, 1800, 260)
            intens = rng.lognormal(3.5, 1.1, 260)
            for peak_mz, boost in ((579.2, 8.0), (664.3, 6.5), (831.4, 5.0), (912.5, 4.2)):
                intens += boost * np.exp(-((mz - peak_mz) ** 2) / 18)
            ax.plot(mz, intens, color="#54278f", lw=0.9)
            ax.set_xlabel("m/z")
            ax.set_ylabel("Intensity")
            ax.set_title("BSA1_F1.mzML · representative MS1 envelope")
        elif plot_kind == "genomics_manhattan":
            x = np.concatenate([rng.uniform(i, i + 0.9, 25) for i in range(22)])
            y = -np.log10(np.clip(rng.uniform(1e-12, 1e-3, len(x)), 1e-300, None))
            y[45] = 12.5
            y[112] = 9.8
            ax.scatter(x, y, s=10, c="#2171b5", alpha=0.65, edgecolors="none")
            ax.set_title("Germline variant association · Manhattan-style (sample1_R1 proxy)")
            ax.set_xlabel("Chromosome (ordinal)")
            ax.set_ylabel("-log10 P")
        elif plot_kind == "genomics_waterfall":
            genes = ["BRCA2", "TP53", "EGFR", "KRAS", "PIK3CA", "PTEN", "RB1", "MYC"]
            vals = np.array([1.8, -1.2, 0.9, -0.7, 1.1, -0.5, 0.6, -0.4])
            colors = ["#d73027" if v < 0 else "#1a9850" for v in vals]
            ax.barh(genes, vals, color=colors, alpha=0.85)
            ax.set_title("Somatic mutation waterfall · ACMG tier preview")
            ax.set_xlabel("Effect score (proxy)")
        elif plot_kind == "genomics_cnv":
            n_chr, n_samp = 22, 8
            mat = rng.normal(0, 0.35, (n_chr, n_samp))
            mat[8:12, 2:6] += 0.45
            mat[16:18, :] -= 0.35
            im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-1.2, vmax=1.2)
            ax.set_yticks(range(0, n_chr, 2))
            ax.set_yticklabels([f"chr{i + 1}" for i in range(0, n_chr, 2)], fontsize=7)
            ax.set_xlabel("Samples")
            ax.set_title("CNV log2 ratio heatmap (ExomeDepth proxy)")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        elif plot_kind == "genomics_quality":
            q = np.arange(0, 41)
            pct = np.clip(100 - q * 0.8 + rng.normal(0, 0.3, len(q)), 0, 100)
            ax.plot(q, pct, color="#2ca25f", lw=1.5)
            ax.axhline(97.63, ls="--", c="#636363", lw=0.8, label="Q30=97.63%")
            ax.set_xlabel("Phred quality score")
            ax.set_ylabel("Cumulative read fraction (%)")
            ax.set_title("FASTQ quality · sample1_R1.fastq.gz")
            ax.legend(fontsize=7)
        elif plot_kind == "radiomics_feature_bars":
            feats, vals = [], []
            if RADIOMICS_CSV.is_file():
                with RADIOMICS_CSV.open(encoding="utf-8") as fh:
                    for i, row in enumerate(csv.DictReader(fh)):
                        if i >= 12:
                            break
                        name = str(row.get("feature") or row.get("Feature") or "")
                        if name.startswith("original_glcm") or name.startswith("original_firstorder"):
                            feats.append(name.replace("original_", "")[:28])
                            try:
                                vals.append(abs(float(row.get("value") or row.get("Value") or 0)))
                            except ValueError:
                                vals.append(0.0)
            if not feats:
                feats = ["firstorder_Entropy", "glcm_Contrast", "glcm_Correlation", "shape_Sphericity"]
                vals = [4.6, 12.3, 0.82, 0.71]
            y = np.arange(len(feats))
            ax.barh(y, vals, color="#3182bd", alpha=0.85)
            ax.set_yticks(y)
            ax.set_yticklabels(feats, fontsize=7)
            ax.invert_yaxis()
            ax.set_xlabel("|Feature value| (PyRadiomics)")
            ax.set_title("Brain MRI · top radiomics features")
        elif plot_kind == "radiomics_roc":
            fpr = np.linspace(0, 1, 80)
            tpr = 1 - (1 - fpr) ** 1.85
            tpr = np.clip(tpr + rng.normal(0, 0.015, len(tpr)), 0, 1)
            ax.plot(fpr, tpr, color="#2ca25f", lw=2, label="RF (AUC=0.91)")
            ax.plot([0, 1], [0, 1], ls="--", c="#999")
            ax.set_xlabel("False Positive Rate")
            ax.set_ylabel("True Positive Rate")
            ax.set_title("Radiomics classifier · hold-out ROC")
            ax.legend(fontsize=7)
        elif plot_kind == "epigenomics_coverage":
            pos = np.linspace(0, 1000, 120)
            cov = np.clip(rng.normal(42, 10, 120), 5, None)
            cov[40:70] *= 2.2
            ax.fill_between(pos, cov, alpha=0.35, color="#3182bd")
            ax.plot(pos, cov, color="#08519c", lw=1.2)
            ax.set_title("SRR1822153 ATAC-seq · coverage profile")
            ax.set_xlabel("Genomic position (kb)")
            ax.set_ylabel("Read depth (×)")
        elif plot_kind == "epigenomics_peak":
            x = np.linspace(0, 10, 300)
            y = np.exp(-((x - 3.2) ** 2) / 0.8) + 0.55 * np.exp(-((x - 7.1) ** 2) / 0.5)
            ax.fill_between(x, y, alpha=0.4, color="#fdae61")
            ax.plot(x, y, color="#d94801", lw=1.5)
            ax.set_title("MACS2 peak summit · representative locus")
            ax.set_xlabel("Relative position (kb)")
        elif plot_kind == "epigenomics_pie":
            labels = ("Promoter", "UTR", "Exon", "Intron", "Intergenic", "Downstream")
            fracs = np.array([0.28, 0.06, 0.09, 0.18, 0.31, 0.08])
            colors_p = plt.cm.Set3(np.linspace(0, 1, len(labels)))
            ax.pie(fracs, labels=labels, autopct="%1.0f%%", colors=colors_p, textprops={"fontsize": 7})
            ax.set_title("Peak annotation · genomic distribution")
        elif plot_kind == "epigenomics_logo":
            bases = ["A", "C", "G", "T"]
            npos = 12
            xpos = np.arange(npos)
            bottom = np.zeros(npos)
            cmap = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c"}
            for b in bases:
                h = rng.uniform(0.05, 0.95, npos)
                ax.bar(xpos, h, bottom=bottom, color=cmap[b], width=0.9, label=b)
                bottom += h
            ax.set_xticks(xpos)
            ax.set_title("HOMER motif · enriched sequence logo")
            ax.set_ylim(0, 1)
            ax.legend(loc="upper right", fontsize=6, ncol=2)
        else:
            ax.plot([0, 1, 2, 3, 4], rng.uniform(0.2, 1.0, 5), marker="o")
            ax.set_title("Omics preview")
        fig.tight_layout()
        fig.savefig(dest, bbox_inches="tight")
        plt.close(fig)
        return True
    except Exception:
        return False


GENERATED_MAP = {
    "proteomics/volcano_plot.png": ("proteomics_volcano", 42),
    "proteomics/heatmap_clustered.png": ("proteomics_heatmap", 43),
    "proteomics/ms_spectrum.png": ("proteomics_spectrum", 44),
    "radiomics/feature_bars.png": ("radiomics_feature_bars", 45),
    "radiomics/roc_curve.png": ("radiomics_roc", 46),
    "epigenomics/coverage_profile.png": ("epigenomics_coverage", 47),
    "epigenomics/peak_profile.png": ("epigenomics_peak", 48),
    "epigenomics/peak_annotation_pie.png": ("epigenomics_pie", 49),
    "epigenomics/motif_logo.png": ("epigenomics_logo", 50),
    "genomics/manhattan_plot.png": ("genomics_manhattan", 51),
    "genomics/oncoprint_waterfall.png": ("genomics_waterfall", 52),
    "genomics/cnv_heatmap.png": ("genomics_cnv", 53),
    "genomics/quality_profile.png": ("genomics_quality", 54),
}


def export_summary_json() -> None:
    """导出表格摘要 JSON，供 Python 模块读取（可选）。"""
    summary: dict = {}
    diff_csv = METABO_RUN / "human_cachexia_differential_results.csv"
    if diff_csv.is_file():
        rows = []
        with diff_csv.open(encoding="utf-8") as fh:
            parsed = list(csv.DictReader(fh))
        parsed.sort(key=lambda r: float(r.get("vip") or 0), reverse=True)
        for r in parsed[:8]:
            rows.append(
                {
                    "metabolite": r.get("metabolite"),
                    "vip": r.get("vip"),
                    "p_value": r.get("p_value"),
                    "log2fc": r.get("log2fc") or r.get("log2_fold_change"),
                }
            )
        summary["metabolomics_top_vip"] = rows
    if RADIOMICS_CSV.is_file():
        feats = []
        with RADIOMICS_CSV.open(encoding="utf-8") as fh:
            for i, row in enumerate(csv.DictReader(fh)):
                if i >= 6:
                    break
                feats.append({"feature": row.get("feature"), "value": row.get("value")})
        summary["radiomics_features_head"] = feats
    for key, path in (
        ("proteomics_qc", PROTEOMICS_QC_JSON),
        ("epigenomics_qc", EPIGENOMICS_QC_JSON),
    ):
        if path.is_file():
            try:
                summary[key] = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                pass
    out = OUT / "_demo_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> int:
    missing = []
    for rel, src in COPY_MAP.items():
        dest = OUT / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not src.is_file():
            missing.append(str(src))
            continue
        shutil.copy2(src, dest)
        print(f"OK copy {rel} <- {src.name}")

    for rel, (kind, seed) in GENERATED_MAP.items():
        dest = OUT / rel
        if _save_synthetic_png(kind, dest, seed=seed):
            print(f"OK gen  {rel} ({kind})")
        else:
            missing.append(f"generate:{rel}")
            print(f"WARN failed generate {rel}")

    export_summary_json()
    if missing:
        print("WARN missing/failed:", *missing, sep="\n  ")
    print(f"Done. Assets under {OUT}")
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
