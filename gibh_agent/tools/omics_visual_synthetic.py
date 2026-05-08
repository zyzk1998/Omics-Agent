"""
秒级合成示意图（matplotlib → PNG base64 data URI），用于组学步骤可视化兜底。

无 matplotlib/numpy 或绘图失败时返回空字符串，由调用方跳过。
"""
from __future__ import annotations

import base64
import io
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def synthetic_png_data_uri(plot_kind: str, seed: int = 42) -> str:
    """
    plot_kind: manhattan | waterfall | coverage | generic | spectra | peak |
               heatmap | bubble | pie | logo | cnv_heatmap
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        logger.debug("omics_visual_synthetic skip (deps): %s", exc)
        return ""

    rng = np.random.default_rng(seed)

    try:
        fig, ax = plt.subplots(figsize=(5.2, 3.0), dpi=110)

        if plot_kind == "manhattan":
            x = np.concatenate([rng.uniform(i, i + 0.9, 25) for i in range(22)])
            y = -np.log10(np.clip(rng.uniform(1e-12, 1e-3, len(x)), 1e-300, None))
            ax.scatter(x, y, s=10, c="#2171b5", alpha=0.65, edgecolors="none")
            ax.set_title("Synthetic GWAS / Manhattan-style preview")
            ax.set_xlabel("Chromosome (ordinal)")
            ax.set_ylabel("-log10 P")

        elif plot_kind == "waterfall":
            genes = [f"G{i:03d}" for i in range(1, 17)]
            vals = np.sort(rng.normal(0, 1.2, len(genes)))
            colors = ["#d73027" if v < 0 else "#1a9850" for v in vals]
            ax.barh(genes, vals, color=colors, alpha=0.85)
            ax.set_title("Synthetic waterfall / gene ranking preview")
            ax.set_xlabel("Effect (arbitrary units)")

        elif plot_kind == "coverage":
            pos = np.linspace(0, 1000, 120)
            cov = np.clip(rng.normal(35, 8, 120), 5, None)
            ax.fill_between(pos, cov, alpha=0.35, color="#3182bd")
            ax.plot(pos, cov, color="#08519c", lw=1.2)
            ax.set_title("Synthetic coverage profile")
            ax.set_xlabel("Position (kb)")
            ax.set_ylabel("Depth (×)")

        elif plot_kind == "spectra":
            mz = np.linspace(300, 1200, 200)
            intens = rng.lognormal(4, 1.2, 200)
            intens[mz < 400] *= 3
            ax.plot(mz, intens, color="#54278f", lw=0.9)
            ax.set_title("Synthetic MS spectrum envelope")
            ax.set_xlabel("m/z")
            ax.set_ylabel("Intensity (a.u.)")

        elif plot_kind == "peak":
            x = np.linspace(0, 10, 300)
            y = np.exp(-((x - 3.2) ** 2) / 0.8) + 0.4 * np.exp(-((x - 7.1) ** 2) / 0.5)
            ax.fill_between(x, y, alpha=0.4, color="#fdae61")
            ax.plot(x, y, color="#d94801", lw=1.5)
            ax.set_title("Synthetic ChIP/ATAC peak profile")

        elif plot_kind == "heatmap":
            data = rng.uniform(-2, 2, (16, 24))
            im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
            ax.set_title("Synthetic clustered heatmap (DE / abundance)")
            ax.set_xlabel("Samples")
            ax.set_ylabel("Features")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        elif plot_kind == "bubble":
            np_cat = 14
            terms = [f"GO:{440 + i:04d}" for i in range(np_cat)]
            logp = rng.uniform(2, 12, np_cat)
            gene_counts = rng.integers(8, 220, np_cat)
            sizes = 30 + gene_counts * 0.8
            logp_max = max(float(np.max(logp)), 1e-6)
            colors = plt.cm.viridis(logp / logp_max)
            ax.scatter(logp, range(np_cat), s=sizes, c=colors, alpha=0.85, edgecolors="none")
            ax.set_yticks(range(np_cat))
            ax.set_yticklabels(terms, fontsize=6)
            ax.set_xlabel("-log10(adj.P)")
            ax.set_title("Synthetic GO/KEGG bubble / dot plot")

        elif plot_kind == "pie":
            labels = ("Promoter", "UTR", "Exon", "Intron", "Intergenic", "Downstream")
            fracs = np.array([0.28, 0.06, 0.09, 0.18, 0.31, 0.08])
            colors_p = plt.cm.Set3(np.linspace(0, 1, len(labels)))
            ax.pie(fracs, labels=labels, autopct="%1.0f%%", colors=colors_p, textprops={"fontsize": 7})
            ax.set_title("Peak genomic distribution (placeholder)")

        elif plot_kind == "logo":
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
            ax.set_title("Synthetic sequence logo (Motif)")
            ax.set_ylim(0, 1)
            ax.legend(loc="upper right", fontsize=6, ncol=2)

        elif plot_kind == "cnv_heatmap":
            n_chr, n_samp = 22, 10
            mat = rng.normal(0, 0.35, (n_chr, n_samp))
            mat[8:12, 3:7] += 0.45
            im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-1.2, vmax=1.2)
            ax.set_yticks(range(0, n_chr, 2))
            ax.set_yticklabels([f"chr{i + 1}" for i in range(0, n_chr, 2)], fontsize=7)
            ax.set_xlabel("Samples")
            ax.set_title("Synthetic CNV profile heatmap (log2 ratio)")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        else:
            ax.plot([0, 1, 2, 3, 4], rng.uniform(0.2, 1.0, 5), marker="o")
            ax.set_title("Synthetic omics preview")

        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        b64 = base64.standard_b64encode(buf.getvalue()).decode("ascii")
        return f"data:image/png;base64,{b64}"
    except Exception as exc:
        logger.debug("synthetic_png_data_uri failed: %s", exc)
        return ""
