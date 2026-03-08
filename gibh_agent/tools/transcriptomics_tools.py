"""
Atomic Transcriptomics tools (pydeseq2).
Bulk RNA-seq differential expression. Requires pydeseq2>=0.4.0.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="run_bulk_deseq2",
    description="Perform differential expression analysis on bulk RNA-seq count matrix using PyDESeq2. Requires counts CSV (genes x samples) and metadata CSV with design factor column.",
    category="Transcriptomics",
    output_type="json",
)
def run_bulk_deseq2(
    counts_csv: str,
    metadata_csv: str,
    design_factor: str,
) -> Dict[str, Any]:
    """
    Run DESeq2-style differential expression via PyDESeq2.

    Args:
        counts_csv: Path to count matrix CSV (rows=genes, columns=samples).
        metadata_csv: Path to sample metadata CSV (must contain design_factor column).
        design_factor: Column name in metadata that defines the experimental groups.

    Returns:
        Dict with status, results table (log2FC, pvalue, padj), and optional error.
    """
    try:
        import pandas as pd
        from pydeseq2.ds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ImportError as e:
        logger.warning("pydeseq2 not installed: %s", e)
        return {"status": "error", "error": "pydeseq2 is required. Install with: pip install pydeseq2>=0.4.0"}
    p_counts = Path(counts_csv)
    p_meta = Path(metadata_csv)
    if not p_counts.is_file():
        return {"status": "error", "error": f"Counts file not found: {counts_csv}"}
    if not p_meta.is_file():
        return {"status": "error", "error": f"Metadata file not found: {metadata_csv}"}
    try:
        counts = pd.read_csv(p_counts, index_col=0)
        metadata = pd.read_csv(p_meta, index_col=0)
        # Ensure sample order matches
        if design_factor not in metadata.columns:
            return {"status": "error", "error": f"Design factor '{design_factor}' not in metadata columns: {list(metadata.columns)}"}
        # Align counts to metadata
        inter = counts.columns.intersection(metadata.index)
        if len(inter) == 0:
            return {"status": "error", "error": "No overlapping samples between counts and metadata."}
        counts = counts[inter].loc[:, metadata.loc[inter].index]
        metadata = metadata.loc[counts.columns]
        inference = DefaultInference()
        dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors=[design_factor], inference=inference)
        dds.deseq2()
        stat = DeseqStats(dds)
        res = stat.summary()
        if res is None or res.empty:
            return {"status": "error", "error": "DESeq2 produced no results (check design and counts)."}
        out_df = res[["log2FoldChange", "pvalue", "padj"]].copy()
        out_df = out_df.rename(columns={"log2FoldChange": "log2FC", "padj": "adj_pvalue"})
        out_df = out_df.dropna(subset=["adj_pvalue"])
        return {
            "status": "success",
            "results": out_df.head(5000).to_dict(orient="index"),
            "design_factor": design_factor,
            "n_genes": len(res),
            "counts_path": str(p_counts.resolve()),
            "metadata_path": str(p_meta.resolve()),
        }
    except Exception as e:
        logger.exception("run_bulk_deseq2 failed: %s", e)
        return {"status": "error", "error": str(e)}
