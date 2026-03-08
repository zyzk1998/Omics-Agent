"""
Atomic Genomics tools (scikit-allel).
VCF conversion and genotype PCA. Requires scikit-allel>=1.3.5.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="vcf_to_dataframe",
    description="Convert a VCF file to a pandas DataFrame (variants x samples genotype matrix and variant annotations).",
    category="Genomics",
    output_type="json",
)
def vcf_to_dataframe(vcf_path: str) -> Dict[str, Any]:
    """
    Convert VCF to pandas DataFrame.

    Args:
        vcf_path: Path to the VCF file.

    Returns:
        Dict with status, dataframe (as records for JSON), shape, and optional error.
    """
    try:
        import allel
        import pandas as pd
    except ImportError as e:
        logger.warning("scikit-allel not installed: %s", e)
        return {"status": "error", "error": "scikit-allel is required. Install with: pip install scikit-allel>=1.3.5"}
    p = Path(vcf_path)
    if not p.is_file():
        return {"status": "error", "error": f"VCF file not found: {vcf_path}"}
    try:
        # Read VCF; extract CHROM, POS, and genotype matrix
        callset = allel.read_vcf(str(p), fields=["CHROM", "POS", "calldata/GT"])
        if not callset:
            return {"status": "error", "error": "VCF produced no data"}
        chrom = callset.get("variants/CHROM")
        pos = callset.get("variants/POS")
        gt = callset.get("calldata/GT")
        # Genotype array: (variants, samples, ploidy) -> count alternate alleles per sample
        import numpy as np
        if gt is None:
            return {"status": "error", "error": "No genotype (GT) data in VCF"}
        gn = allel.GenotypeArray(gt)
        # Alternate allele count per variant per sample: shape (n_variants, n_samples)
        ac = gn.to_n_alt(fill=-1)
        n_variants, n_samples = ac.shape
        # Build variant table
        df_variants = pd.DataFrame({"CHROM": chrom, "POS": pos})
        # Sample names from VCF header if available
        samples = callset.get("samples")
        if samples is not None and len(samples) == n_samples:
            df_geno = pd.DataFrame(ac, columns=[str(s) for s in samples])
        else:
            df_geno = pd.DataFrame(ac, columns=[f"S{i}" for i in range(n_samples)])
        df = pd.concat([df_variants, df_geno], axis=1)
        out = {
            "status": "success",
            "dataframe": df.head(1000).to_dict(orient="records"),
            "shape": {"rows": n_variants, "cols": len(df.columns)},
            "vcf_path": str(p.resolve()),
        }
        return out
    except Exception as e:
        logger.exception("vcf_to_dataframe failed: %s", vcf_path)
        return {"status": "error", "error": str(e), "vcf_path": vcf_path}


@registry.register(
    name="calculate_genomics_pca",
    description="Perform PCA on genotype data from a VCF file. Returns sample coordinates and explained variance.",
    category="Genomics",
    output_type="json",
)
def calculate_genomics_pca(
    vcf_path: str,
    n_components: int = 10,
    scaling: Optional[str] = "patterson",
) -> Dict[str, Any]:
    """
    Perform PCA on genotype data from VCF.

    Args:
        vcf_path: Path to the VCF file.
        n_components: Number of principal components (default 10).
        scaling: 'patterson', 'standard', or None (default 'patterson').

    Returns:
        Dict with status, pca_coordinates, explained_variance_ratio, and optional error.
    """
    try:
        import allel
        import numpy as np
        import pandas as pd
    except ImportError as e:
        logger.warning("scikit-allel not installed: %s", e)
        return {"status": "error", "error": "scikit-allel is required. Install with: pip install scikit-allel>=1.3.5"}
    p = Path(vcf_path)
    if not p.is_file():
        return {"status": "error", "error": f"VCF file not found: {vcf_path}"}
    try:
        callset = allel.read_vcf(str(p), fields=["calldata/GT", "samples"])
        gt = callset.get("calldata/GT")
        if gt is None:
            return {"status": "error", "error": "No genotype (GT) data in VCF"}
        gn = allel.GenotypeArray(gt)
        # Alternate allele count matrix (n_variants, n_samples)
        ac = gn.to_n_alt(fill=-1)
        # Replace -1 (missing) with column mean for PCA
        mask = ac < 0
        if mask.any():
            col_mean = np.ma.masked_array(ac, mask=mask).mean(axis=0)
            ac = ac.astype(float)
            for j in range(ac.shape[1]):
                ac[ac[:, j] < 0, j] = col_mean[j]
        n_samples = ac.shape[1]
        n_comp = min(n_components, n_samples - 1, ac.shape[0])
        if n_comp < 1:
            return {"status": "error", "error": "Not enough samples or variants for PCA"}
        # PCA (transpose so rows = samples, cols = variants)
        coords, model = allel.pca(ac.T, n_components=n_comp, scaling=scaling or "patterson")
        ev = model.explained_variance_ratio_
        samples = callset.get("samples")
        if samples is not None and len(samples) == n_samples:
            sample_ids = [str(s) for s in samples]
        else:
            sample_ids = [f"S{i}" for i in range(n_samples)]
        df = pd.DataFrame(coords, index=sample_ids, columns=[f"PC{i+1}" for i in range(n_comp)])
        return {
            "status": "success",
            "pca_coordinates": df.to_dict(orient="index"),
            "explained_variance_ratio": ev.tolist(),
            "n_components": n_comp,
            "vcf_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("calculate_genomics_pca failed: %s", vcf_path)
        return {"status": "error", "error": str(e), "vcf_path": vcf_path}
