"""
Extended file handlers for modality-specific detection (Spatial, Radiomics, Genomics).
"""
from .extended_handlers import (
    SpatialVisiumHandler,
    RadiomicsHandler,
    GenomicsHandler,
)

__all__ = ["SpatialVisiumHandler", "RadiomicsHandler", "GenomicsHandler"]
