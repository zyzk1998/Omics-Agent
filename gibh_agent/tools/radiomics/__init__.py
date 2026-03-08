"""
Radiomics (medical imaging) tools — PyRadiomics, SimpleITK.
Package-based structure: io, analysis, visualization, modeling.

Runtime self-healing: if pyradiomics is missing, attempt pip install (default then --user for non-root).
"""
import logging
import subprocess
import sys

logger = logging.getLogger(__name__)

def _ensure_radiomics_deps():
    try:
        import radiomics  # noqa: F401
        return
    except ImportError:
        pass
    logger.warning("PyRadiomics not found. Attempting auto-installation...")
    pkgs = ["pyradiomics>=3.0.1", "SimpleITK>=2.3.0"]
    for use_user in (False, True):
        try:
            cmd = [sys.executable, "-m", "pip", "install", "-q"] + (["--user"] if use_user else []) + pkgs
            result = subprocess.run(
                cmd,
                timeout=120,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                stderr = (result.stderr or "").strip() or result.stdout or ""
                logger.warning("pip install %s failed: %s", "--user" if use_user else "default", stderr[:500])
                continue
            logger.info("PyRadiomics installed successfully (%s).", "user" if use_user else "system")
            return
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            logger.warning("Auto-install error: %s. Rebuild image with no-cache or run: pip install pyradiomics", e)
            return
    logger.warning(
        "Auto-install failed (no write permission?). Rebuild with no-cache: ./monitor-lite.sh -> 4) 重建并重启"
    )


_ensure_radiomics_deps()

from . import io
from . import analysis
from . import visualization
from . import modeling

__all__ = ["io", "analysis", "visualization", "modeling"]
