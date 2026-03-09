#!/usr/bin/env python3
"""
Download PyRadiomics official test data (Brain MRI) and convert NRRD -> NIfTI (.nii.gz).
Saves to /app/data/test_data/radiomics/ for pipeline testing.

Usage (from project root):
  PYTHONPATH=. python scripts/download_radiomics_demo.py

Output:
  /app/data/test_data/radiomics/brain_mri.nii.gz
  /app/data/test_data/radiomics/brain_mask.nii.gz
"""
import os
import sys
import urllib.request
from pathlib import Path

try:
    import SimpleITK as sitk
except ImportError:
    print("Error: SimpleITK is required. pip install SimpleITK", file=sys.stderr)
    sys.exit(1)

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = Path(os.getenv("RADIOMICS_DEMO_DIR", "/app/data/test_data/radiomics"))

URLS = {
    "image": "https://github.com/AIM-Harvard/pyradiomics/raw/master/data/brain1_image.nrrd",
    "mask": "https://github.com/AIM-Harvard/pyradiomics/raw/master/data/brain1_label.nrrd",
}
OUTPUT_NAMES = {
    "image": "brain_mri.nii.gz",
    "mask": "brain_mask.nii.gz",
}


def download_file(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "GIBH-Radiomics-Demo/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        dest.write_bytes(resp.read())
    print(f"  -> {dest.stat().st_size} bytes")


def nrrd_to_nii_gz(nrrd_path: Path, nii_path: Path) -> None:
    img = sitk.ReadImage(str(nrrd_path))
    sitk.WriteImage(img, str(nii_path), useCompression=True)
    print(f"Converted {nrrd_path.name} -> {nii_path.name}")


def main():
    # Allow running from host: use local path if /app/data not writable
    if not os.access("/app/data", os.W_OK) and os.access(BASE_DIR / "data", os.W_OK):
        output_dir = BASE_DIR / "data" / "test_data" / "radiomics"
    else:
        output_dir = Path(OUTPUT_DIR)
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = output_dir / ".tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    for key in ("image", "mask"):
        url = URLS[key]
        out_name = OUTPUT_NAMES[key]
        nrrd_path = tmp_dir / url.split("/")[-1]
        nii_path = output_dir / out_name

        if nii_path.exists():
            print(f"Already exists: {nii_path}, skip.")
            continue
        download_file(url, nrrd_path)
        nrrd_to_nii_gz(nrrd_path, nii_path)

    # Cleanup temp NRRD
    for f in tmp_dir.glob("*.nrrd"):
        try:
            f.unlink()
        except OSError:
            pass
    if tmp_dir.exists() and not any(tmp_dir.iterdir()):
        tmp_dir.rmdir()

    print("\nDone. Outputs:")
    print(f"  Image: {output_dir / OUTPUT_NAMES['image']}")
    print(f"  Mask:  {output_dir / OUTPUT_NAMES['mask']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
