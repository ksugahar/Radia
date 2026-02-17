"""
Create a wheel package for Radia without bundled MKL DLLs.

MKL is declared as a pip dependency (Requires-Dist: mkl>=2024.2.0)
instead of being bundled. This reduces the wheel from ~90 MB to ~5 MB.

Usage:
    python create_wheel.py

The wheel will be created in dist/ directory.
"""

import os
import sys
import hashlib
import base64
import zipfile
from pathlib import Path

# Configuration
SRC_DIR = Path("src/radia")
OUTPUT_DIR = Path("dist")

VERSION = "2.0.0"
PACKAGE_NAME = "radia"
PYTHON_TAG = f"cp{sys.version_info.major}{sys.version_info.minor}"
ABI_TAG = PYTHON_TAG
PLATFORM_TAG = "win_amd64"

WHEEL_NAME = f"{PACKAGE_NAME}-{VERSION}-{PYTHON_TAG}-{ABI_TAG}-{PLATFORM_TAG}.whl"

# DLLs to exclude (MKL and Intel runtime - installed via pip dependency)
EXCLUDED_DLLS = {
    "mkl_rt", "mkl_core", "mkl_intel_thread", "mkl_avx2", "mkl_vml_avx2",
    "mkl_def", "mkl_vml_def", "libiomp5md", "libmmd", "svml_dispmd",
}


def sha256_digest(data: bytes) -> str:
    """Compute base64url-encoded SHA-256 digest for RECORD."""
    h = hashlib.sha256(data)
    return "sha256=" + base64.urlsafe_b64encode(h.digest()).decode("ascii").rstrip("=")


def is_excluded_dll(filename: str) -> bool:
    """Check if a DLL should be excluded from the wheel."""
    name_lower = filename.lower()
    if not name_lower.endswith(".dll"):
        return False
    stem = Path(filename).stem.lower()
    # Check against excluded patterns (handle versioned names like mkl_rt.2)
    for excluded in EXCLUDED_DLLS:
        if stem.startswith(excluded.lower()):
            return True
    return False


def collect_files():
    """Collect all files to include in the wheel."""
    files = {}

    if not SRC_DIR.exists():
        raise FileNotFoundError(f"Source directory not found: {SRC_DIR}")

    for f in SRC_DIR.rglob("*"):
        if not f.is_file():
            continue
        if "__pycache__" in str(f):
            continue
        if is_excluded_dll(f.name):
            print(f"  Excluding DLL: {f.name}")
            continue

        # Archive path: radia/xxx (relative to src/)
        rel = f.relative_to(SRC_DIR.parent)
        archive_path = str(rel).replace("\\", "/")
        files[archive_path] = f

    return files


def create_metadata():
    """Create METADATA content."""
    return f"""Metadata-Version: 2.1
Name: {PACKAGE_NAME}
Version: {VERSION}
Summary: Radia 3D Magnetostatics with PEEC, NGSolve Integration, Intel MKL
Home-page: https://github.com/ksugahar/Radia
Author: Oleg Chubar, Pascal Elleaume
License: BSD-style AND MIT
Requires-Python: >=3.12
Requires-Dist: numpy>=1.20
Requires-Dist: mkl>=2024.2.0
Requires-Dist: intel-cmplr-lib-rt>=2024.2.0
Classifier: Development Status :: 5 - Production/Stable
Classifier: Intended Audience :: Science/Research
Classifier: Programming Language :: Python :: 3.12
Classifier: Programming Language :: C++
Classifier: Topic :: Scientific/Engineering :: Physics
Classifier: Operating System :: Microsoft :: Windows
"""


def create_wheel_metadata():
    """Create WHEEL content."""
    return f"""Wheel-Version: 1.0
Generator: create_wheel.py
Root-Is-Purelib: false
Tag: {PYTHON_TAG}-{ABI_TAG}-{PLATFORM_TAG}
"""


def create_top_level():
    """Create top_level.txt content."""
    return "radia\n"


def build_wheel():
    """Build the wheel file."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    wheel_path = OUTPUT_DIR / WHEEL_NAME
    dist_info = f"{PACKAGE_NAME}-{VERSION}.dist-info"

    print(f"Collecting files from {SRC_DIR}...")
    files = collect_files()

    # Count file types
    pyd_count = sum(1 for p in files if p.endswith(".pyd"))
    py_count = sum(1 for p in files if p.endswith(".py"))
    dll_count = sum(1 for p in files if p.endswith(".dll"))
    other_count = len(files) - pyd_count - py_count - dll_count
    print(f"  {pyd_count} .pyd files, {py_count} .py files, {dll_count} .dll files, {other_count} other")

    # Create metadata
    metadata_content = create_metadata().encode("utf-8")
    wheel_content = create_wheel_metadata().encode("utf-8")
    top_level_content = create_top_level().encode("utf-8")

    print(f"Creating wheel: {wheel_path}")
    record_lines = []

    with zipfile.ZipFile(wheel_path, "w", zipfile.ZIP_DEFLATED) as whl:
        # Write all collected files
        for archive_path, source in sorted(files.items()):
            data = source.read_bytes()
            whl.writestr(archive_path, data)
            digest = sha256_digest(data)
            record_lines.append(f"{archive_path},{digest},{len(data)}")

        # Write dist-info files
        for name, content in [
            ("METADATA", metadata_content),
            ("WHEEL", wheel_content),
            ("top_level.txt", top_level_content),
        ]:
            path = f"{dist_info}/{name}"
            whl.writestr(path, content)
            digest = sha256_digest(content)
            record_lines.append(f"{path},{digest},{len(content)}")

        # Write RECORD (itself has no hash)
        record_path = f"{dist_info}/RECORD"
        record_lines.append(f"{record_path},,")
        record_content = "\n".join(record_lines) + "\n"
        whl.writestr(record_path, record_content)

    size_mb = wheel_path.stat().st_size / (1024 * 1024)
    print(f"\nWheel created: {wheel_path}")
    print(f"Size: {size_mb:.1f} MB")
    print(f"Tag: {PYTHON_TAG}-{ABI_TAG}-{PLATFORM_TAG}")
    print(f"Files: {len(files) + 4}")  # +4 for dist-info files

    # Verify no MKL DLLs leaked into the wheel
    with zipfile.ZipFile(wheel_path, "r") as whl:
        for info in whl.infolist():
            if is_excluded_dll(info.filename.split("/")[-1]):
                print(f"\nWARNING: Excluded DLL found in wheel: {info.filename}")

    return wheel_path


if __name__ == "__main__":
    wheel_path = build_wheel()
    print(f"\nTo install: pip install {wheel_path}")
