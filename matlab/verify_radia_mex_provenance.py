"""Verify the native build sidecar before MATLAB loads radia_mex."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


SCHEMA = "radia.native-build-provenance.v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(mex: Path, *, checkout: Path | None = None) -> str:
    mex = mex.resolve()
    sidecar = Path(f"{mex}.build.json")
    if not mex.is_file() or not sidecar.is_file():
        raise RuntimeError("radia_mex or its native build provenance is missing; rebuild from clean main")
    data = json.loads(sidecar.read_text(encoding="utf-8"))
    source_commit = data.get("source_commit", "")
    if (data.get("schema") != SCHEMA or data.get("source_dirty") is not False
            or not re.fullmatch(r"[0-9a-f]{40}", source_commit)
            or data.get("binary_name") != mex.name
            or data.get("binary_bytes") != mex.stat().st_size
            or data.get("binary_sha256") != sha256(mex)):
        raise RuntimeError("radia_mex does not match its clean native build provenance")
    if checkout is not None and (checkout / ".git").exists():
        trust = f"safe.directory={checkout.resolve().as_posix()}"
        result = subprocess.run(["git", "-c", trust, "-C", str(checkout),
                                 "rev-parse", "HEAD"],
                                check=True, capture_output=True, text=True)
        if result.stdout.strip().lower() != source_commit:
            raise RuntimeError("radia_mex was built from a different source commit")
        result = subprocess.run(["git", "-c", trust, "-C", str(checkout),
                                 "status", "--porcelain=v1", "--untracked-files=all",
                                 "--", "src", "matlab", "Build.ps1", "CMakeLists.txt",
                                 "tools/native_build_provenance.ps1"],
                                check=True, capture_output=True, text=True)
        if result.stdout.strip():
            raise RuntimeError("native source changed after radia_mex was built")
    return source_commit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mex", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    checkout = root if args.mex.resolve().parent == (root / "matlab").resolve() else None
    print("RADIA_MEX_COMMIT:" + verify(args.mex, checkout=checkout))


if __name__ == "__main__":
    main()
