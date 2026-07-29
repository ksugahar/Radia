#!/usr/bin/env python
"""Verify a Radia IH Simulink archive and run its extracted MATLAB smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath


REQUIRED_MEMBERS = {
    "manifest.json",
    "matlab/IH_VERSION",
    "matlab/install_radia_ih.m",
    "matlab/radia_ih.slx",
    "matlab/verify_radia_ih_release.m",
    "matlab/+radia/+simulink/validateIHNativeConfig.m",
    "matlab/radia_ih_eddy_sfun.mexw64",
    "matlab/radia_ih_thermal_sfun.mexw64",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _console_safe(value: str, encoding: str | None = None) -> str:
    target = encoding or sys.stdout.encoding or "utf-8"
    return value.encode(target, errors="backslashreplace").decode(target)


def _safe_member(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def verify_archive(archive: Path) -> dict:
    if not archive.is_file():
        raise FileNotFoundError(f"IH release archive does not exist: {archive}")
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        unsafe = sorted(name for name in names if not _safe_member(name))
        if unsafe:
            raise RuntimeError(f"Unsafe ZIP members: {', '.join(unsafe)}")
        missing = sorted(REQUIRED_MEMBERS - names)
        if missing:
            raise RuntimeError(f"IH release archive is incomplete: {', '.join(missing)}")
        manifest = json.loads(bundle.read("manifest.json"))
        if manifest.get("schema") != "radia.simulink.ih-release-manifest.v1":
            raise RuntimeError("Unsupported IH release manifest schema")
        if manifest.get("release_channel") != "preview":
            raise RuntimeError("The first IH release must declare preview channel")
        if manifest.get("matlab_release") != "R2026a" or \
                manifest.get("platform") != "win64" or \
                manifest.get("mex_extension") != "mexw64":
            raise RuntimeError("IH release runtime compatibility is invalid")
        if manifest.get("backend") != "native-mex-sfunction" or \
                manifest.get("python_fallback") is not False:
            raise RuntimeError("IH release backend contract is invalid")
        if manifest.get("operator_assembly") != "preassembled":
            raise RuntimeError("IH release must declare its assembly boundary")

        declared = set()
        for item in manifest.get("files", []):
            name = item.get("path", "")
            if name in declared or name not in names or not _safe_member(name):
                raise RuntimeError(f"Invalid manifest member: {name}")
            payload = bundle.read(name)
            if len(payload) != item.get("size") or digest(payload) != item.get("sha256"):
                raise RuntimeError(f"Manifest hash or size mismatch: {name}")
            declared.add(name)
        payload_members = names - {"manifest.json"}
        if declared != payload_members:
            drift = sorted(declared ^ payload_members)
            raise RuntimeError(f"Manifest file inventory drift: {', '.join(drift)}")
    return manifest


def run_matlab_smoke(archive: Path, matlab: Path, timeout: int = 300) -> str:
    if not matlab.is_file():
        raise FileNotFoundError(f"MATLAB executable does not exist: {matlab}")
    scratch = Path(r"C:\temp")
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="radia-ih-verify-", dir=scratch) as temporary:
        root = Path(temporary)
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(root)
        matlab_root = str((root / "matlab").resolve()).replace("'", "''")
        expression = (
            f"addpath('{matlab_root}','-begin');"
            "report=verify_radia_ih_release();assert(report.passed);"
        )
        result = subprocess.run(
            [str(matlab), "-batch", expression],
            capture_output=True,
            timeout=timeout,
        )
        # MATLAB emits UTF-8 in batch mode even when the Windows process
        # locale is cp932.  Keep subprocess in binary mode so Python never
        # starts a locale-decoding reader thread that can fail before the
        # ASCII release marker is inspected.
        output = (result.stdout or b"").decode("utf-8", errors="backslashreplace") + \
            (result.stderr or b"").decode("utf-8", errors="backslashreplace")
        if result.returncode != 0:
            raise RuntimeError(
                f"Extracted IH MATLAB smoke failed with exit {result.returncode}:\n{output}"
            )
        if "RADIA_IH_RELEASE_OK" not in output:
            raise RuntimeError("MATLAB smoke did not emit RADIA_IH_RELEASE_OK")
        return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--matlab", type=Path)
    parser.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    archive = args.archive.resolve()
    manifest = verify_archive(archive)
    if not args.manifest_only:
        if args.matlab is None:
            parser.error("--matlab is required unless --manifest-only is used")
        matlab_output = run_matlab_smoke(
            archive, args.matlab.resolve(), args.timeout)
        print(_console_safe(matlab_output.rstrip()))
    print(json.dumps({
        "status": "passed",
        "package": manifest["package"],
        "version": manifest["version"],
        "commit": manifest["commit"],
        "sha256": digest(archive.read_bytes()),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
