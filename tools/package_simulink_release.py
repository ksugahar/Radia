#!/usr/bin/env python
"""Build a self-verifying Radia IH Simulink release candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_MEX = ("radia_ih_eddy_sfun.mexw64", "radia_ih_thermal_sfun.mexw64")
REQUIRED_MODELS = ("radia_ih.slx",)
PACKAGE_FILES = (
    "IH_RELEASE.md",
    "IH_VERSION",
    "install_radia_ih.m",
    "radia_ih.slx",
    "verify_radia_ih_release.m",
    "+radia/+simulink/buildIHNativeModel.m",
    "+radia/+simulink/configureIHNativeModel.m",
    "+radia/+simulink/loadIHNativeConfig.m",
    "+radia/+simulink/makeIHNativeConfig.m",
    "+radia/+simulink/makeIHNativeSmokeConfig.m",
    "+radia/+simulink/openIH.m",
    "+radia/+simulink/requireIHNativeRuntime.m",
    "+radia/+simulink/validateIHNativeConfig.m",
    "+radia/+simulink/validateVolFiles.m",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def version() -> str:
    value = (ROOT / "matlab" / "IH_VERSION").read_text(encoding="ascii").strip()
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(?:-[0-9A-Za-z.-]+)?", value):
        raise RuntimeError("matlab/IH_VERSION is not a valid release version")
    return value


def radia_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not match:
        raise RuntimeError("Cannot read Radia version from pyproject.toml")
    return match.group(1)


def commit() -> str:
    github_sha = os.environ.get("GITHUB_SHA", "").strip()
    if re.fullmatch(r"[0-9a-fA-F]{40}", github_sha):
        return github_sha.lower()
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        raise RuntimeError(
            "Cannot determine the release source commit; set GITHUB_SHA or "
            "run the packager inside a Git worktree"
        ) from error
    value = result.stdout.strip()
    if not re.fullmatch(r"[0-9a-fA-F]{40}", value):
        raise RuntimeError("git rev-parse HEAD returned an invalid commit")
    return value.lower()


def validate_native_binary(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required native IH binary is missing: {path}")
    if path.stat().st_size < 1024:
        raise RuntimeError(f"Native IH binary is unexpectedly small: {path}")
    with path.open("rb") as stream:
        header = stream.read(64)
        if header[:2] != b"MZ":
            raise RuntimeError(f"Native IH binary is not a Windows PE file: {path}")
        pe_offset = int.from_bytes(header[60:64], "little")
        stream.seek(pe_offset)
        if stream.read(4) != b"PE\0\0":
            raise RuntimeError(f"Native IH binary has no PE signature: {path}")
        machine = int.from_bytes(stream.read(2), "little")
        if machine != 0x8664:
            raise RuntimeError(f"Native IH binary is not Windows x64: {path}")


def build_package(mex_dir: Path, output_dir: Path) -> tuple[Path, Path]:
    matlab_source = ROOT / "matlab"
    for model in REQUIRED_MODELS:
        if not (matlab_source / model).is_file():
            raise FileNotFoundError(f"Required Simulink model is missing: {model}")
    for name in REQUIRED_MEX:
        validate_native_binary(mex_dir / name)

    output_dir.mkdir(parents=True, exist_ok=True)
    package_version = version()
    archive = output_dir / f"radia-ih-simulink-v{package_version}.zip"
    sums = output_dir / "SHA256SUMS.txt"
    external_manifest = output_dir / "manifest.json"

    with tempfile.TemporaryDirectory(prefix="radia-ih-release-", dir=output_dir) as tmp:
        stage = Path(tmp)
        for relative in PACKAGE_FILES:
            source = matlab_source / relative
            if not source.is_file():
                raise FileNotFoundError(
                    f"Required IH release support file is missing: {relative}"
                )
            destination = stage / "matlab" / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        shutil.copy2(ROOT / "LICENSE", stage / "LICENSE")
        for name in REQUIRED_MEX:
            shutil.copy2(mex_dir / name, stage / "matlab" / name)

        files = []
        for path in sorted(stage.rglob("*")):
            if path.is_file():
                files.append({
                    "path": path.relative_to(stage).as_posix(),
                    "size": path.stat().st_size,
                    "sha256": sha256(path),
                })
        manifest = {
            "schema": "radia.simulink.ih-release-manifest.v1",
            "package": "radia-ih-simulink",
            "version": package_version,
            "radia_version": radia_version(),
            "commit": commit(),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "release_channel": "preview",
            "matlab_release": "R2026a",
            "platform": "win64",
            "mex_extension": "mexw64",
            "backend": "native-mex-sfunction",
            "operator_assembly": "preassembled",
            "entry_model": "matlab/radia_ih.slx",
            "required_mex": [f"matlab/{name}" for name in REQUIRED_MEX],
            "python_fallback": False,
            "dt_order": "eddy;transport(theta_prev,theta_now);thermal",
            "files": files,
        }
        manifest_text = json.dumps(manifest, indent=2) + "\n"
        (stage / "manifest.json").write_text(manifest_text, encoding="utf-8")
        external_manifest.write_text(manifest_text, encoding="utf-8")

        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    bundle.write(path, path.relative_to(stage).as_posix())

    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        required = {"manifest.json", *(f"matlab/{name}" for name in REQUIRED_MEX),
                    *(f"matlab/{name}" for name in REQUIRED_MODELS)}
        missing = sorted(required - names)
        if missing:
            raise RuntimeError(f"Release archive is incomplete: {', '.join(missing)}")
    sums.write_text(
        f"{sha256(archive)}  {archive.name}\n"
        f"{sha256(external_manifest)}  {external_manifest.name}\n",
        encoding="ascii",
    )
    return archive, sums


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mex-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    archive, sums = build_package(args.mex_dir.resolve(), args.output_dir.resolve())
    print(archive)
    print(args.output_dir.resolve() / "manifest.json")
    print(sums)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
