#!/usr/bin/env python
"""Build a self-verifying IH preview or full Radia Simulink release."""

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
FULL_REQUIRED_MEX = (*REQUIRED_MEX, "radia_mex.mexw64")
FULL_REQUIRED_MODELS = (
    "radia_simulink_library.slx",
    "radia_ih.slx",
    "radia_streamfunction_optimization.slx",
)
FULL_RUNTIME_DLLS = (
    "mkl_avx2.2.dll",
    "mkl_core.2.dll",
    "mkl_def.2.dll",
    "mkl_intel_thread.2.dll",
    "mkl_rt.2.dll",
)
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
        raise FileNotFoundError(f"Required native release binary is missing: {path}")
    if path.stat().st_size < 1024:
        raise RuntimeError(f"Native release binary is unexpectedly small: {path}")
    with path.open("rb") as stream:
        header = stream.read(64)
        if header[:2] != b"MZ":
            raise RuntimeError(f"Native release binary is not a Windows PE file: {path}")
        pe_offset = int.from_bytes(header[60:64], "little")
        stream.seek(pe_offset)
        if stream.read(4) != b"PE\0\0":
            raise RuntimeError(f"Native release binary has no PE signature: {path}")
        machine = int.from_bytes(stream.read(2), "little")
        if machine != 0x8664:
            raise RuntimeError(f"Native release binary is not Windows x64: {path}")


def release_matlab_files() -> tuple[Path, ...]:
    result = subprocess.run(
        [
            "git",
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            "matlab",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    files = tuple(
        Path(line.strip())
        for line in result.stdout.splitlines()
        if line.strip()
    )
    if not files or any(not (ROOT / path).is_file() for path in files):
        raise RuntimeError("Cannot enumerate the MATLAB release surface")
    return files


def build_package(
    mex_dir: Path,
    output_dir: Path,
    *,
    full_library: bool = False,
) -> tuple[Path, Path]:
    matlab_source = ROOT / "matlab"
    required_models = FULL_REQUIRED_MODELS if full_library else REQUIRED_MODELS
    required_mex = FULL_REQUIRED_MEX if full_library else REQUIRED_MEX
    for model in required_models:
        if not (matlab_source / model).is_file():
            raise FileNotFoundError(f"Required Simulink model is missing: {model}")
    for name in required_mex:
        validate_native_binary(mex_dir / name)
    if full_library:
        for name in FULL_RUNTIME_DLLS:
            validate_native_binary(mex_dir / name)

    output_dir.mkdir(parents=True, exist_ok=True)
    package_version = radia_version() if full_library else version()
    package_name = (
        "radia-simulink-library" if full_library else "radia-ih-simulink"
    )
    archive = output_dir / f"{package_name}-v{package_version}.zip"
    sums = output_dir / "SHA256SUMS.txt"
    external_manifest = output_dir / "manifest.json"

    with tempfile.TemporaryDirectory(prefix="radia-ih-release-", dir=output_dir) as tmp:
        stage = Path(tmp)
        if full_library:
            source_files = release_matlab_files()
        else:
            source_files = tuple(Path("matlab") / path for path in PACKAGE_FILES)
        for relative in source_files:
            source = ROOT / relative
            if not source.is_file():
                raise FileNotFoundError(f"Release support file is missing: {relative}")
            destination = stage / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        shutil.copy2(ROOT / "LICENSE", stage / "LICENSE")
        for name in required_mex:
            shutil.copy2(mex_dir / name, stage / "matlab" / name)
        if full_library:
            for name in FULL_RUNTIME_DLLS:
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
            "schema": (
                "radia.simulink.library-release-manifest.v1"
                if full_library
                else "radia.simulink.ih-release-manifest.v1"
            ),
            "package": package_name,
            "version": package_version,
            "radia_version": radia_version(),
            "commit": commit(),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "release_channel": "production" if full_library else "preview",
            "matlab_release": "R2026a",
            "platform": "win64",
            "mex_extension": "mexw64",
            "backend": (
                "native-mex-sfunction-and-mex-handle"
                if full_library
                else "native-mex-sfunction"
            ),
            "operator_assembly": (
                "application-specific" if full_library else "preassembled"
            ),
            "entry_model": (
                "matlab/radia_simulink_library.slx"
                if full_library
                else "matlab/radia_ih.slx"
            ),
            "required_mex": [f"matlab/{name}" for name in required_mex],
            "required_runtime_dll": (
                [f"matlab/{name}" for name in FULL_RUNTIME_DLLS]
                if full_library
                else []
            ),
            "python_per_step": False,
            "python_runtime_required_for_headless_application_blocks": full_library,
            "dt_order": "eddy;transport(theta_prev,theta_now);thermal",
            "files": files,
        }
        if full_library:
            manifest["application_batch_backend"] = (
                "python-headless-or-native-as-declared-by-block"
            )
            manifest["python_fallback_per_step"] = False
            manifest["library_groups"] = [
                "Applications",
                "Material Models",
                "LTspice",
                "Coupling",
                "Reduced Models",
                "Optimization",
            ]
            manifest["verification_entry"] = (
                "matlab/verify_radia_simulink_release.m"
            )
        else:
            manifest["python_fallback"] = False
        manifest_text = json.dumps(manifest, indent=2) + "\n"
        (stage / "manifest.json").write_text(manifest_text, encoding="utf-8")
        external_manifest.write_text(manifest_text, encoding="utf-8")

        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    bundle.write(path, path.relative_to(stage).as_posix())

    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        required = {
            "manifest.json",
            *(f"matlab/{name}" for name in required_mex),
            *(f"matlab/{name}" for name in required_models),
        }
        if full_library:
            required.update(f"matlab/{name}" for name in FULL_RUNTIME_DLLS)
            required.update(
                {
                    "matlab/install_radia_simulink.m",
                    "matlab/verify_radia_simulink_release.m",
                }
            )
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
    parser.add_argument("--full-library", action="store_true")
    args = parser.parse_args()
    archive, sums = build_package(
        args.mex_dir.resolve(),
        args.output_dir.resolve(),
        full_library=args.full_library,
    )
    print(archive)
    print(args.output_dir.resolve() / "manifest.json")
    print(sums)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
