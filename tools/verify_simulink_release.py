#!/usr/bin/env python
"""Verify a Radia Simulink archive and run its extracted MATLAB smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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
    "matlab/radia_ih_eddy_sfun.m",
    "matlab/radia_ih_thermal_sfun.m",
    "matlab/verify_radia_ih_release.m",
    "matlab/+radia/setup.m",
    "matlab/+radia/+internal/pythonProcessPath.m",
    "matlab/+radia/+simulink/ihEddySFunction.m",
    "matlab/+radia/+simulink/ihThermalSFunction.m",
    "matlab/+radia/+simulink/validateIHNativeConfig.m",
    "matlab/radia_mex.mexw64",
    "matlab/mkl_avx2.2.dll",
    "matlab/mkl_core.2.dll",
    "matlab/mkl_def.2.dll",
    "matlab/mkl_intel_thread.2.dll",
    "matlab/mkl_rt.2.dll",
}
LEGACY_REQUIRED_MEMBERS = {
    "manifest.json",
    "matlab/IH_VERSION",
    "matlab/install_radia_ih.m",
    "matlab/radia_ih.slx",
    "matlab/verify_radia_ih_release.m",
    "matlab/+radia/+simulink/validateIHNativeConfig.m",
    "matlab/radia_ih_eddy_sfun.mexw64",
    "matlab/radia_ih_thermal_sfun.mexw64",
}
FULL_REQUIRED_MEMBERS = {
    "manifest.json",
    "matlab/install_radia_simulink.m",
    "matlab/verify_radia_simulink_release.m",
    "matlab/radia_simulink_library.slx",
    "matlab/radia_ih.slx",
    "matlab/radia_maglev.slx",
    "matlab/radia_streamfunction_optimization.slx",
    "matlab/radia_mex.mexw64",
    "matlab/radia_ih_eddy_sfun.m",
    "matlab/radia_ih_thermal_sfun.m",
    "matlab/+radia/+simulink/ihEddySFunction.m",
    "matlab/+radia/+simulink/ihThermalSFunction.m",
    "matlab/mkl_avx2.2.dll",
    "matlab/mkl_core.2.dll",
    "matlab/mkl_def.2.dll",
    "matlab/mkl_intel_thread.2.dll",
    "matlab/mkl_rt.2.dll",
}
LEGACY_FULL_REQUIRED_MEMBERS = {
    "manifest.json",
    "matlab/install_radia_simulink.m",
    "matlab/verify_radia_simulink_release.m",
    "matlab/radia_simulink_library.slx",
    "matlab/radia_ih.slx",
    "matlab/radia_streamfunction_optimization.slx",
    "matlab/radia_mex.mexw64",
    "matlab/radia_ih_eddy_sfun.mexw64",
    "matlab/radia_ih_thermal_sfun.mexw64",
    "matlab/mkl_avx2.2.dll",
    "matlab/mkl_core.2.dll",
    "matlab/mkl_def.2.dll",
    "matlab/mkl_intel_thread.2.dll",
    "matlab/mkl_rt.2.dll",
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
        raise FileNotFoundError(f"Simulink release archive does not exist: {archive}")
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        unsafe = sorted(name for name in names if not _safe_member(name))
        if unsafe:
            raise RuntimeError(f"Unsafe ZIP members: {', '.join(unsafe)}")
        manifest = json.loads(bundle.read("manifest.json"))
        schema = manifest.get("schema")
        preview_v1 = schema == "radia.simulink.ih-release-manifest.v1"
        preview_v2 = schema == "radia.simulink.ih-release-manifest.v2"
        full_v1 = schema == "radia.simulink.library-release-manifest.v1"
        full_v2 = schema == "radia.simulink.library-release-manifest.v2"
        if preview_v1:
            required_members = LEGACY_REQUIRED_MEMBERS
        elif preview_v2:
            required_members = REQUIRED_MEMBERS
        elif full_v1:
            required_members = LEGACY_FULL_REQUIRED_MEMBERS
        elif full_v2:
            required_members = FULL_REQUIRED_MEMBERS
        else:
            raise RuntimeError("Unsupported Simulink release manifest schema")
        if not re.fullmatch(r"[0-9a-f]{40}", str(manifest.get("commit", ""))):
            raise RuntimeError("Simulink release manifest commit is invalid")
        missing = sorted(required_members - names)
        if missing:
            raise RuntimeError(
                f"Simulink release archive is incomplete: {', '.join(missing)}"
            )
        if manifest.get("matlab_release") != "R2026a" or \
                manifest.get("platform") != "win64" or \
                manifest.get("mex_extension") != "mexw64":
            raise RuntimeError("Simulink release runtime compatibility is invalid")
        if manifest.get("required_matlab_products") != ["MATLAB", "Simulink"]:
            raise RuntimeError("Simulink release product requirements are invalid")
        if preview_v1 or preview_v2:
            if manifest.get("release_channel") != "preview":
                raise RuntimeError("The first IH release must declare preview channel")
            expected_backend = (
                "native-mex-sfunction"
                if preview_v1
                else "matlab-level2+radia-mex-handles"
            )
            if manifest.get("backend") != expected_backend or \
                    manifest.get("python_fallback") is not False:
                raise RuntimeError("IH release backend contract is invalid")
            if manifest.get("operator_assembly") != "preassembled":
                raise RuntimeError("IH release must declare its assembly boundary")
            if preview_v2:
                _verify_level2_ih_contract(manifest)
        else:
            if manifest.get("release_channel") != "production":
                raise RuntimeError("The full library must declare production channel")
            if manifest.get("package") != "radia-simulink-library" or \
                    manifest.get("entry_model") != \
                    "matlab/radia_simulink_library.slx":
                raise RuntimeError("The full library entry contract is invalid")
            expected_backend = (
                "native-mex-sfunction-and-mex-handle"
                if full_v1 else "application-specific"
            )
            if manifest.get("backend") != expected_backend or \
                    manifest.get("python_per_step") is not False or \
                    manifest.get("python_fallback_per_step") is not False:
                raise RuntimeError("The full library backend contract is invalid")
            if full_v2:
                _verify_level2_ih_contract(manifest)
            elif set(manifest.get("required_mex", [])) != {
                    "matlab/radia_mex.mexw64",
                    "matlab/radia_ih_eddy_sfun.mexw64",
                    "matlab/radia_ih_thermal_sfun.mexw64",
                }:
                    raise RuntimeError("The full library MEX inventory is invalid")
            toolbox_requirements = manifest.get("feature_toolbox_requirements", {})
            if toolbox_requirements != {
                "adjoint_topology_optimization": ["Optimization Toolbox"],
                "stream_function_topology_optimization": [
                    "Optimization Toolbox"
                ],
            }:
                raise RuntimeError("The full library toolbox contract is invalid")

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


def _verify_level2_ih_contract(manifest: dict) -> None:
    if manifest.get("ih_backend") != "matlab-level2+radia-mex-handles":
        raise RuntimeError("The IH Level-2 backend declaration is invalid")
    if set(manifest.get("required_mex", [])) != {
        "matlab/radia_mex.mexw64",
    }:
        raise RuntimeError("The IH standalone MEX inventory is invalid")
    if set(manifest.get("required_matlab_sfunctions", [])) != {
        "matlab/radia_ih_eddy_sfun.m",
        "matlab/radia_ih_thermal_sfun.m",
        "matlab/+radia/+simulink/ihEddySFunction.m",
        "matlab/+radia/+simulink/ihThermalSFunction.m",
    }:
        raise RuntimeError("The IH Level-2 S-Function inventory is invalid")
    if manifest.get("standalone_mex_debug_api") is not True:
        raise RuntimeError("The IH standalone MEX debug contract is invalid")
    if manifest.get("python_runtime_required_for_native_mex") is not True or \
            manifest.get("python_per_step") is not False:
        raise RuntimeError("The IH native runtime dependency contract is invalid")


def run_matlab_smoke(archive: Path, matlab: Path, timeout: int = 300) -> str:
    if not matlab.is_file():
        raise FileNotFoundError(f"MATLAB executable does not exist: {matlab}")
    scratch = Path(r"C:\temp")
    scratch.mkdir(parents=True, exist_ok=True)
    manifest = verify_archive(archive)
    full_library = (
        manifest.get("schema") in {
            "radia.simulink.library-release-manifest.v1",
            "radia.simulink.library-release-manifest.v2",
        }
    )
    verification_function = (
        "verify_radia_simulink_release"
        if full_library
        else "verify_radia_ih_release"
    )
    success_marker = (
        "RADIA_SIMULINK_RELEASE_OK" if full_library else "RADIA_IH_RELEASE_OK"
    )
    with tempfile.TemporaryDirectory(
            prefix="radia-ih-verify-", dir=scratch) as temporary:
        root = Path(temporary)
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(root)
        matlab_root = str((root / "matlab").resolve()).replace("'", "''")
        expression = (
            f"addpath('{matlab_root}','-begin');"
            f"report={verification_function}();assert(report.passed);"
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
                f"Extracted Simulink MATLAB smoke failed with exit "
                f"{result.returncode}:\n{output}"
            )
        if success_marker not in output:
            raise RuntimeError(
                f"MATLAB smoke did not emit {success_marker}"
            )
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
