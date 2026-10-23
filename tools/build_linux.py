#!/usr/bin/env python3
"""Build CPU Python extensions on Linux with an explicit Intel LLVM toolchain.

Run with the isolated job environment's Python. Dependency pins come from
pyproject.toml; this entry point neither installs packages nor loads modules.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import tomllib

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ("_radia_pybind", "peec_matrices", "sparsesolv_ngsolve", "axifem")


def check_versions(root=ROOT):
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    versions = {}
    for name in ("ngsolve", "netgen-mesher"):
        requirement = next(r for r in project["project"]["dependencies"]
                           if r.startswith(name + "=="))
        expected = requirement.split("==", 1)[1]
        actual = importlib.metadata.version(name)
        if actual != expected:
            raise RuntimeError(f"{name}: installed {actual}, required {expected}; use a fresh job venv")
        versions[name] = actual
    return versions


def tool(name):
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Missing tool: {name}. Load the intended toolchain first; no compiler fallback is used")
    return executable


def run(command, **kwargs):
    print("+ " + " ".join(map(str, command)), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def source_identity(root):
    pathspec = ["--", ".", ":(exclude)packages/cubit-mesh-export/src/cubit_mesh_export/cubit_mesh_export.ccm"]
    def git(*args):
        return run(["git", "-C", root, *args], capture_output=True,
                   text=True, encoding="utf-8").stdout.rstrip("\n")
    sha = git("rev-parse", "HEAD").strip().lower()
    if not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise RuntimeError(f"Invalid source SHA: {sha!r}")
    status = git("status", "--porcelain=v1", "--untracked-files=all", *pathspec)
    diff = git("diff", "--binary", "HEAD", *pathspec)
    return {"source_commit": sha, "source_dirty": bool(status),
            "source_change_fingerprint_sha256": hashlib.sha256(
                (status + "\0" + diff).encode("utf-8")).hexdigest()}


def write_provenance(binary, identity):
    if identity["source_dirty"]:
        return False
    manifest = {"schema": "radia.native-build-provenance.v1", **identity,
                "binary_name": binary.name, "binary_bytes": binary.stat().st_size,
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    Path(str(binary) + ".build.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True,
                        help="job-owned build directory, separate from other toolchains")
    parser.add_argument("--jobs", type=int, default=1,
                        help="admitted parallel build jobs (default: 1)")
    parser.add_argument("--cc", default=os.environ.get("CC", "icx"))
    parser.add_argument("--cxx", default=os.environ.get("CXX", "icpx"))
    parser.add_argument("--targets", nargs="+", choices=TARGETS, default=list(TARGETS))
    parser.add_argument("--require-clean-source", action="store_true",
                        help="reject dirty source before configuring an acceptance build")
    args = parser.parse_args(argv)
    if sys.platform != "linux":
        parser.error("Linux only; use Build.ps1 on Windows")
    if sys.version_info[:2] != (3, 12):
        parser.error("Use Python 3.12 in an isolated job venv")
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    build = args.build_dir.resolve()
    if build == ROOT or build.is_relative_to(ROOT / "src"):
        parser.error("--build-dir must not replace the source directory")
    build.mkdir(parents=True, exist_ok=True)
    report_path = build / "linux-build.json"
    report = {"schema": "radia.linux-build.v1", "passed": False,
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "source": str(ROOT), "python": sys.executable,
              "platform": platform.platform(), "jobs": args.jobs,
              "targets": args.targets, "artifacts": {},
              "environment": {k: os.environ.get(k) for k in
                              ("MKLROOT", "MKL_THREADING_LAYER", "OMP_NUM_THREADS", "NGS_NUM_THREADS")},
              "affinity": sorted(os.sched_getaffinity(0))}
    try:
        start = source_identity(ROOT)
        report["source_identity"] = start
        if args.require_clean_source and start["source_dirty"]:
            raise RuntimeError("Clean source is required for this build")
        report["versions"] = check_versions()
        cc, cxx, cmake = tool(args.cc), tool(args.cxx), tool("cmake")
        tool("ninja")
        compiler_version = run([cxx, "--version"], capture_output=True, text=True).stdout
        c_version = run([cc, "--version"], capture_output=True, text=True).stdout
        if not all(re.search(r"Intel|clang", version, flags=re.I)
                   for version in (c_version, compiler_version)):
            raise RuntimeError("Select Intel LLVM icx/icpx or Clang explicitly; GCC is not this build entry point's toolchain")
        report["compiler"] = {"cc": cc, "cxx": cxx, "version": compiler_version,
                              "c_version": c_version}
        mkl = Path(os.environ.get("MKLROOT", ""))
        if not os.environ.get("MKLROOT") or not (mkl / "include/mkl.h").is_file() or not (mkl / "lib/libmkl_rt.so").is_file():
            raise RuntimeError("Set MKLROOT to oneMKL containing include/mkl.h and lib/libmkl_rt.so")
        for target in args.targets:
            Path(str(ROOT / "src/radia" / f"{target}.so") + ".build.json").unlink(missing_ok=True)
        run([cmake, "-S", ROOT, "-B", build, "-G", "Ninja",
             "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_C_COMPILER={cc}",
             f"-DCMAKE_CXX_COMPILER={cxx}", f"-DPython3_EXECUTABLE={sys.executable}"])
        run([cmake, "--build", build, "--target", *args.targets, "--parallel", args.jobs])
        # Import from this checkout, never an unrelated installed Radia.
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONNOUSERSITE="1")
        modules = ["radia" if t == "_radia_pybind" else f"radia.{t}" for t in args.targets]
        run([sys.executable, "-c", "import importlib; " +
             "; ".join(f"print(importlib.import_module({m!r}).__file__)" for m in modules)],
            env=env, cwd=build)
        for target in args.targets:
            artifact = ROOT / "src/radia" / f"{target}.so"
            report["artifacts"][target] = {"path": str(artifact),
                                          "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
                                          "dependencies": run(["ldd", artifact], capture_output=True, text=True).stdout}
        if source_identity(ROOT) != start:
            raise RuntimeError("Source identity changed during native build; provenance was not issued")
        for target in args.targets:
            artifact = ROOT / "src/radia" / f"{target}.so"
            report["artifacts"][target]["provenance_issued"] = write_provenance(artifact, start)
        if start["source_dirty"]:
            print("Dirty source checkout: native acceptance manifests were not issued")
        report["passed"] = True
    except (OSError, RuntimeError, subprocess.CalledProcessError, importlib.metadata.PackageNotFoundError) as error:
        report["error"] = str(error)
        print(f"Build failed: {error}", file=sys.stderr)
    finally:
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Build report: {report_path}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
