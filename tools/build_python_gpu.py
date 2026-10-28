"""Build optional Radia pybind CUDA modules from one shared numerical core.

No packages are installed, no CUDA context is created, and no numerical test
is run. Windows uses MSVC 14.44 and MATLAB-bundled CUDA 12.8; Linux uses an
explicit CUDA 12.8 nvcc and host compiler. The runtime remains pip-owned.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import sysconfig
import time
import tomllib

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "src/core/gpu"
SHARED_SOURCES = ("csr_pcg_host.hpp", "csr_pcg_cuda.cuh", "csr_pcg_python.cpp",
           "csr_pcg_python_cuda.cu", "ams_amg_kernels.cuh", "resident_amg_host.hpp",
           "resident_amg_cuda.cu", "resident_amg_python.cpp", "snapshot.hpp",
           "ams_core.cuh")
SOURCES = SHARED_SOURCES + ("ams_python.cpp", "ams_python_cuda.cu")
COMPONENTS = (("_gpu_csr", "csr_pcg_python.cpp", "csr_pcg_python_cuda.cu"),
              ("_gpu_resident_amg", "resident_amg_python.cpp", "resident_amg_cuda.cu"),
              ("_gpu_ams", "ams_python.cpp", "ams_python_cuda.cu"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def capture(command):
    return subprocess.check_output(list(map(str, command)), text=True, encoding="utf-8", errors="replace").strip()


def identity():
    commit = capture(["git", "-C", ROOT, "rev-parse", "HEAD"])
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise RuntimeError("Invalid source commit")
    if capture(["git", "-C", ROOT, "status", "--porcelain", "--untracked-files=all"]):
        raise RuntimeError("GPU build requires a clean exact-source checkout")
    return commit


def check_toolchain(cuda_root, host, architecture, windows):
    nvcc = cuda_root / ("bin/nvcc.exe" if windows else "bin/nvcc")
    version = capture([nvcc, "--version"])
    if not re.search(r"release 12\.8,", version):
        raise RuntimeError("CUDA toolkit must be 12.8; no compiler override or fallback")
    if architecture not in ({"75"} if windows else {"100"}):
        raise RuntimeError("This acceptance build requires Windows sm_75 or Linux sm_100")
    if windows:
        compiler = subprocess.run([str(host)], capture_output=True, text=True)
        host_version = compiler.stdout + compiler.stderr
        if not re.search(r"19\.44\.", host_version):
            raise RuntimeError("MSVC 14.44 compiler environment is required")
    else:
        host_version = capture([host, "--version"])
        if not re.search(r"(?:g\+\+|gcc|GCC).*11\.5", host_version):
            raise RuntimeError("Explicit reviewed GCC 11.5 host compiler required")
    for key in ("NVCC_PREPEND_FLAGS", "NVCC_APPEND_FLAGS", "NVCC_CCBIN"):
        if os.environ.get(key):
            raise RuntimeError(f"Unreviewed compiler environment: {key}")
    return nvcc, version, host_version


def headers():
    # The installed NGSolve/Netgen ABI owns pybind headers.
    import netgen
    import ngsolve
    requirements = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    for package in ("ngsolve", "netgen-mesher"):
        expected = next(r.split("==", 1)[1] for r in requirements if r.startswith(package + "=="))
        if metadata.version(package) != expected:
            raise RuntimeError(f"{package} must match project metadata: {expected}")
    include = Path(netgen.__file__).resolve().parent / "include"
    if not (include / "pybind11/pybind11.h").is_file():
        raise RuntimeError("Installed Netgen pybind development headers missing")
    return include


def admission(path, digest, commit):
    if sha(path) != digest:
        raise RuntimeError("Admission hash mismatch")
    record = json.loads(path.read_text(encoding="utf-8-sig"))
    if record.get("source_sha") != commit or not isinstance(record.get("owner"), str) or not record.get("owner"):
        raise RuntimeError("Admission must name this exact source and owner")
    if record.get("gpu_build_approved") is not True:
        raise RuntimeError("Individual GPU build admission missing")
    return record


def execute(command, out, label, commands):
    commands.append(list(map(str, command)))
    with (out / (label + ".log")).open("w", encoding="utf-8") as log:
        subprocess.run(list(map(str, command)), cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cuda-root", type=Path, required=True)
    parser.add_argument("--host-compiler", type=Path, required=True)
    parser.add_argument("--architecture", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--matlab-root", type=Path)
    args = parser.parse_args()
    commit = identity()
    record = admission(args.admission, args.admission_sha256, commit)
    windows = platform.system() == "Windows"
    if not windows and platform.system() != "Linux":
        raise RuntimeError("Unsupported build OS")
    cuda = args.cuda_root.resolve()
    if windows and (args.matlab_root is None or cuda != (args.matlab_root.resolve() / "sys/cuda/win64/cuda")):
        raise RuntimeError("Windows must use the explicitly named MATLAB-bundled CUDA")
    nvcc, nvcc_version, host_version = check_toolchain(cuda, args.host_compiler, args.architecture, windows)
    include = headers()
    out = args.output.resolve()
    if out == ROOT or ROOT in out.parents:
        raise RuntimeError("Build output must be outside the source checkout")
    out.mkdir(parents=True, exist_ok=False)
    source_hashes = {str((CORE / name).relative_to(ROOT)).replace("\\", "/"): sha(CORE / name) for name in SOURCES}
    source_hashes["tools/build_python_gpu.py"] = sha(__file__)
    header_hashes = {p.relative_to(include).as_posix(): sha(p) for p in (include / "pybind11").rglob("*") if p.is_file()}
    commands = []
    for module, cpp, cu in COMPONENTS:
        binary = out / (module + sysconfig.get_config_var("EXT_SUFFIX"))
        obj = out / (module + (".obj" if windows else ".o"))
        cuda_command = [nvcc, "-c", "-std=c++17", "-O2", "-ccbin", args.host_compiler,
                        "-gencode=arch=compute_" + args.architecture + ",code=sm_" + args.architecture,
                        "-I" + str(CORE), "-o", obj, CORE / cu]
        if not windows:
            cuda_command += ["-Xcompiler", "-fPIC"]
        execute(cuda_command, out, module + "-cuda", commands)
        if windows:
            command = [args.host_compiler, "/nologo", "/LD", "/EHsc", "/MD", "/O2", "/fp:precise", "/std:c++20", "/utf-8",
                       "/I" + str(include), "/I" + sysconfig.get_path("include"),
                       "/Fo" + str(out / (module + "-host.obj")), CORE / cpp, obj,
                       "/link", "/OUT:" + str(binary), "/IMPLIB:" + str(out / (module + ".lib")),
                       Path(sys.base_prefix) / "libs" / f"python{sys.version_info.major}{sys.version_info.minor}.lib",
                       cuda / "lib/x64/cudart.lib"]
            if module in ("_gpu_resident_amg", "_gpu_ams"):
                command += [cuda / "lib/x64/cublas.lib"]
        else:
            # Link against exact pip runtime files; never embed a site-toolkit rpath.
            spec = __import__("importlib.util", fromlist=["util"])
            runtime_spec = spec.spec_from_file_location("radia_build_runtime", ROOT / "src/radia/_cuda_runtime.py")
            runtime = spec.module_from_spec(runtime_spec)
            runtime_spec.loader.exec_module(runtime)
            libraries = runtime._library_paths(runtime._pins())
            command = [args.host_compiler, "-shared", "-fPIC", "-O2", "-std=c++20", "-I" + str(include),
                       "-I" + sysconfig.get_path("include"), CORE / cpp, obj, "-o", binary,
                       libraries["libcudart.so.12"]["path"]]
            if module in ("_gpu_resident_amg", "_gpu_ams"):
                command += [libraries["libcublas.so.12"]["path"]]
        execute(command, out, module + "-host", commands)
        if identity() != commit or any(sha(ROOT / name) != digest for name, digest in source_hashes.items()):
            raise RuntimeError("Source changed during build")
        private_receipt = Path(str(binary) + ".private-build.json")
        private_receipt.write_text(json.dumps(dict(source_commit=commit, binary_sha256=sha(binary),
            admission=record, commands=commands[:], source_root=str(ROOT), cuda_root=str(cuda),
            host_compiler_path=str(args.host_compiler), pybind_include=str(include)), indent=2) + "\n", encoding="utf-8")
        receipt = dict(schema="radia.native-build-provenance.v1", source_commit=commit, source_dirty=False,
                       binary_name=binary.name, binary_sha256=sha(binary), source_hashes=source_hashes,
                       pybind_header_hashes=header_hashes, python=sys.version, platform=platform.platform(),
                       nvcc=nvcc_version, host_compiler=host_version, architecture=args.architecture,
                       private_receipt_sha256=sha(private_receipt), runtime_import_checked=False,
                       generated_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        Path(str(binary) + ".build.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    # Operational originals stay private; wheels select only the minimal receipt.
    (out / "admission-private.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
