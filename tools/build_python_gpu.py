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
HDIV = CORE / "hdiv"
SHARED_SOURCES = ("csr_pcg_host.hpp", "csr_pcg_cuda.cuh", "csr_pcg_python.cpp",
           "csr_pcg_python_cuda.cu", "ams_amg_kernels.cuh", "resident_amg_host.hpp",
           "resident_amg_cuda.cu", "resident_amg_python.cpp", "snapshot.hpp",
           "ams_core.cuh")
SOURCES = SHARED_SOURCES + ("ams_python.cpp", "ams_python_cuda.cu")
COMPONENTS = (("_gpu_csr", "csr_pcg_python.cpp", "csr_pcg_python_cuda.cu"),
              ("_gpu_resident_amg", "resident_amg_python.cpp", "resident_amg_cuda.cu"),
              ("_gpu_ams", "ams_python.cpp", "ams_python_cuda.cu"),
              ("_hdiv_cuda", "hdiv/hdiv_resident_python.cpp", "hdiv/hdiv_resident_python_cuda.cu"))


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


def check_toolchain(cuda_root, host, architectures, windows):
    nvcc = cuda_root / ("bin/nvcc.exe" if windows else "bin/nvcc")
    version = capture([nvcc, "--version"])
    if not re.search(r"release 12\.8,", version):
        raise RuntimeError("CUDA toolkit must be 12.8; no compiler override or fallback")
    expected = ["75", "80"] if windows else ["100"]
    if list(architectures) != expected:
        raise RuntimeError("Explicit architectures must be Windows sm_75+sm_80 or Linux sm_100")
    if windows:
        compiler = subprocess.run([str(host)], capture_output=True, text=True, encoding="utf-8", errors="replace")
        host_version = compiler.stdout + compiler.stderr
        if not re.search(r"19\.44\.", host_version):
            raise RuntimeError("MSVC 14.44 compiler environment is required")
    else:
        host_version = capture([host, "--version"])
        if not re.search(r"(?:g\+\+|gcc|GCC).*11\.5", host_version):
            raise RuntimeError("Explicit reviewed GCC 11.5 host compiler required")
    for key in ("NVCC_PREPEND_FLAGS", "NVCC_APPEND_FLAGS", "NVCC_CCBIN", "CL", "_CL_", "CFLAGS", "CXXFLAGS", "CPPFLAGS", "CUDAFLAGS"):
        if os.environ.get(key):
            raise RuntimeError(f"Unreviewed compiler environment: {key}")
    return nvcc, version, host_version


def windows_cuda_layout(cuda_root, matlab_root, nvcc_version):
    """Bind the bundled SDK and explicitly locate its relocated libdevice.

    Keep the profile's working CICC/include configuration; profile bypass is
    not required by this source change and needs actual compile evidence.
    Returned paths/hashes belong only in the private build receipt.
    """
    cuda_root, matlab_root = Path(cuda_root).resolve(), Path(matlab_root).resolve()
    if cuda_root != (matlab_root / "sys/cuda/win64/cuda").resolve():
        raise RuntimeError("Windows CUDA SDK must be the explicit MATLAB-bundled root")
    files = {"nvcc": cuda_root / "bin/nvcc.exe",
             "profile": cuda_root / "bin/nvcc.profile",
             "cicc": cuda_root / "nvvm/bin/cicc.exe",
             "libdevice": matlab_root / "bin/win64/libdevice.10.bc",
             "cuda_header": cuda_root / "include/cuda.h"}
    for name, path in files.items():
        if not path.is_file():
            raise RuntimeError(f"Missing MATLAB-bundled CUDA SDK input: {name}")
    release = re.search(r"release (\d+)\.(\d+),", nvcc_version)
    header = re.search(r"^\s*#\s*define\s+CUDA_VERSION\s+(\d+)\b",
                       files["cuda_header"].read_text(encoding="utf-8"), re.MULTILINE)
    if not release or not header:
        raise RuntimeError("Cannot determine bundled CUDA compiler/header ABI")
    expected = int(release[1])*1000 + int(release[2])*10
    if int(header[1]) != expected:
        raise RuntimeError("Bundled CUDA header does not match the selected nvcc release")
    flags = ["--libdevice-directory", str(files["libdevice"].parent)]
    identity = {"header_version": int(header[1]),
                "files": {name: {"path": str(path), "sha256": sha(path)}
                          for name, path in files.items()}}
    return flags, identity


def headers():
    # The installed NGSolve/Netgen ABI owns pybind headers.
    import netgen
    import ngsolve
    requirements = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    for package in ("ngsolve", "netgen-mesher"):
        expected = next(r.split("==", 1)[1] for r in requirements if r.startswith(package + "=="))
        if metadata.version(package) != expected:
            raise RuntimeError(f"{package} must match project metadata: {expected}")
    candidates = {Path(netgen.__file__).resolve().parent / "include"}
    dist = metadata.distribution("netgen-mesher")
    for item in dist.files or ():
        if str(item).replace("\\", "/").endswith("/pybind11/pybind11.h"):
            candidates.add(Path(dist.locate_file(item)).resolve().parents[1])
    found = {p.resolve() for p in candidates if (p / "pybind11/pybind11.h").is_file()}
    if len(found) != 1:
        raise RuntimeError("Missing or ambiguous installed Netgen pybind development headers")
    return found.pop()


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


def installed_config(package, filename):
    """Resolve exact installed SDK metadata, never a site-global build prefix."""
    dist = metadata.distribution(package)
    matches = {Path(dist.locate_file(f)).resolve() for f in dist.files or ()
               if Path(str(f)).name == filename and Path(dist.locate_file(f)).is_file()}
    if len(matches) != 1:
        raise RuntimeError(f"Missing or ambiguous {package} CMake package: {filename}")
    return matches.pop()


def cuda_libraries(cuda, windows):
    if windows:
        return {name: str(cuda / "lib/x64" / (name + ".lib"))
                for name in ("cudart", "cublas", "cusolver")}
    import importlib.util
    spec = importlib.util.spec_from_file_location("radia_build_runtime", ROOT / "src/radia/_cuda_runtime.py")
    runtime = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runtime)
    libraries = runtime._library_paths(runtime._pins())
    return {key: libraries[name]["path"] for key, name in
            [("cudart", "libcudart.so.12"), ("cublas", "libcublas.so.12"), ("cusolver", "libcusolver.so.11")]}


def hdiv_source_hashes():
    """Bind nested CMake fragments as well as the numerical/host sources."""
    return {p.relative_to(ROOT).as_posix(): sha(p)
            for p in HDIV.rglob("*") if p.is_file()}


def build_hdiv_host(args, binary, cuda_object, include, cuda, out, commands):
    """NGSolve owns host compile definitions/includes/link ABI; nvcc stays separate."""
    import shutil
    cmake = shutil.which("cmake")
    if not cmake:
        raise RuntimeError("CMake is required for NGSolve-owned HDiv host staging")
    ng = installed_config("ngsolve", "NGSolveConfig.cmake")
    netgen = installed_config("netgen-mesher", "NetgenConfig.cmake")
    libs = cuda_libraries(cuda, platform.system() == "Windows")
    build = out / "hdiv-host"
    # Imported SDK macros may evaluate these values as CMake source again.
    # Backslashes in Windows paths then become escapes (for example, \h).
    # Keep subprocess argv separate and render CMake paths with forward slashes.
    cmake_path = lambda value: Path(value).as_posix()
    configure = [cmake, "-S", cmake_path(HDIV), "-B", cmake_path(build), "-G", "Ninja",
                 "-DCMAKE_CXX_COMPILER=" + cmake_path(args.host_compiler), "-DCMAKE_BUILD_TYPE=Release",
                 "-DNGSolve_DIR=" + cmake_path(ng.parent), "-DNetgen_DIR=" + cmake_path(netgen.parent),
                 "-DPython3_EXECUTABLE=" + cmake_path(sys.executable), "-DRADIA_PYBIND_INCLUDE=" + cmake_path(include),
                 "-DRADIA_CUDA_OBJECT=" + cmake_path(cuda_object), "-DRADIA_PYTHON_SUFFIX=" + sysconfig.get_config_var("EXT_SUFFIX"),
                 "-DRADIA_OUTPUT=" + cmake_path(out)]
    configure += ["-DRADIA_" + key.upper() + "=" + cmake_path(path) for key, path in libs.items()]
    execute(configure, out, "_hdiv_cuda-host-configure", commands)
    execute([cmake, "--build", build, "--parallel", "2"], out, "_hdiv_cuda-host", commands)
    if not binary.is_file():
        raise RuntimeError("HDiv host build did not produce the exact Python ABI artifact")
    (out / "hdiv-sdk-private.json").write_text(json.dumps({str(p): sha(p) for p in (ng, netgen)}, indent=2) + "\n")
    archive = out / ("radia_hdiv_kernels.lib" if platform.system() == "Windows" else "libradia_hdiv_kernels.a")
    if not archive.is_file():
        raise RuntimeError("HDiv host build did not produce the shared numerical archive")
    with archive.open("rb") as stream:
        if stream.read(8) != b"!<arch>\n" or archive.stat().st_size <= 8:
            raise RuntimeError("Invalid HDiv numerical archive")
    # Build artifacts only: never add an archive or CUDA runtime to the wheel.
    return {"cuda_object": {"name": cuda_object.name, "sha256": sha(cuda_object)},
            "static_archive": {"name": archive.name, "sha256": sha(archive)}}



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cuda-root", type=Path, required=True)
    parser.add_argument("--host-compiler", type=Path, required=True)
    parser.add_argument("--architectures", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--matlab-root", type=Path)
    parser.add_argument("--components", nargs="+", choices=[row[0] for row in COMPONENTS],
                        default=[row[0] for row in COMPONENTS], help="Explicit subset to rebuild from this exact commit")
    args = parser.parse_args()
    commit = identity()
    record = admission(args.admission, args.admission_sha256, commit)
    windows = platform.system() == "Windows"
    if not windows and platform.system() != "Linux":
        raise RuntimeError("Unsupported build OS")
    cuda = args.cuda_root.resolve()
    if windows and (args.matlab_root is None or cuda != (args.matlab_root.resolve() / "sys/cuda/win64/cuda")):
        raise RuntimeError("Windows must use the explicitly named MATLAB-bundled CUDA")
    nvcc, nvcc_version, host_version = check_toolchain(cuda, args.host_compiler, args.architectures, windows)
    cuda_sdk_flags, cuda_sdk_identity = (windows_cuda_layout(cuda, args.matlab_root, nvcc_version)
                                         if windows else ([], None))
    include = headers()
    out = args.output.resolve()
    if out == ROOT or ROOT in out.parents:
        raise RuntimeError("Build output must be outside the source checkout")
    out.mkdir(parents=True, exist_ok=False)
    source_hashes = {str((CORE / name).relative_to(ROOT)).replace("\\", "/"): sha(CORE / name) for name in SOURCES}
    # HDiv sources have one owner and are also consumed by the MEX gateways.
    source_hashes.update(hdiv_source_hashes())
    source_hashes["tools/build_python_gpu.py"] = sha(__file__)
    header_hashes = {p.relative_to(include).as_posix(): sha(p) for p in (include / "pybind11").rglob("*") if p.is_file()}
    commands = []
    for module, cpp, cu in COMPONENTS:
        if module not in args.components:
            continue
        binary = out / (module + sysconfig.get_config_var("EXT_SUFFIX"))
        obj = out / (module + (".obj" if windows else ".o"))
        cuda_command = [nvcc, "-c", "-std=c++17", "-O2", "-ccbin", args.host_compiler,
                        "-I" + str(CORE), "-o", obj, CORE / cu]
        cuda_command += cuda_sdk_flags
        cuda_command += ["-gencode=arch=compute_" + arch + ",code=sm_" + arch for arch in args.architectures]
        if windows:
            cuda_command += ["-Xcompiler", "/MD"]
        if module == "_hdiv_cuda":
            cuda_command += ["--fmad=false", "-I" + str(ROOT / "src/core")]
        if not windows:
            cuda_command += ["-Xcompiler", "-fPIC"]
        execute(cuda_command, out, module + "-cuda", commands)
        if module == "_hdiv_cuda":
            hdiv_artifacts = build_hdiv_host(args, binary, obj, include, cuda, out, commands)
        elif windows:
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
        if module != "_hdiv_cuda":
            execute(command, out, module + "-host", commands)
        if identity() != commit or any(sha(ROOT / name) != digest for name, digest in source_hashes.items()):
            raise RuntimeError("Source changed during build")
        if windows and any(not Path(item["path"]).is_file() or sha(item["path"]) != item["sha256"]
                           for item in cuda_sdk_identity["files"].values()):
            raise RuntimeError("Bundled CUDA SDK inputs changed during build")
        private_receipt = Path(str(binary) + ".private-build.json")
        private_receipt.write_text(json.dumps(dict(source_commit=commit, binary_sha256=sha(binary),
            admission=record, commands=commands[:], source_root=str(ROOT), cuda_root=str(cuda),
            host_compiler_path=str(args.host_compiler), pybind_include=str(include),
            ngsolve_sdk_sha256=sha(out / "hdiv-sdk-private.json") if module == "_hdiv_cuda" else None,
            hdiv_kernel_artifacts=hdiv_artifacts if module == "_hdiv_cuda" else None,
            cuda_sdk_identity=cuda_sdk_identity), indent=2) + "\n", encoding="utf-8")
        receipt = dict(schema="radia.native-build-provenance.v1", source_commit=commit, source_dirty=False,
                       binary_name=binary.name, binary_sha256=sha(binary), source_hashes=source_hashes,
                       pybind_header_hashes=header_hashes, python=sys.version, platform=platform.platform(),
                       nvcc=nvcc_version, host_compiler=host_version, architectures=args.architectures, ptx_architectures=[],
                       fault_injection=False, device_fault_injection=False, lifecycle_diagnostic=False,
                       private_receipt_sha256=sha(private_receipt), runtime_import_checked=False,
                       generated_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        Path(str(binary) + ".build.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    # Operational originals stay private; wheels select only the minimal receipt.
    (out / "admission-private.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
