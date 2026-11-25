"""Build optional Radia pybind CUDA modules from one shared numerical core.

No packages are installed, no CUDA context is created, and no numerical test
is run. Windows uses MSVC 14.44 and MATLAB-bundled CUDA 12.8; Linux uses an
explicit CUDA 12.8 nvcc and host compiler. The runtime remains pip-owned.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import sys
import sysconfig
import time
import tomllib
import zipfile

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "src/core/gpu"
HDIV = CORE / "hdiv"
SHARED_SOURCES = ("csr_pcg_host.hpp", "csr_pcg_cuda.cuh", "csr_pcg_python.cpp",
           "csr_pcg_python_cuda.cu", "ams_amg_kernels.cuh", "resident_amg_host.hpp",
           "resident_amg_cuda.cu", "resident_amg_python.cpp", "snapshot.hpp",
           "ams_core.cuh", "ams_refresh.cuh")
# Radia-only GPU inputs: not shared with the radia-simulink MEX builds.
SOURCES = SHARED_SOURCES + ("ams_python.cpp", "ams_python_cuda.cu", "peec_fin_host.hpp",
                            "peec_fin_python.cpp", "peec_fin_cuda.cu", "../peec_oriented_box.hpp")
COMPONENTS = (("_gpu_csr", "csr_pcg_python.cpp", "csr_pcg_python_cuda.cu"),
              ("_gpu_resident_amg", "resident_amg_python.cpp", "resident_amg_cuda.cu"),
              ("_gpu_ams", "ams_python.cpp", "ams_python_cuda.cu"),
              ("_hdiv_cuda", "hdiv/hdiv_resident_python.cpp", "hdiv/hdiv_resident_python_cuda.cu"),
              ("_fem_cuda", "fem/fem_python.cpp", "fem/fem_python_cuda.cu"),
              ("_peec_fin_cuda", "peec_fin_python.cpp", "peec_fin_cuda.cu"))


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
    expected = ["80"] if windows else ["100"]
    if list(architectures) != expected:
        raise RuntimeError("Explicit architectures must be Windows sm_80 or Linux sm_100")
    if windows:
        compiler = subprocess.run([str(host)], capture_output=True, text=True, encoding="utf-8", errors="replace")
        host_version = compiler.stdout + compiler.stderr
        if not re.search(r"19\.44\.", host_version):
            raise RuntimeError("MSVC 14.44 compiler environment is required")
    else:
        host_version = capture([host, "--version"])
        if not re.search(r"(?:g\+\+|gcc|GCC).*11\.5", host_version):
            raise RuntimeError("Explicit reviewed GCC 11.5 host compiler required")
    for key in ("NVCC_PREPEND_FLAGS", "NVCC_APPEND_FLAGS", "NVCC_CCBIN", "CL", "_CL_", "LINK", "_LINK_", "CFLAGS", "CXXFLAGS", "CPPFLAGS", "CUDAFLAGS"):
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


def verify_windows_fork():
    from ci_fork_wheels import LOCK, verify_installed

    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    verify_installed(lock)
    return lock


def verify_linux_fork(lock_path):
    """Verify wheel bytes and the installed Linux SDK before using its headers."""
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    packages = {"ngsolve": "ngsolve", "netgen-mesher": "netgen"}
    if lock.get("schema") != 1 or set(lock["wheels"]) != packages.keys():
        raise ValueError("Expected a Linux fork lock for NGSolve and Netgen")
    prefix = Path(sys.prefix).resolve()
    if prefix == Path(sys.base_prefix).resolve():
        raise RuntimeError("Linux fork build requires an isolated virtual environment")
    for package, module in packages.items():
        record = lock["wheels"][package]
        version = record["version"]
        source = record["source_commit"]
        version_source = re.search(r"\+radia\.([0-9a-f]{7,40})(?:\.|$)", version)
        if (not re.fullmatch(r"[0-9a-f]{40}", source) or not version_source
                or not source.startswith(version_source[1])):
            raise ValueError(f"Invalid Linux fork source identity: {package}")
        wheel = lock_path.parent / record["file"]
        if (Path(record["file"]).name != record["file"] or wheel.suffix != ".whl"
                or sha(wheel) != record["sha256"]):
            raise RuntimeError(f"Linux fork wheel hash mismatch: {package}")
        dist = metadata.distribution(package)
        if dist.version != version:
            raise RuntimeError(f"Incorrect installed Linux fork version: {package}")
        spec = importlib.util.find_spec(module)
        origin = Path(spec.origin).resolve() if spec and spec.origin else None
        expected_origin = Path(dist.locate_file(f"{module}/__init__.py")).resolve()
        if origin != expected_origin or not expected_origin.is_relative_to(prefix):
            raise RuntimeError(f"Incorrect installed Linux fork import origin: {module}")
        native_count = 0
        with zipfile.ZipFile(wheel) as archive:
            for member in archive.infolist():
                path = PurePosixPath(member.filename)
                if member.is_dir() or path.parts[0].endswith(".dist-info"):
                    continue
                if path.is_absolute() or ".." in path.parts or chr(92) in member.filename:
                    raise ValueError(f"Invalid Linux fork wheel member: {member.filename}")
                if path.parts[0].endswith(".data"):
                    if len(path.parts) < 3 or path.parts[1] not in {"purelib", "platlib"}:
                        # Wheel scripts and headers outside site-packages are not SDK inputs here.
                        if len(path.parts) > 1 and path.parts[1] == "scripts":
                            continue
                        raise ValueError(f"Unsupported Linux fork wheel layout: {member.filename}")
                    path = PurePosixPath(*path.parts[2:])
                installed = Path(dist.locate_file(str(path))).resolve()
                if not installed.is_relative_to(prefix) or not installed.is_file():
                    raise RuntimeError(f"Missing installed Linux fork member: {path}")
                with archive.open(member) as stream:
                    expected = hashlib.file_digest(stream, "sha256").hexdigest()
                if sha(installed) != expected:
                    raise RuntimeError(f"Installed Linux fork member hash mismatch: {path}")
                native_count += bool(re.search(r"\.so(?:\.|$)", path.name))
        if not native_count:
            raise RuntimeError(f"No native libraries in Linux fork wheel: {package}")
    return lock


def headers(lock=None):
    # The installed NGSolve/Netgen ABI owns pybind headers.
    import netgen
    import ngsolve
    requirements = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    for package in ("ngsolve", "netgen-mesher"):
        expected = next(r.split("==", 1)[1] for r in requirements if r.startswith(package + "=="))
        actual = metadata.version(package)
        if lock is not None:
            actual = actual.split("+", 1)[0]
        if actual != expected:
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
    parser.add_argument("--linux-fork-lock", type=Path, help="Hash-pinned Linux fork wheel manifest")
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
    if windows and args.linux_fork_lock:
        parser.error("--linux-fork-lock is Linux-only")
    lock = verify_windows_fork() if windows else (
        verify_linux_fork(args.linux_fork_lock) if args.linux_fork_lock else None
    )
    include = headers(lock)
    out = args.output.resolve()
    if out == ROOT or ROOT in out.parents:
        raise RuntimeError("Build output must be outside the source checkout")
    out.mkdir(parents=True, exist_ok=False)
    source_hashes = {(CORE / name).resolve().relative_to(ROOT.resolve()).as_posix(): sha(CORE / name) for name in SOURCES}
    # HDiv sources have one owner and are also consumed by the MEX gateways.
    source_hashes.update(hdiv_source_hashes())
    source_hashes.update({p.relative_to(ROOT).as_posix(): sha(p) for p in (CORE / "fem").rglob("*") if p.is_file()})
    source_hashes["src/core/gpu/ams_amg_kernels.cuh"] = sha(CORE / "ams_amg_kernels.cuh")
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
            cuda_command.append("-gencode=arch=compute_80,code=compute_80")
        if windows:
            cuda_command += ["-Xcompiler", "/MD"]
        if module == "_hdiv_cuda":
            cuda_command += ["--fmad=false", "-I" + str(ROOT / "src/core")]
        if module == "_fem_cuda":
            cuda_command += ["--fmad=false"]
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
            if module in ("_gpu_resident_amg", "_gpu_ams", "_fem_cuda"):
                command += [cuda / "lib/x64/cublas.lib"]
            if module == "_gpu_ams":
                command += [cuda / "lib/x64/cusolver.lib"]
            if module == "_peec_fin_cuda":
                command += [cuda / "lib/x64/cublas.lib", cuda / "lib/x64/cusolver.lib"]
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
            if module in ("_gpu_resident_amg", "_gpu_ams", "_fem_cuda"):
                command += [libraries["libcublas.so.12"]["path"]]
            if module == "_gpu_ams":
                command += [libraries["libcusolver.so.11"]["path"]]
            if module == "_peec_fin_cuda":
                command += [libraries["libcublas.so.12"]["path"], libraries["libcusolver.so.11"]["path"]]
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
                       nvcc=nvcc_version, host_compiler=host_version, architectures=args.architectures, ptx_architectures=["80"] if windows else [],
                       fault_injection=False, device_fault_injection=False, lifecycle_diagnostic=False,
                       private_receipt_sha256=sha(private_receipt), runtime_import_checked=False,
                       generated_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        Path(str(binary) + ".build.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    # Operational originals stay private; wheels select only the minimal receipt.
    (out / "admission-private.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
