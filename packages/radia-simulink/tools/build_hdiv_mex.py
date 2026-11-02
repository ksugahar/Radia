"""Plain MEX for the HDiv device kernels, linked from the reviewed nvcc object.

The CUDA object is the one ``tools/build_python_gpu.py`` compiled for the
``_hdiv_cuda`` pybind11 module at the same clean commit; it is hash-checked
against that build's private receipt and never recompiled here. CMake and
``matlab_add_mex`` link it with MSVC: no mexcuda, NGSolve or Python link.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "packages/radia-simulink/native/hdiv_device"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu-build", type=Path, required=True, help="tools/build_python_gpu.py output with _hdiv_cuda")
    parser.add_argument("--matlab-root", type=Path, required=True)
    parser.add_argument("--host-compiler", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != "Windows":
        raise RuntimeError("radia-simulink MEX builds are Windows-only")
    commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    if subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip():
        raise RuntimeError("Clean source required")
    gpu = args.gpu_build.resolve()
    private = json.loads(next(gpu.glob("_hdiv_cuda*.pyd.private-build.json")).read_text(encoding="utf-8"))
    obj = gpu / private["hdiv_kernel_artifacts"]["cuda_object"]["name"]
    if private["source_commit"] != commit or sha(obj) != private["hdiv_kernel_artifacts"]["cuda_object"]["sha256"]:
        raise RuntimeError("The nvcc object does not belong to this exact commit")
    cuda = args.matlab_root.resolve() / "sys/cuda/win64/cuda"
    if Path(private["cuda_root"]).resolve() != cuda:
        raise RuntimeError("The MEX links the same MATLAB-bundled CUDA as the nvcc object")
    out = args.output.resolve()
    if out == ROOT or ROOT in out.parents:
        raise RuntimeError("Build output must be outside the source checkout")
    out.mkdir(parents=True, exist_ok=False)
    cmake = shutil.which("cmake")
    if not cmake:
        raise RuntimeError("CMake is required for the plain MEX link")
    posix = lambda value: Path(value).as_posix()
    commands = [[cmake, "-S", posix(PROJECT), "-B", posix(out / "build"), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
                 "-DCMAKE_CXX_COMPILER=" + posix(args.host_compiler), "-DMatlab_ROOT_DIR=" + posix(args.matlab_root),
                 "-DRADIA_CUDA_OBJECT=" + posix(obj), "-DRADIA_OUTPUT=" + posix(out)]
                + ["-DRADIA_%s=%s" % (name.upper(), posix(cuda / "lib/x64" / (name + ".lib")))
                   for name in ("cudart", "cublas", "cusolver")],
                [cmake, "--build", posix(out / "build"), "--parallel", "2"]]
    for index, command in enumerate(commands):
        with (out / ("cmake-%d.log" % index)).open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = out / "radia_hdiv_device.mexw64"
    if not binary.is_file():
        raise RuntimeError("The plain MEX link did not produce radia_hdiv_device.mexw64")
    if subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise RuntimeError("Source changed during build")
    receipt = dict(schema="radia.hdiv-mex-build.v1", source_commit=commit, source_dirty=False,
                   binary_name=binary.name, binary_sha256=sha(binary),
                   cuda_object=dict(name=obj.name, sha256=sha(obj)), python_binary_sha256=private["binary_sha256"],
                   source_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(PROJECT.rglob("*")) if p.is_file()},
                   commands=[list(map(str, c)) for c in commands], mexcuda=False,
                   generated_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    Path(str(binary) + ".build.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
