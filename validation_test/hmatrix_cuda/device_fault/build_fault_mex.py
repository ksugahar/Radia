"""Test-only fault MEX linked from the same nvcc object as the shipped plain MEX; never shipped.

Its own binary name and receipt schema keep it apart from `radia_hdiv_device`.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu-build", type=Path, required=True, help="tools/build_python_gpu.py output with _hdiv_cuda")
    parser.add_argument("--matlab-root", type=Path, required=True)
    parser.add_argument("--host-compiler", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary-source-sha", help="Kernel commit when it differs from HEAD; the native closure must then be checked by the caller")
    args = parser.parse_args()
    head = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    if subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip():
        raise RuntimeError("Clean source required")
    gpu = args.gpu_build.resolve()
    private = json.loads(next(gpu.glob("_hdiv_cuda*.pyd.private-build.json")).read_text(encoding="utf-8"))
    obj = gpu / private["hdiv_kernel_artifacts"]["cuda_object"]["name"]
    if private["source_commit"] != (args.binary_source_sha or head) or sha(obj) != private["hdiv_kernel_artifacts"]["cuda_object"]["sha256"]:
        raise RuntimeError("The nvcc object does not belong to the expected commit")
    cuda = args.matlab_root.resolve() / "sys/cuda/win64/cuda"
    if Path(private["cuda_root"]).resolve() != cuda:
        raise RuntimeError("The fault MEX links the same MATLAB-bundled CUDA as the nvcc object")
    out = args.output.resolve()
    if out == ROOT or ROOT in out.parents:
        raise RuntimeError("Build output must be outside the source checkout")
    out.mkdir(parents=True, exist_ok=False)
    cmake = shutil.which("cmake")
    if not cmake:
        raise RuntimeError("CMake is required")
    posix = lambda value: Path(value).as_posix()
    commands = [[cmake, "-S", posix(HERE), "-B", posix(out / "build"), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
                 "-DCMAKE_CXX_COMPILER=" + posix(args.host_compiler), "-DMatlab_ROOT_DIR=" + posix(args.matlab_root),
                 "-DRADIA_CUDA_OBJECT=" + posix(obj), "-DRADIA_OUTPUT=" + posix(out)]
                + ["-DRADIA_%s=%s" % (name.upper(), posix(cuda / "lib/x64" / (name + ".lib")))
                   for name in ("cudart", "cublas", "cusolver")],
                [cmake, "--build", posix(out / "build"), "--parallel", "2"]]
    for index, command in enumerate(commands):
        with (out / ("cmake-%d.log" % index)).open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = out / "radia_hdiv_device_fault.mexw64"
    if not binary.is_file():
        raise RuntimeError("The link did not produce radia_hdiv_device_fault.mexw64")
    receipt = dict(schema="radia.hdiv-fault-mex-build.v1", source_commit=head, kernel_commit=private["source_commit"],
                   fault_injection=True, binary_name=binary.name, binary_sha256=sha(binary),
                   cuda_object=dict(name=obj.name, sha256=sha(obj)),
                   commands=[list(map(str, c)) for c in commands],
                   generated_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    Path(str(binary) + ".build.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
