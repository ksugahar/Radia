"""Install or verify the exact public FE SDK, and smoke-test a release wheel."""
import argparse
import hashlib
import importlib
import importlib.metadata as metadata
import json
import math
from pathlib import Path
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
FE_PACKAGES = {"ngsolve": "ngsolve", "netgen-mesher": "netgen"}
CPU_MODULES = ("_radia_pybind", "peec_matrices", "sparsesolv_ngsolve", "axifem")
GPU_MODULES = ("_gpu_csr", "_gpu_resident_amg", "_gpu_ams", "_hdiv_cuda", "_fem_cuda", "_peec_fin_cuda")


def pins():
    requirements = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    return {name: next(item.split("==", 1)[1] for item in requirements if item.startswith(name + "=="))
            for name in FE_PACKAGES}


def verify():
    prefix = Path(sys.prefix).resolve()
    if prefix == Path(sys.base_prefix).resolve():
        raise RuntimeError("A release FE verification requires an isolated venv")
    result = {}
    for name, expected in pins().items():
        actual = metadata.version(name)
        if actual != expected:
            raise RuntimeError(f"{name} resolved to {actual}; required public version {expected}")
        package = importlib.import_module(FE_PACKAGES[name])
        origin = Path(package.__file__).resolve()
        if not origin.is_relative_to(prefix):
            raise RuntimeError(f"{name} resolved outside the venv: {origin}")
        result[name] = {"version": actual, "origin": str(origin)}
    return result


def smoke(wheel, gpu):
    # Do not preinstall FE packages here. This exercises ordinary pip resolution.
    requirement = str(wheel.resolve()) + ("[gpu]" if gpu else "")
    python = [sys.executable, "-I"]
    subprocess.run([*python, "-m", "pip", "install", "--only-binary=:all:", requirement], check=True)
    subprocess.run([*python, "-m", "pip", "check"], check=True)
    subprocess.run([*python, str(Path(__file__).resolve()), "import", "--wheel", str(wheel.resolve()),
                    *(["--gpu"] if gpu else [])], check=True)


def import_wheel(wheel, gpu):
    result = verify()
    prefix = Path(sys.prefix).resolve()
    for name in CPU_MODULES:
        module = importlib.import_module("radia." + name)
        origin = Path(module.__file__).resolve()
        if not origin.is_relative_to(prefix) or origin.suffix.lower() not in {".pyd", ".so"}:
            raise RuntimeError(f"Unexpected native module origin for {name}: {origin}")
        result[name] = str(origin)
    if gpu:
        from radia.gpu import _native
        for name in GPU_MODULES:
            origin = Path(_native(name).__file__).resolve()
            if not origin.is_relative_to(prefix):
                raise RuntimeError(f"Unexpected native module origin for {name}: {origin}")
            result[name] = str(origin)
    import radia
    magnet = radia.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1])
    field = float(radia.Fld(magnet, "bz", [0, 0, 0.02]))
    if not math.isfinite(field) or field == 0:
        raise RuntimeError(f"Native field smoke calculation failed: {field}")
    with wheel.open("rb") as stream:
        result["wheel_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
    result["smoke_bz"] = field
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("install", "verify", "smoke", "import"))
    parser.add_argument("--wheel", type=Path)
    parser.add_argument("--gpu", action="store_true")
    args = parser.parse_args()
    if args.action in {"smoke", "import"}:
        if args.wheel is None:
            parser.error("--wheel is required for wheel acceptance")
        (smoke if args.action == "smoke" else import_wheel)(args.wheel, args.gpu)
        return
    if args.action == "install":
        if sys.prefix == sys.base_prefix:
            raise RuntimeError("Public SDK install requires an isolated venv")
        subprocess.run([sys.executable, "-I", "-m", "pip", "install", "--only-binary=:all:",
                        *(f"{name}=={version}" for name, version in pins().items())], check=True)
    print(json.dumps(verify()))


if __name__ == "__main__":
    main()
