"""Run-local Windows DLL resolution for MATLAB AMS and FMM validation."""

import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[2]


def _require_dll(directory, name):
    path = directory / name
    if not path.is_file():
        raise FileNotFoundError(f"Missing MATLAB runtime dependency {name}: {path}")
    return path


def matlab_dll_directories():
    prefix = Path(sys.prefix).resolve()
    base = Path(sys.base_prefix).resolve()
    if prefix == base:
        raise RuntimeError("MATLAB validation requires the isolated fork venv")
    public = os.environ.get("RADIA_MATLAB_PUBLIC_FE") == "1"
    if public:
        runpy.run_path(str(ROOT / "tools/verify_public_fe_sdk.py"))["verify"]()
    else:
        verifier = runpy.run_path(str(ROOT / "tools/ci_fork_wheels.py"))
        lock = json.loads(verifier["LOCK"].read_text(encoding="utf-8"))
        verifier["verify_installed"](lock)
    net = Path(importlib.util.find_spec("netgen").origin).resolve().parent
    mkl = prefix / "Library" / "bin"
    occt = prefix / "bin"
    for distribution, name, directory in (
        ("mkl", "mkl_rt.3.dll", mkl),
        ("netgen-occt", "TKOffset.dll", occt),
    ):
        path = _require_dll(directory, name)
        package = metadata.distribution(distribution)
        if not any(Path(package.locate_file(item)).resolve() == path.resolve()
                   for item in package.files or () if item.name.lower() == name.lower()):
            raise RuntimeError(f"MATLAB dependency {name} is outside its installed package record")
    _require_dll(base, f"python{sys.version_info.major}{sys.version_info.minor}.dll")
    openblas = net.parent / "ngsolve_openblas"
    return (net, mkl, occt, base, openblas) if public and openblas.is_dir() else (net, mkl, occt, base)


def start_matlab(engine, options):
    directories = matlab_dll_directories()
    print("MATLAB child DLL directories: " + json.dumps(list(map(str, directories))), flush=True)
    original = os.environ.get("PATH")
    os.environ["PATH"] = os.pathsep.join([*map(str, directories), original or ""])
    try:
        return engine.start_matlab(options)
    finally:
        if original is None:
            os.environ.pop("PATH", None)
        else:
            os.environ["PATH"] = original


def probe_ams_mex(engine, path):
    """Exercise the loader in the new MATLAB process before any numerical call."""
    path = Path(path).resolve()
    _require_dll(path.parent, path.name)
    expected = ("[x,info]=radia_cuda_ams(" if path.stem.endswith("_validation")
                else "A resident HCurl command is required")
    engine.addpath(str(path.parent), "-begin", nargout=0)
    try:
        if Path(engine.which(path.stem)).resolve() != path:
            raise RuntimeError(f"MATLAB resolved a different MEX for {path}")
        try:
            engine.feval(path.stem, nargout=0)
        except Exception as error:
            if expected not in str(error):
                raise RuntimeError(f"MATLAB could not load {path.name}: {error}") from error
        else:
            raise RuntimeError(f"{path.name} did not return its expected no-input diagnostic")
    finally:
        engine.rmpath(str(path.parent), nargout=0)
