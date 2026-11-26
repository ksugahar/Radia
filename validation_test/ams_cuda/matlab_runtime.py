"""Run-local Windows DLL resolution for MATLAB AMS and FMM validation."""

from concurrent.futures import TimeoutError as FutureTimeout
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import runpy
import sys
from threading import Timer

import psutil

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
        raise RuntimeError("MATLAB validation requires the isolated FE SDK venv")
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


START_TIMEOUT = 120
PATH_TIMEOUT = 30
MEX_TIMEOUT = 60


def _matlab_children():
    children = set()
    for process in psutil.Process().children(recursive=True):
        try:
            if process.name().casefold() == "matlab.exe":
                children.add(process)
        except psutil.NoSuchProcess:
            continue
    return children


def _call(phase, seconds, invoke, owned):
    print(f"BEGIN {phase}", flush=True)

    def stop():
        print(f"TIMEOUT {phase} after {seconds} seconds", flush=True)
        for process in owned():
            try:
                process.kill()
            except psutil.NoSuchProcess:
                continue

    def watchdog():
        try:
            stop()
        finally:
            # The native Engine client itself can be stuck even after MATLAB exits.
            os._exit(124)

    timer = Timer(seconds + 5, watchdog)
    timer.daemon = True
    timer.start()
    future = None
    timeout_type = getattr(sys.modules.get("matlab.engine"), "TimeoutError", FutureTimeout)
    try:
        future = invoke()
        value = future.result(timeout=seconds)
    except (FutureTimeout, timeout_type) as exc:
        stop()
        if phase == "MATLAB engine start" and future is not None:
            future.cancel()
        raise TimeoutError(f"MATLAB phase timed out: {phase}") from exc
    except Exception:
        print(f"ERROR {phase}", flush=True)
        raise
    finally:
        timer.cancel()
    print(f"END {phase}", flush=True)
    return value


def start_matlab(engine, options):
    directories = matlab_dll_directories()
    existing = _matlab_children()
    owned = lambda: _matlab_children() - existing
    session = _call("MATLAB engine start", START_TIMEOUT,
                    lambda: engine.start_matlab(options, background=True), owned)
    try:
        original = _call("MATLAB PATH read", PATH_TIMEOUT,
                         lambda: session.getenv("PATH", background=True), owned)
        selected = os.pathsep.join([*map(str, directories), original])
        _call("MATLAB PATH set", PATH_TIMEOUT,
              lambda: session.setenv("PATH", selected, nargout=0, background=True), owned)
        observed = _call("MATLAB PATH readback", PATH_TIMEOUT,
                         lambda: session.getenv("PATH", background=True), owned)
        if observed != selected:
            raise RuntimeError("MATLAB process did not retain the requested DLL search path")
        print("MATLAB DLL search identity: " + json.dumps({
            "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
            "selected_sha256": hashlib.sha256(selected.encode()).hexdigest(),
            "directories": list(map(str, directories)),
        }), flush=True)
        return session
    except Exception:
        for process in owned():
            try:
                process.kill()
            except psutil.NoSuchProcess:
                continue
        raise


def probe_ams_mex(engine, path):
    """Exercise the loader in the new MATLAB process before any numerical call."""
    path = Path(path).resolve()
    _require_dll(path.parent, path.name)
    expected = ("[x,info]=radia_cuda_ams(" if path.stem.endswith("_validation")
                else "A resident HCurl command is required")
    candidates = _matlab_children()
    pid = _call(f"{path.name} MATLAB process identity", PATH_TIMEOUT,
                lambda: engine.feature("getpid", background=True),
                lambda: candidates if len(candidates) == 1 else ())
    owned = lambda: {process for process in candidates if process.pid == int(pid)}
    if len(owned()) != 1:
        raise RuntimeError(f"Could not establish ownership of MATLAB process {pid}")
    print(f"BEGIN AMS MEX probe {path.name}", flush=True)
    _call(f"{path.name} addpath", PATH_TIMEOUT,
          lambda: engine.addpath(str(path.parent), "-begin", nargout=0, background=True), owned)
    resolved = _call(f"{path.name} origin", PATH_TIMEOUT,
                     lambda: engine.which(path.stem, background=True), owned)
    if Path(resolved).resolve() != path:
        raise RuntimeError(f"MATLAB resolved a different MEX for {path}")
    try:
        _call(f"{path.name} feval", MEX_TIMEOUT,
              lambda: engine.feval(path.stem, nargout=0, background=True), owned)
    except TimeoutError:
        raise
    except Exception as error:
        if expected not in str(error):
            raise RuntimeError(f"MATLAB could not load {path.name}: {error}") from error
    else:
        raise RuntimeError(f"{path.name} did not return its expected no-input diagnostic")
    _call(f"{path.name} rmpath", PATH_TIMEOUT,
          lambda: engine.rmpath(str(path.parent), nargout=0, background=True), owned)
    print(f"END AMS MEX probe {path.name}", flush=True)
