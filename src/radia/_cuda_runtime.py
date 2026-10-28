"""Lazy, exact-pin CUDA runtime admission for explicit Python GPU calls."""
from importlib import metadata
import base64
import hashlib
import os
from pathlib import Path
import re
import sys
import threading

# Versions live only in project.optional-dependencies.gpu / installed METADATA.
_LIBRARIES = (
    ("nvidia-nvjitlink-cu12", ("libnvJitLink.so.12",)),
    ("nvidia-cuda-runtime-cu12", ("libcudart.so.12",)),
    ("nvidia-cublas-cu12", ("libcublasLt.so.12", "libcublas.so.12")),
    ("nvidia-cusparse-cu12", ("libcusparse.so.12",)),
    ("nvidia-cusolver-cu12", ("libcusolver.so.11",)),
)
_lock = threading.Lock()
_handles = []
_identity = None


def _pins():
    requirements = metadata.requires("radia") or ()
    result = {}
    for package, _ in _LIBRARIES:
        matches = [re.match(re.escape(package) + r"==([0-9.]+)(?:\s*;|$)", line)
                   for line in requirements]
        pins = {m[1] for m in matches if m}
        if len(pins) != 1:
            raise ImportError(f"radia GPU metadata requires one exact pin for {package}")
        result[package] = pins.pop()
    return result


def _maps():
    paths = set()
    for line in Path("/proc/self/maps").read_text().splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) == 6 and parts[5].startswith("/"):
            paths.add(Path(parts[5]).resolve())
    return paths


def _library_paths(pins, library_specs=None):
    found = {}
    for package, libraries in (_LIBRARIES if library_specs is None else library_specs):
        dist = metadata.distribution(package)
        if dist.version != pins[package]:
            raise ImportError(f"{package}: required {pins[package]}, found {dist.version}")
        for name in libraries:
            files = [f for f in dist.files or () if Path(str(f)).name == name]
            if len(files) != 1:
                raise ImportError(f"Missing or ambiguous pip CUDA library: {name}")
            entry = files[0]
            path = Path(dist.locate_file(entry)).resolve()
            expected = entry.hash
            if not path.is_file() or expected is None or expected.mode != "sha256":
                raise ImportError(f"CUDA library has no SHA256 RECORD: {name}")
            with path.open("rb") as stream:
                actual = hashlib.file_digest(stream, "sha256").digest()
            encoded = base64.urlsafe_b64encode(actual).decode().rstrip("=")
            if encoded != expected.value:
                raise ImportError(f"CUDA library RECORD mismatch: {name}")
            found[name] = {"path": str(path), "sha256": actual.hex(),
                           "distribution": package, "version": dist.version}
    return found


def _check_loaded(libraries, *, require_all=False):
    paths = _maps()
    for name, item in libraries.items():
        expected = Path(item["path"])
        prefix = name.split(".so.")[0] + ".so"
        conflicts = [p for p in paths if p.name.startswith(prefix) and p != expected]
        if conflicts:
            raise ImportError(f"Conflicting CUDA library already loaded: {conflicts[0]}")
        if require_all and expected not in paths:
            raise ImportError(f"Pip CUDA library load could not be verified: {expected}")


def prepare_cuda_runtime():
    """Load only hash-checked pip ELF libraries after explicit GPU selection.

    The site driver remains external. No MATLAB, toolkit search, install or CPU
    fallback is performed. Loading is not a GPU numerical acceptance verdict.
    """
    if sys.platform == "win32":
        from ._cuda_runtime_windows import prepare_cuda_runtime as prepare_windows
        return prepare_windows()
    if sys.platform != "linux":
        raise RuntimeError("Python CUDA requires Linux or Windows")
    global _identity
    with _lock:
        if _identity is not None:
            current = _library_paths(_pins())
            if current != _identity["libraries"]:
                raise ImportError("Pip CUDA runtime identity changed after loading")
            _check_loaded(current, require_all=True)
            return {"pins": dict(_identity["pins"]),
                    "libraries": {k: dict(v) for k, v in _identity["libraries"].items()}}
        pins = _pins()
        libraries = _library_paths(pins)
        _check_loaded(libraries)
        import ctypes
        # Keep successful loads even if a later load fails; re-entry still
        # verifies the actual mapped bytes and never switches providers.
        for item in libraries.values():
            _handles.append(ctypes.CDLL(item["path"], mode=os.RTLD_GLOBAL | os.RTLD_NOW))
        _check_loaded(libraries, require_all=True)
        _identity = {"pins": pins, "libraries": libraries}
        return {"pins": dict(pins), "libraries": {k: dict(v) for k, v in libraries.items()}}
