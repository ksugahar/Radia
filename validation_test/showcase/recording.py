"""Metadata and isolated scratch locations for future showcase calculations."""
from pathlib import Path
import hashlib
import os
import platform
import sys
import tempfile
from importlib.metadata import version, PackageNotFoundError


def runtime_metadata(driver, helpers=(), threads=None):
    versions = {}
    for package in ("radia", "ngsolve", "numpy", "scipy", "cadquery-ocp", "vtk"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            continue
    # Prefer the loaded module version to installed distribution metadata.
    module_sources = {}
    for name, module in list(sys.modules.items()):
        if name in ("radia", "ngsolve", "numpy", "scipy"):
            if getattr(module, "__version__", None) is not None:
                versions[name] = str(module.__version__)
        if name.startswith("radia_mcp.radia_ngsolve"):
            source = getattr(module, "__file__", None)
            if source and Path(source).is_file():
                module_sources[name] = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    root = Path(driver).resolve().parents[2]
    return dict(python=platform.python_version(), versions=versions,
                threads=threads, module_source_sha256=module_sources, driver_sha256=hashlib.sha256(Path(driver).read_bytes()).hexdigest(),
                source_sha256={Path(p).resolve().relative_to(root).as_posix():
                               hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in helpers})


def scratch_directory(requested=None):
    if requested is not None:
        path = Path(requested).resolve()
        root = Path(__file__).resolve().parents[2]
        if path == root or root in path.parents:
            raise ValueError("Scratch directory must be outside the repository")
        path.mkdir(parents=True, exist_ok=True)
        return path
    base = Path("C:/temp") if os.name == "nt" else Path(tempfile.gettempdir())
    base.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="radia-showcase-", dir=base))
