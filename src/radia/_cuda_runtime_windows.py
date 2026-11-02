"""Lazy Windows pip CUDA directories and actual loaded-DLL identity checks."""
from importlib import metadata
import os
from pathlib import Path
import re
import sys
import threading

_lock = threading.Lock()
_handles = []
_identity = None


def _loaded_modules():
    # Loader inventory only; never call a CUDA API or create a context.
    import ctypes
    from ctypes import wintypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.K32EnumProcessModules.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.HMODULE), wintypes.DWORD, ctypes.POINTER(wintypes.DWORD)]
    kernel.K32EnumProcessModules.restype = wintypes.BOOL
    kernel.GetModuleFileNameW.argtypes = [wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
    modules = (wintypes.HMODULE * 4096)()
    needed = wintypes.DWORD()
    if not kernel.K32EnumProcessModules(kernel.GetCurrentProcess(), modules, ctypes.sizeof(modules), ctypes.byref(needed)) or needed.value > ctypes.sizeof(modules):
        raise ImportError("Cannot inventory the actual Windows CUDA dependency closure")
    paths = set()
    for index in range(needed.value // ctypes.sizeof(wintypes.HMODULE)):
        buffer = ctypes.create_unicode_buffer(32768)
        count = kernel.GetModuleFileNameW(modules[index], buffer, len(buffer))
        if not count or count >= len(buffer):
            raise ImportError("Cannot verify an actual loaded Windows DLL")
        paths.add(Path(buffer.value).resolve())
    return paths


def _check_loaded(libraries):
    family = lambda name: re.sub(r"[0-9]+", "", name.casefold())
    expected = {family(name): Path(item["path"]) for name, item in libraries.items()}
    for path in _loaded_modules():
        selected = expected.get(family(path.name))
        if selected is not None and path != selected:
            raise ImportError(f"Conflicting CUDA DLL already loaded: {path.name}")


def _libraries(pins):
    from ._cuda_runtime import _library_paths
    specs = []
    for package in pins:
        try:
            dist = metadata.distribution(package)
        except metadata.PackageNotFoundError as error:
            raise ImportError("Install radia[gpu] for the exact CUDA runtime") from error
        names = tuple(Path(str(f)).name for f in dist.files or () if str(f).casefold().endswith(".dll"))
        if not names:
            raise ImportError(f"Pip CUDA DLL payload missing: {package}")
        specs.append((package, names))
    return _library_paths(pins, specs)


def prepare_cuda_runtime():
    """Verify pip RECORD bytes and actual DLL paths on cold and warm calls.

    This registers directories; loading remains the native module's task. A
    second call after native import verifies its actual Windows dependencies.
    """
    if sys.platform != "win32":
        raise RuntimeError("Windows DLL loader used on a different OS")
    from ._cuda_runtime import _pins
    global _identity
    with _lock:
        pins = _pins()
        libraries = _libraries(pins)
        if _identity is not None and libraries != _identity:
            raise ImportError("Pip CUDA runtime identity changed after directory registration")
        _check_loaded(libraries)
        if not _handles:
            pending = []
            try:
                for directory in sorted({str(Path(v["path"]).parent) for v in libraries.values()}):
                    pending.append(os.add_dll_directory(directory))
            except BaseException:
                for handle in reversed(pending):
                    handle.close()
                raise
            _handles.extend(pending)
        _identity = libraries
        return {"pins": pins, "libraries": libraries}
