"""Lazy Windows CUDA runtime search paths; no MATLAB or compiler discovery."""
from importlib import metadata
import os
import sys
import threading

_VERSIONS = {
    "nvidia-cuda-runtime-cu12": "12.8.57",
    "nvidia-cublas-cu12": "12.8.3.14",
    "nvidia-cusolver-cu12": "11.7.2.55",
    "nvidia-cusparse-cu12": "12.5.7.53",
    "nvidia-nvjitlink-cu12": "12.8.61",
}
_lock = threading.Lock()
_handles = []


def _loaded_path(name):
    # Kernel loader inspection only: never call a CUDA API during admission.
    import ctypes
    from ctypes import wintypes
    from pathlib import Path
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetModuleHandleW.argtypes = [wintypes.LPCWSTR]
    kernel.GetModuleHandleW.restype = wintypes.HMODULE
    module = kernel.GetModuleHandleW(name)
    if not module:
        return None
    kernel.GetModuleFileNameW.argtypes = [wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
    buffer = ctypes.create_unicode_buffer(32768)
    count = kernel.GetModuleFileNameW(module, buffer, len(buffer))
    if not count or count >= len(buffer):
        raise ImportError("Cannot verify the already loaded CUDA runtime")
    return Path(buffer.value).resolve()


def prepare_cuda_runtime():
    """Register qualified dependency DLL directories after explicit GPU selection.

    Keep directory handles for the Python process lifetime. This registers paths;
    it neither creates a CUDA context nor searches an installed MATLAB runtime.
    Missing or different dependencies raise, with no alternate backend.
    """
    if sys.platform != "win32":
        raise RuntimeError("The qualified Python CUDA binding requires Windows")
    with _lock:
        if _handles:
            return
        directories = set()
        for name, expected in _VERSIONS.items():
            try:
                dist = metadata.distribution(name)
            except metadata.PackageNotFoundError as error:
                raise ImportError("Install radia-simulink[gpu] for the CUDA runtime") from error
            if dist.version != expected:
                raise ImportError(f"CUDA runtime requires {name}=={expected}; found {dist.version}")
            files = [dist.locate_file(f).resolve() for f in dist.files or ()
                     if str(f).lower().endswith(".dll")]
            if not files or not all(f.is_file() for f in files):
                raise ImportError(f"CUDA runtime DLL payload missing for {name}")
            for file in files:
                loaded = _loaded_path(file.name)
                if loaded is not None and loaded != file:
                    raise ImportError(f"A different CUDA runtime is already loaded: {file.name}")
            directories.update(f.parent for f in files)
        pending = []
        try:
            for directory in sorted(directories):
                pending.append(os.add_dll_directory(str(directory)))
        except BaseException:
            for handle in reversed(pending):
                handle.close()
            raise
        _handles.extend(pending)
