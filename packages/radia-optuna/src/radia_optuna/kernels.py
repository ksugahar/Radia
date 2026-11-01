"""Private explicit backend loader; no native import at package import time."""
import importlib
from importlib import metadata
import os
import sys
from pathlib import Path

_dll_handles=[]
def _loaded_path(name):
    import ctypes
    from ctypes import wintypes
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.GetModuleHandleW.argtypes=[wintypes.LPCWSTR]
    kernel.GetModuleHandleW.restype=wintypes.HMODULE
    module=kernel.GetModuleHandleW(name)
    if not module:return None
    kernel.GetModuleFileNameW.argtypes=[wintypes.HMODULE,wintypes.LPWSTR,wintypes.DWORD]
    output=ctypes.create_unicode_buffer(32768)
    if not kernel.GetModuleFileNameW(module,output,len(output)):
        raise ImportError('Cannot identify already loaded CUDA runtime')
    return Path(output.value).resolve()
def load(backend):
    if backend not in ('cpu','gpu'):
        raise ValueError('An explicit cpu or gpu backend is required')
    if sys.platform!='win32':
        raise RuntimeError('The qualified Optuna native ABI is Windows-only')
    # nvcc fatbinary registration can map cudart on DLL load, even on CPU.
    # Mapping the approved runtime is distinct from creating a CUDA context.
    # MATLAB-only users need not install this optional Python primitive binding.
    if not _dll_handles:
        dist=metadata.distribution('nvidia-cuda-runtime-cu12')
        if dist.version!='12.8.57':
            raise ImportError('Optuna native binding requires CUDA runtime12.8.57')
        files=[Path(dist.locate_file(f)).resolve() for f in dist.files or ()
               if str(f).lower().endswith('.dll')]
        if not files or not all(f.is_file() for f in files):raise ImportError('CUDA runtime DLLs are missing')
        for file in files:
            loaded=_loaded_path(file.name)
            if loaded is not None and loaded!=file:
                raise ImportError('Different CUDA runtime already loaded: '+file.name)
        pending=[]
        try:
            for directory in sorted({f.parent for f in files}):
                pending.append(os.add_dll_directory(str(directory)))
        except BaseException:
            for handle in reversed(pending):handle.close()
            raise
        _dll_handles.extend(pending)
    return importlib.import_module('radia_optuna._optuna_kernels')
