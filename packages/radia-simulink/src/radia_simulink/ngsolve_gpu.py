"""Explicit access to NGSolve's optional CUDA device; never accept host fallback.

This adapter does not compile kernels or enable NVRTC. Radia's mexcuda-built
kernels remain a separate native route. CPU workflows must not import this
module to establish their independent CPU reference.
"""
from __future__ import annotations


def cuda_device(*, require_float64=True):
    """Return the upstream CUDA device or raise, without selecting a substitute.

    NGSolve's common API can return a host reference device on a GPU machine.
    Check the actual backend before requesting a device. Kernel compilation
    and runtime/toolchain redistribution require separate acceptance.
    """
    import importlib
    gpu = importlib.import_module("ngsolve.gpu")
    backend = str(gpu.backend).lower()
    if backend != "cuda":
        raise RuntimeError(f"NGSolve CUDA backend required; actual backend is {backend!r}")
    device = gpu.GetGPUDevice()
    if require_float64 and not device.has_float64:
        raise RuntimeError("NGSolve CUDA device does not support required FP64")
    return device
