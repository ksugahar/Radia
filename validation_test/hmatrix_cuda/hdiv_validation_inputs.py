"""Explicit, checked inputs for independent HDiv validation runners."""
from pathlib import Path
import base64
import hashlib
import json
import os


def checked_mesh(mesh, contract, mesh_sha256, contract_sha256, receipt, receipt_sha256,
                 provenance, report, tet_only=True):
    """Verify the producer's strict receipt without installing Cubit on this host."""
    from radia.hdiv_gpu import _verify_mesh_receipt

    mesh_bytes = Path(mesh).read_bytes()
    contract_bytes = Path(contract).read_bytes()
    for name, payload, expected in (("mesh", mesh_bytes, mesh_sha256),
                                     ("label contract", contract_bytes, contract_sha256)):
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError("Checked input hash mismatch: " + name)
    evidence = _verify_mesh_receipt(
        mesh_bytes, contract_bytes, Path(receipt).read_bytes(), receipt_sha256,
        json.loads(Path(provenance).read_text(encoding="utf-8-sig")))
    report_bytes = base64.b64decode(evidence["report_base64"], validate=True)
    result = json.loads(report_bytes)
    quality = result["quality"]
    if tet_only:
        count = quality.get("volume_element_count")
        if type(count) is not int or count <= 0 or quality.get("tetrahedron_count") != count:
            raise ValueError("Checked mesh is not exclusively tetrahedral")
    Path(report).write_bytes(report_bytes)
    return dict(mesh_sha256=mesh_sha256, contract_sha256=contract_sha256,
                receipt_sha256=receipt_sha256, checker=evidence["checker"], report=result)


def dependency_directories(paths):
    if not paths:
        raise ValueError("Explicit dependency directories are required")
    directories = [Path(p).resolve(strict=True) for p in paths]
    if any(not p.is_dir() for p in directories):
        raise ValueError("Dependency path must be a directory")
    os.environ["PATH"] = os.pathsep.join(map(str, directories)) + os.pathsep + os.environ["PATH"]
    return directories, [os.add_dll_directory(str(p)) for p in directories] if os.name == "nt" else []


def process_affinity(pid=None):
    if os.name != "nt":
        return None
    import ctypes
    kernel=ctypes.WinDLL("kernel32",use_last_error=True)
    kernel.OpenProcess.restype=ctypes.c_void_p
    kernel.OpenProcess.argtypes=[ctypes.c_uint32,ctypes.c_bool,ctypes.c_uint32]
    handle=kernel.OpenProcess(0x1000,False,os.getpid() if pid is None else pid)
    if not handle:raise ctypes.WinError(ctypes.get_last_error())
    kernel.GetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_size_t),ctypes.POINTER(ctypes.c_size_t)]
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    try:
        process,system=ctypes.c_size_t(),ctypes.c_size_t()
        if not kernel.GetProcessAffinityMask(handle,ctypes.byref(process),ctypes.byref(system)):
            raise ctypes.WinError(ctypes.get_last_error())
        return process.value
    finally:kernel.CloseHandle(handle)


def thread_record(count):
    if count < 1:
        raise ValueError("CPU threads must be positive")
    affinity = None
    if os.name == "nt":
        import ctypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
        process, system = ctypes.c_size_t(), ctypes.c_size_t()
        if not kernel.GetProcessAffinityMask(kernel.GetCurrentProcess(), ctypes.byref(process), ctypes.byref(system)):
            raise ctypes.WinError(ctypes.get_last_error())
        affinity = process.value
        if count > affinity.bit_count():
            raise ValueError("Requested NGSolve threads exceed admitted CPU affinity")
    return dict(cpu_ngsolve=count, cpu_blas=1, matlab_single_comp_thread=True,
                process_affinity_mask=affinity)


def oracle_vectors(fixture):
    import numpy as np
    expected, rhs = (np.asarray(fixture[name]).reshape(-1) for name in ("expected", "rhs"))
    system, means = fixture["system"], fixture["mean_operator"]
    if expected.shape != rhs.shape or system.shape != (rhs.size, rhs.size) or means.shape != (3, rhs.size):
        raise ValueError("Independent oracle vector/system/mean shapes mismatch")
    return expected, system, rhs


def cpu_receipt(path, native_path):
    receipt=json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if receipt.get("source_dirty") is not False or receipt.get("binary_sha256") != hashlib.sha256(Path(native_path).read_bytes()).hexdigest():
        raise ValueError("Loaded CPU native binary does not match clean build receipt")
    return receipt


def single_thread_blas():
    from threadpoolctl import threadpool_limits
    limit=threadpool_limits(limits=1,user_api="blas")
    return limit, observed_blas()


def observed_blas():
    from threadpoolctl import threadpool_info
    observed=threadpool_info()
    if any(p["num_threads"] != 1 for p in observed if p["user_api"] == "blas"):
        raise RuntimeError("BLAS single-thread policy did not take effect")
    return observed


def sorted_oracle_csr(value):
    """Explicit producer ordering only; never merge duplicate coordinates."""
    matrix=value.tocsr(copy=True)
    matrix.sort_indices()
    if not matrix.nnz or not matrix.has_canonical_format:
        raise ValueError("Oracle CSR is empty or contains duplicate coordinates; no coalescing allowed")
    return matrix
