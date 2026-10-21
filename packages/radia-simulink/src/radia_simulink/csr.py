"""Explicit NGSolve CPU or shared native CUDA solve for real SPD matrices.

The caller owns NGSolve assembly, numbering and any TaskManager region.
GPU support is optional and imported only after explicit selection. No Engine.
"""
from __future__ import annotations


def solve_spd(matrix, rhs, *, backend, tolerance=1e-10, max_iterations=10000):
    """Return an owned FP64 solution array and original-system residual metadata.

    CPU uses NGSolve sparsecholesky, independently of shared CUDA kernels.
    GPU copies the assembled CSR system to the native Jacobi-PCG implementation.
    No implicit fallback or dtype conversion is provided.
    """
    if backend not in ("cpu", "gpu"):
        raise ValueError("backend must be explicitly cpu or gpu")
    import math
    import numpy as np
    if not isinstance(max_iterations, int) or isinstance(max_iterations, bool) or not 1 <= max_iterations <= 100000:
        raise ValueError("max_iterations must be an integer in [1,100000]")
    if isinstance(tolerance, bool) or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6:
        raise ValueError("tolerance must be finite in (0,1e-6]")
    b = np.asarray(rhs.FV().NumPy())
    if b.ndim != 1 or b.dtype != np.dtype("float64") or not np.isfinite(b).all():
        raise ValueError("rhs must be an NGSolve real FP64 vector with finite entries")
    if getattr(matrix, "is_complex", False):
        raise ValueError("matrix must be real FP64")
    if matrix.height != b.size or matrix.width != b.size or b.size == 0:
        raise ValueError("matrix must be square and match rhs")
    if backend == "cpu":
        inverse = matrix.Inverse(inverse="sparsecholesky")
        native = rhs.CreateVector()
        native.data = inverse * rhs
        solution = np.array(native.FV().NumPy(), dtype=np.float64, copy=True)
        info = dict(backend="cpu", solver="NGSolve sparsecholesky", iterations=None)
    else:
        # Import errors propagate; never route a failed GPU request to CPU.
        from ._cuda_runtime import prepare_cuda_runtime
        prepare_cuda_runtime()
        from . import _csr_cuda
        values, columns, row_pointers = matrix.CSR()
        values = np.asarray(values)
        if values.dtype != np.dtype("float64"):
            raise ValueError("matrix CSR values must be real FP64")
        raw_col, raw_row = np.asarray(columns), np.asarray(row_pointers)
        limit = np.iinfo(np.int32).max
        if b.size > limit or values.size > limit:
            raise OverflowError("CSR dimensions exceed native int32 limits")
        for indices in (raw_col, raw_row):
            if indices.dtype.kind not in "iu" or (indices.size and (indices.min() < 0 or indices.max() > limit)):
                raise OverflowError("CSR indices exceed native int32 limits")
        row = np.ascontiguousarray(raw_row, dtype=np.int32)
        col = np.ascontiguousarray(raw_col, dtype=np.int32)
        if row.size != b.size+1 or col.size != values.size or row[0] != 0 or row[-1] != values.size or (np.diff(row) < 0).any() or (col >= b.size).any():
            raise ValueError("invalid CSR layout")
        diagonal = np.empty(b.size, dtype=np.float64)
        for i in range(b.size):
            diagonal[i] = values[row[i]:row[i+1]][col[row[i]:row[i+1]] == i].sum()
        if not np.isfinite(diagonal).all() or (diagonal <= 0).any():
            raise ValueError("GPU Jacobi requires a finite positive diagonal")
        solution, info = _csr_cuda.solve(row, col, np.ascontiguousarray(values),
            np.ascontiguousarray(b), np.ascontiguousarray(1.0/diagonal), tolerance, max_iterations)
        native = rhs.CreateVector()
        native.FV().NumPy()[:] = solution
    residual = rhs.CreateVector()
    residual.data = rhs - matrix * native
    r = np.asarray(residual.FV().NumPy())
    scale = np.linalg.norm(b)
    relative = float(np.linalg.norm(r) / scale) if scale else float(np.linalg.norm(r))
    if not math.isfinite(relative) or relative > min(1e-6, max(10*tolerance, 1e-12)):
        raise RuntimeError("true original-system residual exceeds acceptance")
    return solution, dict(info, true_relative_residual=relative)
