"""Explicit sparse CPU/GPU operations; importing this module does not load CUDA.

NGSolve owns assembly and CPU hierarchy construction. GPU execution consumes
an owned algebraic snapshot, never a MATLAB adapter or a CPU fallback.
"""
from __future__ import annotations

import math
from numbers import Integral, Real

__all__ = ["solve_spd", "AMG", "AMS", "ResidentAMG"]


def _controls(tolerance, max_iterations):
    if (isinstance(tolerance, bool) or not isinstance(tolerance, Real)
            or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6):
        raise ValueError("tolerance must be finite in (0, 1e-6]")
    if (isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral)
            or not 1 <= max_iterations <= 100000):
        raise ValueError("max_iterations must be an integer in [1, 100000]")
    return float(tolerance), int(max_iterations)


def _vector(value, n):
    import numpy as np
    data = np.asarray(value)
    if (data.dtype != np.dtype("float64") or data.shape != (n,)
            or not np.isfinite(data).all()):
        raise ValueError("a finite float64 vector of the system dimension is required")
    return np.array(data, order="C", copy=True)


def _csr(matrix):
    import numpy as np
    if getattr(matrix, "is_complex", False) or matrix.height != matrix.width:
        raise ValueError("a real square matrix is required")
    n = matrix.height
    values, columns, offsets = (np.asarray(v) for v in matrix.CSR())
    limit = np.iinfo(np.int32).max
    if not 0 < n <= limit or values.size > limit:
        raise OverflowError("CSR dimensions exceed native int32 capacity")
    if values.ndim != 1 or values.dtype != np.dtype("float64") or not np.isfinite(values).all():
        raise ValueError("CSR values must be a finite float64 vector")
    for indices in (columns, offsets):
        if (indices.ndim != 1 or indices.dtype.kind not in "iu"
                or (indices.size and (indices.min() < 0 or indices.max() > limit))):
            raise ValueError("CSR indices must be nonnegative int32-compatible vectors")
    if (offsets.size != n + 1 or columns.size != values.size
            or offsets[0] != 0 or offsets[-1] != values.size
            or (np.diff(offsets.astype(np.int64)) < 0).any()
            or (columns >= n).any()):
        raise ValueError("invalid CSR layout")
    return dict(rows=int(n), cols=int(n), rowptr=np.array(offsets, dtype=np.int32),
                columns=np.array(columns, dtype=np.int32), values=np.array(values, copy=True))


def _residual(system, x, rhs):
    import numpy as np
    offsets, columns, values = (system[k] for k in ("rowptr", "columns", "values"))
    residual = np.empty(rhs.size, dtype=np.float64)
    for i in range(rhs.size):
        start, end = offsets[i:i + 2]
        residual[i] = values[start:end] @ x[columns[start:end]] - rhs[i]
    scale = np.linalg.norm(rhs)
    return float(np.linalg.norm(residual) / scale) if scale else float(np.linalg.norm(residual))


def _matrix_residual(matrix, template, x):
    import numpy as np
    native = template.CreateVector()
    native.FV().NumPy()[:] = x
    residual = template.CreateVector()
    residual.data = matrix * native - template
    scale = np.linalg.norm(template.FV().NumPy())
    norm = np.linalg.norm(residual.FV().NumPy())
    return float(norm / scale) if scale else float(norm)


def _checked_result(system, rhs, result, backend, tolerance, *, solve=True,
                    original_residual=None):
    x, info = result
    x = _vector(x, rhs.size)
    if not isinstance(info, dict) or info.get("backend") != backend:
        raise RuntimeError("native operation did not report the requested backend")
    info = dict(info)
    if solve:
        relative = original_residual(x) if original_residual else _residual(system, x, rhs)
        if not math.isfinite(relative) or relative > min(1e-6, max(10 * tolerance, 1e-12)):
            raise RuntimeError("true original-system residual exceeds acceptance")
        info.update(true_relative_residual=relative, residual_checked=True)
    else:
        info.update(operation="preconditioner_apply", residual_checked=False)
    return x, info


def _native(name):
    from importlib import import_module
    import hashlib
    import json
    from pathlib import Path
    import sysconfig
    import re
    root = Path(__file__).resolve().parent
    binary = root / (name + sysconfig.get_config_var("EXT_SUFFIX"))
    receipt_path = Path(str(binary) + ".build.json")
    if not binary.is_file() or not receipt_path.is_file():
        raise ImportError(f"Missing Radia GPU artifact or build provenance: {name}")
    receipt = json.loads(receipt_path.read_text())
    import radia._radia_pybind as core
    core_receipt = json.loads(Path(str(Path(core.__file__)) + ".build.json").read_text())
    for record, path in [(receipt, binary), (core_receipt, Path(core.__file__))]:
        if (record.get("schema") != "radia.native-build-provenance.v1"
                or not re.fullmatch(r"[0-9a-f]{40}", record.get("source_commit", ""))
                or record.get("source_dirty") is not False
                or record.get("binary_name") != path.name
                or record.get("binary_sha256") != hashlib.sha256(path.read_bytes()).hexdigest()):
            raise ImportError("Radia native build provenance mismatch")
    if receipt.get("source_commit") != core_receipt.get("source_commit"):
        raise ImportError("Radia CPU/GPU source revisions differ")
    from ._cuda_runtime import prepare_cuda_runtime
    prepare_cuda_runtime()
    module = import_module("." + name, __package__)
    prepare_cuda_runtime()  # Check the actual post-import DLL/ELF closure too.
    return module



def solve_spd(matrix, rhs, *, backend, tolerance=1e-10, max_iterations=10000,
              preconditioner=None):
    """Solve a real SPD NGSolve matrix with an explicitly selected backend.

    CPU uses NGSolve sparsecholesky under the caller's TaskManager ownership.
    GPU uses Jacobi-PCG or an explicitly prepared AMG/AMS hierarchy. A failed
    GPU request raises. Preconditioners are bound to their snapshot operator.
    """
    if backend not in ("cpu", "gpu"):
        raise ValueError("backend must be explicitly cpu or gpu")
    tolerance, max_iterations = _controls(tolerance, max_iterations)
    if (getattr(matrix, "is_complex", False) or matrix.height != matrix.width
            or matrix.height <= 0):
        raise ValueError("a real square matrix is required")
    system = _csr(matrix) if backend == "gpu" else {"rows": matrix.height}
    b = _vector(rhs.FV().NumPy(), matrix.height)
    if backend == "cpu":
        if preconditioner is not None:
            raise ValueError("the independent CPU route uses NGSolve sparsecholesky")
        inverse = matrix.Inverse(inverse="sparsecholesky")
        x = rhs.CreateVector()
        x.data = inverse * rhs
        result = (x.FV().NumPy(), dict(backend="cpu", solver="NGSolve sparsecholesky"))
    elif preconditioner is not None:
        if not isinstance(preconditioner, _Resident):
            raise TypeError("preconditioner must be a Radia GPU AMG or AMS owner")
        import numpy as np
        if any(not np.array_equal(system[k], preconditioner._system[k])
               for k in ("rows", "cols", "rowptr", "columns", "values")):
            raise ValueError("the matrix differs from the GPU hierarchy snapshot")
        result = preconditioner.solve(b, tolerance=tolerance, max_iterations=max_iterations)
    else:
        import numpy as np
        diagonal = np.empty(b.size, dtype=np.float64)
        for i in range(b.size):
            start, end = system["rowptr"][i:i + 2]
            diagonal[i] = system["values"][start:end][system["columns"][start:end] == i].sum()
        if not np.isfinite(diagonal).all() or (diagonal <= 0).any():
            raise ValueError("GPU Jacobi-PCG requires a finite positive diagonal")
        native = _native("_gpu_csr")
        result = native.solve(system["rowptr"], system["columns"], system["values"],
                              b, np.ascontiguousarray(1. / diagonal), tolerance, max_iterations)
    # NGSolve's operator owns storage semantics (including symmetric storage).
    # Certifying only the exported CSR could accept a truncated triangle.
    return _checked_result(system, b, result, backend, tolerance,
                           original_residual=lambda x: _matrix_residual(matrix, rhs, x))


class _Resident:
    def __init__(self, hierarchy, *, backend, kind):
        if backend != "gpu":
            raise ValueError("AMG/AMS requires explicit backend='gpu'")
        self._owner = None
        snapshot = hierarchy._accelerator_snapshot()
        if snapshot.get("schema") != "radia." + kind + "-accelerator.v1":
            raise ValueError("unsupported " + kind.upper() + " accelerator snapshot")
        system = snapshot["A"]
        # Own the host operator used for independently recomputing residuals.
        import numpy as np
        self._system = {k: np.array(v, copy=True) if hasattr(v, "dtype") else v
                        for k, v in system.items()}
        for value in self._system.values():
            if isinstance(value, np.ndarray):
                value.flags.writeable = False
        self._kind = kind
        self._owner = _native("_gpu_ams").Resident(snapshot)

    def _require_open(self):
        if self._owner is None:
            raise RuntimeError("GPU hierarchy is closed")
        return self._owner

    def solve(self, rhs, *, tolerance=1e-12, max_iterations=10000):
        """Solve the snapshot operator; every matrix/vector operation is on GPU."""
        owner = self._require_open()
        tolerance, max_iterations = _controls(tolerance, max_iterations)
        b = _vector(rhs, self._system["rows"])
        result = owner.solve(b, tolerance, max_iterations)
        return _checked_result(self._system, b, result, "gpu", tolerance)

    def apply(self, rhs):
        """Apply one preconditioner cycle. This does not certify a solved system."""
        owner = self._require_open()
        b = _vector(rhs, self._system["rows"])
        return _checked_result(self._system, b, owner.apply(b), "gpu", 0., solve=False)

    def close(self):
        owner = self._owner
        if owner is not None:
            owner.close()
            self._owner = None

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()


class AMG(_Resident):
    """Upload an already constructed Radia CompactAMG hierarchy to CUDA."""
    def __init__(self, hierarchy, *, backend):
        super().__init__(hierarchy, backend=backend, kind="amg")


class AMS(_Resident):
    """Upload an already constructed real FP64 Radia AMS hierarchy to CUDA.

    Only the existing AMS cycles 1 and 7 are supported. NGSolve still owns
    FE setup and orientation. This interface does not introduce HDiv-MMM.
    """
    def __init__(self, hierarchy, *, backend):
        super().__init__(hierarchy, backend=backend, kind="ams")


class ResidentAMG:
    """Own a checked native GPU hierarchy; construction stays with the caller.

    No CPU solve or hierarchy construction is substituted for GPU execution.
    The original-system residual and lifecycle contracts are owned by the
    shared native implementation.
    """

    def __init__(self, levels, coarse_inverse, *, backend):
        if backend != "gpu":
            raise ValueError('ResidentAMG requires explicit backend="gpu"')
        self._owner = _native("_gpu_resident_amg").Resident(levels, coarse_inverse)

    def solve(self, rhs, *, tolerance=1e-12, max_iterations=10000):
        import math
        from numbers import Integral, Real
        if (isinstance(tolerance, bool) or not isinstance(tolerance, Real)
                or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6):
            raise ValueError("AMG tolerance must be finite in (0, 1e-6]")
        if (isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral)
                or not 1 <= max_iterations <= 100000):
            raise ValueError("AMG iteration limit must be an integer in [1, 100000]")
        if self._owner is None:
            raise RuntimeError("Resident AMG owner is closed")
        return self._owner.solve(rhs, float(tolerance), int(max_iterations))

    def close(self):
        owner = getattr(self, "_owner", None)
        if owner is not None:
            owner.close()
            self._owner = None

    def __enter__(self):
        if self._owner is None:
            raise RuntimeError("Resident AMG owner is closed")
        return self

    def __exit__(self, *_):
        self.close()
