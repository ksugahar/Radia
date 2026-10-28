"""Guards: explicit dispatch, original residuals and owned AMG/AMS lifetimes.
Assumes: native execution is mocked; these are API contracts, not GPU accuracy evidence.
"""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def gpu():
    spec = importlib.util.spec_from_file_location("gpu_entry", ROOT / "src/radia/gpu.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Vector:
    def __init__(self, values):
        self.values = np.asarray(values)
    def FV(self):
        return self
    def NumPy(self):
        return self.values
    def CreateVector(self):
        return Vector(np.zeros_like(self.values))
    @property
    def data(self):
        return self.values
    @data.setter
    def data(self, values):
        self.values[:] = values.values if isinstance(values, Vector) else values
    def __sub__(self, other):
        return Vector(self.values - other.values)


class Matrix:
    height = width = 2
    is_complex = False
    def __init__(self):
        self.values = np.array([4., -1., -1., 4.])
        self.columns = np.array([0, 1, 0, 1], dtype=np.int32)
        self.rows = np.array([0, 2, 4], dtype=np.int32)
    def CSR(self):
        return self.values, self.columns, self.rows
    def __mul__(self, vector):
        return Vector(self.values.reshape(2, 2) @ vector.values)
    def Inverse(self, **kwargs):
        assert kwargs == {"inverse": "sparsecholesky"}
        class Inverse:
            def __mul__(self, vector):
                return np.linalg.solve(np.array([[4., -1.], [-1., 4.]]), vector.values)
        return Inverse()


def test_cpu_does_not_load_cuda(gpu, monkeypatch):
    monkeypatch.setattr(gpu, "_native", lambda _: pytest.fail("CPU loaded CUDA"))
    x, info = gpu.solve_spd(Matrix(), Vector([1., 2.]), backend="cpu", tolerance=1e-12)
    np.testing.assert_allclose(x, [0.4, 0.6], rtol=1e-14)
    assert info["backend"] == "cpu" and info["true_relative_residual"] < 1e-14


def test_gpu_dispatch_owns_contiguous_typed_input(gpu, monkeypatch):
    matrix = Matrix()
    def solve(row, col, val, rhs, diagonal, tolerance, limit):
        assert row.dtype == col.dtype == np.dtype("int32")
        assert val.dtype == rhs.dtype == diagonal.dtype == np.dtype("float64")
        assert all(a.flags.c_contiguous for a in (row, col, val, rhs, diagonal))
        assert tolerance == 1e-12 and limit == 1000
        np.testing.assert_array_equal(diagonal, [0.25, 0.25])
        assert not np.shares_memory(val, matrix.values)
        return np.array([0.4, 0.6]), dict(backend="gpu")
    monkeypatch.setattr(gpu, "_native", lambda name: SimpleNamespace(solve=solve))
    x, info = gpu.solve_spd(matrix, Vector([1., 2.]), backend="gpu",
                            tolerance=1e-12, max_iterations=1000)
    np.testing.assert_allclose(x, [0.4, 0.6])
    assert info["residual_checked"] is True


@pytest.mark.parametrize("values,info", [
    ([0., 0.], {"backend": "gpu", "true_relative_residual": 0.}),
    ([0.4, 0.6], {"backend": "cpu"}),
    ([float("nan"), 0.6], {"backend": "gpu"}),
])
def test_bad_native_return_is_not_accepted(gpu, monkeypatch, values, info):
    monkeypatch.setattr(gpu, "_native", lambda _: SimpleNamespace(
        solve=lambda *_: (np.array(values), info)))
    with pytest.raises((RuntimeError, ValueError)):
        gpu.solve_spd(Matrix(), Vector([1., 2.]), backend="gpu", tolerance=1e-12)


def test_missing_gpu_module_never_calls_cpu_inverse(gpu, monkeypatch):
    matrix = Matrix()
    matrix.Inverse = lambda **_: pytest.fail("GPU fell back to CPU")
    def missing(_):
        raise ImportError("unqualified GPU artifact")
    monkeypatch.setattr(gpu, "_native", missing)
    with pytest.raises(ImportError, match="unqualified"):
        gpu.solve_spd(matrix, Vector([1., 2.]), backend="gpu")


def test_original_operator_detects_incomplete_symmetric_csr(gpu, monkeypatch):
    class Triangle(Matrix):
        def CSR(self):
            return (np.array([4., -1., 4.]), np.array([0, 0, 1], dtype=np.int32),
                    np.array([0, 1, 3], dtype=np.int32))
    monkeypatch.setattr(gpu, "_native", lambda _: SimpleNamespace(solve=lambda *_:
        (np.array([0.25, 0.5625]), dict(backend="gpu", true_relative_residual=0.))))
    with pytest.raises(RuntimeError, match="original-system residual"):
        gpu.solve_spd(Triangle(), Vector([1., 2.]), backend="gpu")
    # The independent CPU route still respects the actual NGSolve operator.
    x, _ = gpu.solve_spd(Triangle(), Vector([1., 2.]), backend="cpu")
    np.testing.assert_allclose(x, [0.4, 0.6])


@pytest.mark.parametrize("change", ["negative offset", "bad column", "nonfinite", "wrong dtype"])
def test_bad_csr_fails_before_loading_cuda(gpu, monkeypatch, change):
    matrix = Matrix()
    if change == "negative offset": matrix.rows[1] = -1
    elif change == "bad column": matrix.columns[1] = 2
    elif change == "nonfinite": matrix.values[0] = np.inf
    else: matrix.values = matrix.values.astype(np.float32)
    monkeypatch.setattr(gpu, "_native", lambda _: pytest.fail("invalid input loaded CUDA"))
    with pytest.raises((ValueError, OverflowError)):
        gpu.solve_spd(matrix, Vector([1., 2.]), backend="gpu")


@pytest.mark.parametrize("kind", ["amg", "ams"])
def test_hierarchy_snapshot_lifetime_apply_and_solve(gpu, monkeypatch, kind):
    system = gpu._csr(Matrix())
    snapshot = dict(schema="radia." + kind + "-accelerator.v1", A=system)
    calls = []
    class Native:
        def apply(self, b):
            calls.append("apply")
            return 2 * b, dict(backend="gpu")
        def solve(self, b, tolerance, limit):
            calls.append("solve")
            return np.array([0.4, 0.6]), dict(backend="gpu")
        def close(self): calls.append("close")
    monkeypatch.setattr(gpu, "_native", lambda name: SimpleNamespace(Resident=lambda data: Native()))
    hierarchy = SimpleNamespace(_accelerator_snapshot=lambda: snapshot)
    cls = gpu.AMG if kind == "amg" else gpu.AMS
    with cls(hierarchy, backend="gpu") as owner:
        system["values"][:] = 99.  # Host snapshot ownership must be independent.
        applied, info = owner.apply(np.array([1., 2.]))
        np.testing.assert_array_equal(applied, [2., 4.])
        assert info["operation"] == "preconditioner_apply" and info["residual_checked"] is False
        x, info = gpu.solve_spd(Matrix(), Vector([1., 2.]), backend="gpu", preconditioner=owner)
        np.testing.assert_allclose(x, [0.4, 0.6])
        assert info["true_relative_residual"] < 1e-14
        changed = Matrix();changed.values[0] = 5.
        with pytest.raises(ValueError, match="differs"):
            gpu.solve_spd(changed, Vector([1., 2.]), backend="gpu", preconditioner=owner)
    owner.close()
    assert calls == ["apply", "solve", "close"]
    with pytest.raises(RuntimeError, match="closed"):
        owner.solve(np.ones(2))


@pytest.mark.parametrize("tolerance,limit", [(True, 1), (np.nan, 1), (1e-12, True), (1e-12, 0)])
def test_bad_controls_fail_before_runtime(gpu, monkeypatch, tolerance, limit):
    monkeypatch.setattr(gpu, "_native", lambda _: pytest.fail("invalid controls loaded CUDA"))
    with pytest.raises(ValueError):
        gpu.solve_spd(Matrix(), Vector([1., 2.]), backend="gpu",
                      tolerance=tolerance, max_iterations=limit)


def test_vectorized_csr_residual_empty_rows_and_duplicate_entries(gpu):
    system = dict(rowptr=np.array([0, 2, 2, 3], dtype=np.int32),
                  columns=np.array([0, 0, 2], dtype=np.int32), values=np.array([1., 2., 4.]))
    x = np.array([2., 9., 3.])
    rhs = np.array([6., 0., 12.])
    assert gpu._residual(system, x, rhs) == 0.
    assert gpu._residual(system, x, np.zeros(3)) == pytest.approx(np.linalg.norm(rhs))
