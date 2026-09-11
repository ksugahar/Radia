"""Missing native symbols must not silently select a NumPy substitute."""
import numpy as np
import pytest
import radia._radia_pybind as native
from radia.vim import _eddy_hybrid as eddy


@pytest.mark.parametrize("name", ["_HybridVIMSolve", "_HybridVIMSchurComplement",
    "_ProjectedBaseMatrix", "_ChargeGramHMatrix", "_ReducedBlockMatrix",
    "_SkinImpedance", "_SIBCAdmittanceTail", "_SIBCSchurTerminationImpedance",
    "_SIBCSchurTerminationAdmittance", "_EVRSTMethodAlgebra"])
def test_missing_native_kernel_raises(monkeypatch, name):
    monkeypatch.delattr(native, name, raising=False)
    with pytest.raises(RuntimeError, match=name):
        eddy._radia_cpp_kernel(name)


@pytest.mark.parametrize("backend", ["auto", "cpp"])
def test_evrs_does_not_fall_back(monkeypatch, backend):
    monkeypatch.delattr(native, "_EVRSTMethodAlgebra", raising=False)
    matrices = [np.eye(1)] * 7
    with pytest.raises(RuntimeError, match="_EVRSTMethodAlgebra"):
        eddy.EVRSTMethodAlgebra(*matrices, backend=backend)


def test_explicit_python_backend_stays_explicit(monkeypatch):
    monkeypatch.delattr(native, "_EVRSTMethodAlgebra", raising=False)
    expected = object()
    monkeypatch.setattr(eddy, "_evrs_tmethod_algebra_numpy", lambda *args: expected)
    assert eddy.EVRSTMethodAlgebra(*([np.eye(1)] * 7), backend="python") is expected


@pytest.mark.parametrize("name,call", [
    ("_SkinImpedance", lambda: eddy.SkinImpedance(1j, 1.0)),
    ("_SIBCAdmittanceTail", lambda: eddy.SIBCAdmittanceTail(1j, 1.0, 1.0)),
    ("_SIBCSchurTerminationImpedance", lambda: eddy.SIBCSchurTerminationImpedance(1j, 1.0)),
    ("_SIBCSchurTerminationAdmittance", lambda: eddy.SIBCSchurTerminationAdmittance(1j, 1.0)),
])
def test_public_sibc_helpers_require_native(monkeypatch, name, call):
    monkeypatch.delattr(native, name, raising=False)
    with pytest.raises(RuntimeError, match=name):
        call()
