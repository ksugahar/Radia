"""Fast regression checks for the NGSolve MATLAB MEX Python oracle fixture."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_HELPER = ROOT / "tests" / "matlab" / "ngsolve_mex_python_reference.py"
VALIDATION_CASE = (
    ROOT / "validation_test" / "ngsolve_integration" / "matlab_mex_native_case.json"
)


def _load_reference_helper():
    spec = importlib.util.spec_from_file_location("ngsolve_matlab_oracle", REFERENCE_HELPER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_python_oracle_covers_native_matlab_ngsolve_contract():
    helper = _load_reference_helper()
    fixture = helper.build_reference(ROOT / "tests" / "fixtures" / "beam" / "affine_field_tetra.vol")

    assert fixture["fixture_schema"] == "radia.ngsolve-matlab-mex-parity.v1"
    assert fixture["mesh_dimension"] == 3
    assert fixture["h1_mass_ndof"] == 10
    assert fixture["hcurl_mass_ndof"] == 14
    assert fixture["hdiv_mass_ndof"] == 30
    for name in (
        "h1_mass_matrix",
        "h1_stiffness_matrix",
        "hcurl_mass_matrix",
        "hcurl_curlcurl_matrix",
        "hdiv_mass_matrix",
        "hdiv_divdiv_matrix",
    ):
        matrix = fixture[name]
        assert matrix.shape[0] == matrix.shape[1]
        assert matrix.nnz > 0
        assert np.count_nonzero(matrix.data == 0.0) == 0
        assert fixture[name.removesuffix("_matrix") + "_nnz"] == matrix.nnz
        assert np.all(np.isfinite(matrix.data))
        assert np.allclose(matrix.toarray(), matrix.toarray().T, rtol=1e-13, atol=1e-14)

    assert fixture["coordinates"].shape == (3, 3)
    assert fixture["constant_vector"].shape == (3, 3)
    assert fixture["h1_grid_evaluation"].shape == (3,)
    assert np.linalg.norm(fixture["h1_solver_residual"], ord=np.inf) < 1e-12


def test_validation_case_declares_native_handle_and_large_scale_boundaries():
    case = json.loads(VALIDATION_CASE.read_text(encoding="utf-8"))
    assert case["schema"] == "radia.ngsolve-matlab-native-validation-case.v1"
    assert case["python_oracle"]["boundary"] == "public NGSolve Python API"
    assert case["matlab_boundary"]["backend"] == "radia_mex through radia.ngsolve handles"
    assert case["next_scale"]["matlab_role"].startswith("submit, monitor")
