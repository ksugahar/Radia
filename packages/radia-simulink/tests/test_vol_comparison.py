"""Numerical acceptance rejects invalid solver evidence without MATLAB."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "gpu_vol_comparison",
    Path(__file__).resolve().parents[3] / "validation_test/radia_simulink/compare_vol_gpu.py",
)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_matching_solution_and_finite_residual_pass():
    assert comparison.validate_solution([1., 2.], [1., 2.], 1e-12) == 0


@pytest.mark.parametrize("residual", [np.nan, np.inf, -np.inf, -1e-12, 1e-5, True])
def test_invalid_true_residual_fails(residual):
    with pytest.raises(RuntimeError, match="true residual"):
        comparison.validate_solution([1., 2.], [1., 2.], residual)


@pytest.mark.parametrize("reference,solution", [
    ([], []), ([0.], [0.]), ([1., 2.], [1.]),
    ([np.nan], [1.]), ([1.], [np.inf]), ([1.], [1.1]),
])
def test_invalid_or_inaccurate_vectors_fail(reference, solution):
    with pytest.raises(RuntimeError):
        comparison.validate_solution(reference, solution, 1e-12)
