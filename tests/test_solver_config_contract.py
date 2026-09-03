"""Runtime contract for the retained global nonlinear solver settings."""

import pytest

import radia


EXPECTED_KEYS = {
    "relax_param",
    "keep_magnetization",
    "newton_method",
    "newton_damping",
    "newton_damping_max_iter",
    "newton_damping_min_omega",
    "b_input_newton",
    "b_input_hantila",
    "hantila_alpha",
    "hantila_relax",
}
RETIRED_GLOBAL_OPTIONS = (
    "hacapk_eps",
    "hacapk_leaf",
    "hacapk_eta",
    "hmatrix_eps",
    "bicgstab_tol",
)


def test_get_solver_config_reports_only_retained_global_state():
    assert set(radia.GetSolverConfig()) == EXPECTED_KEYS


@pytest.mark.parametrize("option", RETIRED_GLOBAL_OPTIONS)
def test_retired_kernel_options_fail_loudly(option):
    with pytest.raises(ValueError, match=f"unknown SolverConfig option: {option}"):
        radia.SolverConfig(**{option: 1.0})
