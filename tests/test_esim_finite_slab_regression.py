"""Surrogate test of the adopted deep-skin ESIM slab.

Guards: changes to depth-mesh grading and conservative conductance assembly
that alter real surface impedance or the diffusion phase; both conductance
x1.001 and sign-reversal mutations were detected by the phase sensor.
Assumes: surrogate test of adoption record
induction_heating/20261023_esim_slab_impedance (convergence level 1), and
surrogate test of adoption record
induction_heating/20261023_esim_slab_phase (convergence level 1).
N25 on the fixed smooth grading, constant real mu and sigma, exp(+j omega t),
a/delta about 105. The baseline depends on depth_mesh grading and the numerical
library route; changed grading requires reviewed regeneration, not tolerance
relaxation. This checks historical outputs, not analytic accuracy. It does not
certify finite-thickness effects, adaptive refinement, nonlinear BH, full coils,
GPU or FEM. Runtime is metadata, never an assertion.
"""
from pathlib import Path

from radia.esim_cell_problem import ESIMFiniteSlabSolver
from tools.regression_lock import assert_regression, load_baseline


def test_fixed_slab_impedance_and_phase_match_adopted_history():
    conditions = {
        "geometry": "slab", "half_depth_m": 0.005,
        "sigma_S_per_m": 5e6, "frequency_Hz": 222e3,
        "mu_r": 100, "H_peak_A_per_m": 1000, "initial_nodes": 25,
        "mesh": "fixed physical grading; no adaptive outer refinement",
        "production_entry_point": "radia.esim_cell_problem.ESIMFiniteSlabSolver._solve_on_mesh",
    }
    solver = ESIMFiniteSlabSolver(
        conditions["half_depth_m"], sigma=conditions["sigma_S_per_m"],
        frequency=conditions["frequency_Hz"], mu_r=conditions["mu_r"],
        geometry=conditions["geometry"], n_nodes=conditions["initial_nodes"],
    )
    result = solver._solve_on_mesh(
        conditions["H_peak_A_per_m"], tol=1e-6, max_iter=200, relaxation=0.5
    )
    phase = float(result["H_solution"][1].imag / conditions["H_peak_A_per_m"])
    assert phase < 0, (
        f"ESIM diffusion phase: actual Im(H1)/H0={phase}; expected <0 for "
        "exp(+j omega t). This surrogate depends on depth_mesh mesher grading."
    )
    baseline = load_baseline(
        Path(__file__).parent / "fixtures" / "esim_slab_A25_cf77cd6e0.json"
    )
    try:
        assert_regression({"Re_Z_ohm": float(result["Z"].real),
                           "Im_H1_over_H0": phase}, baseline, conditions=conditions)
    except AssertionError as exc:
        raise AssertionError(
            "ESIM A25 surrogate depends on depth_mesh mesher grading and the "
            f"numerical library route; review baseline changes explicitly. {exc}"
        ) from exc
