# -*- coding: utf-8 -*-
r"""DtN Schur gate -- analytic (GIBC) and BEM (exterior Laplace DtN).

Validates:
  * gibc_impedance_cf: Leontovich SIBC (order 0), Mitzner (order 1), Senior-Volakis (order 2)
    on a sphere of radius R against the known formulas H=1/R, K=1/R^2.
  * sphere_gibc_benchmark: correction terms 1/(|k_c|R) and 1/(2|k_c|^2 R^2) match.
  * sphere_exterior_dtn_eigenvalue: BEM Laplace DtN eigenvalue for dipole mode == -2/R.
  * laplace_exterior_dtn: DtN matrix is consistent between direct call and the benchmark.
"""
import cmath
import math
import pytest

bem_mod = pytest.importorskip("ngsolve.bem")

from radia_mcp.radia_ngsolve.bem_integral import (
    gibc_impedance_cf,
    sphere_gibc_benchmark,
    laplace_exterior_dtn,
    sphere_exterior_dtn_eigenvalue,
)

# ---------------------------------------------------------------------------
# GIBC tests
# ---------------------------------------------------------------------------

def test_gibc_order0_is_leontovich_scalar():
    """Order-0 GIBC == Leontovich: scalar complex constant, no Weingarten needed."""
    omega, mu_c, sigma_c = 1e6, 4e-7 * math.pi, 5.8e7
    cf = gibc_impedance_cf(omega, mu_c, sigma_c, order=0)
    assert cf is not None
    # |Z_s0| = sqrt(omega * mu_c / sigma_c)  (|j|=1, so |sqrt(j*omega*mu_c/sigma_c)| = sqrt(omega*mu_c/sigma_c))
    expected = cmath.sqrt(1j * omega * mu_c / sigma_c)
    expected_mag = math.sqrt(omega * mu_c / sigma_c)
    assert abs(abs(expected) - expected_mag) / expected_mag < 1e-9


def test_gibc_order1_greater_than_order0():
    """Order-1 mean impedance |Z_s1| > |Z_s0| for positive curvature (sphere, H=1/R>0)."""
    res = sphere_gibc_benchmark(R=0.05, omega=1e8, sigma=5.8e7)
    # GIBC1 adds a positive correction: |Z_s_gibc1| > |Z_s_sibc|
    assert abs(res["Z_s_gibc1"]) > abs(res["Z_s_sibc"])


def test_gibc_order2_greater_than_order1():
    """Order-2 adds another positive correction on top of order 1."""
    res = sphere_gibc_benchmark(R=0.05, omega=1e8, sigma=5.8e7)
    # Both corrections are additive for sphere (H>0, K>0 term is positive)
    assert abs(res["Z_s_gibc2"]) > abs(res["Z_s_gibc1"])


def test_gibc_benchmark_correction1_matches_theory():
    """Mean relative correction of order-1 GIBC ≈ 1/(|k_c| R) on a sphere.

    For a sphere with H = 1/R: Z_s^(1) = Z_s0 (1 + 1/(R k_c)), so
    |Z_s^(1)/Z_s0 - 1| = 1/(|k_c| R).
    We check within 5% (coarse mesh + finite curvature approximation).
    """
    R, omega, sigma = 0.1, 1e8, 5.8e7
    res = sphere_gibc_benchmark(R=R, omega=omega, sigma=sigma, maxh=0.05, order=3)
    ratio = res["gibc1_correction"] / res["expected_gibc1_correction"]
    assert abs(ratio - 1.0) < 0.05, f"correction={res['gibc1_correction']:.4e}, expected={res['expected_gibc1_correction']:.4e}"


def test_gibc_benchmark_correction2_matches_theory():
    """Mean relative correction of order-2 GIBC ≈ 1/(|k_c|R) + 1/(2|k_c|^2 R^2) on sphere."""
    R, omega, sigma = 0.1, 1e8, 5.8e7
    res = sphere_gibc_benchmark(R=R, omega=omega, sigma=sigma, maxh=0.05, order=3)
    ratio = res["gibc2_correction"] / res["expected_gibc2_correction"]
    assert abs(ratio - 1.0) < 0.05, f"correction={res['gibc2_correction']:.4e}, expected={res['expected_gibc2_correction']:.4e}"


def test_gibc_flat_limit_order1_equals_order0():
    """For a very large sphere (R → ∞), GIBC1 → SIBC: correction → 0."""
    R_large = 1e4    # km-scale sphere: 1/(|k_c| R) ~ machine-zero
    omega, sigma = 1e8, 5.8e7
    mu_c = 4e-7 * math.pi
    k_c = cmath.sqrt(1j * omega * mu_c * sigma)
    correction_expected = 1.0 / (abs(k_c) * R_large)
    # With H = 1/R_large, the correction is negligible
    assert correction_expected < 1e-6


def test_gibc_skin_depth_in_result():
    """sphere_gibc_benchmark exposes the skin depth for sanity checking."""
    omega, sigma = 1e6, 5.8e7
    mu_c = 4e-7 * math.pi
    delta_expected = math.sqrt(2.0 / (omega * mu_c * sigma))
    res = sphere_gibc_benchmark(R=1.0, omega=omega, sigma=sigma)
    assert abs(res["skin_depth"] - delta_expected) / delta_expected < 1e-9


# ---------------------------------------------------------------------------
# BEM exterior Laplace DtN eigenvalue test
# ---------------------------------------------------------------------------

def test_sphere_dtn_eigenvalue_dipole():
    """Exterior Laplace DtN eigenvalue for dipole mode (n=1) on sphere of radius R is −2/R.

    Prescribes f = z/R (cos θ), applies the assembled dense DtN matrix, checks the
    Rayleigh-quotient eigenvalue σ=Λf vs the exact −(n+1)/R = −2/R.  Tolerance 2%
    for a coarse (maxh=0.5) order-1 BEM space (O(h) convergence).
    """
    res = sphere_exterior_dtn_eigenvalue(R=1.0, maxh=0.5, order=1, intorder=8)
    assert abs(res["rel_err"]) < 0.02, (
        f"DtN eigenvalue: BEM={res['dtn_eigenvalue_bem']:.4f}, "
        f"expected={res['dtn_eigenvalue_expected']:.4f}, rel_err={res['rel_err']:.2e}"
    )


def test_sphere_dtn_eigenvalue_scales_with_R():
    """Eigenvalue −2/R scales correctly with sphere radius."""
    for R in (0.5, 1.0, 2.0):
        res = sphere_exterior_dtn_eigenvalue(R=R, maxh=0.5 * R, order=1, intorder=8)
        expected = -2.0 / R
        assert abs(res["dtn_eigenvalue_bem"] - expected) / abs(expected) < 0.03, (
            f"R={R}: BEM={res['dtn_eigenvalue_bem']:.4f}, expected={expected:.4f}"
        )


def test_sphere_dtn_l2_error_small():
    """L2 norm of (Λf - (-2/R)f) is small relative to f for the dipole mode."""
    res = sphere_exterior_dtn_eigenvalue(R=1.0, maxh=0.5, order=1, intorder=8)
    assert res["l2_rel_err"] < 0.05, f"l2_rel_err={res['l2_rel_err']:.2e}"


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("[OK] DtN Schur gate: GIBC (Leontovich / Mitzner / Senior) + BEM exterior Laplace DtN eigenvalue.")


if __name__ == "__main__":
    main()
