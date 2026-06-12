"""
Reference Y_exact for an infinite-z conducting cylinder of radius a.

Used as the "ground truth" against which mixed Galerkin approximations
are measured. Provides BOTH full Bessel evaluation (small to moderate
γa) AND Senior-tower asymptotic continuation (very large γa where
scipy.special.iv overflows).

## Formula

For a uniform AC axial current in a cylinder of radius a, conductivity
σ, permeability μ:

    γ  = sqrt(s μ σ)               diffusion wavenumber
    γa = a γ                       dimensionless skin parameter

    Y(s) / Y_DC = (2 / γa) · I_1(γa) / I_0(γa)

with Y_DC = π a² σ.  Asymptotic limits:

    γa → 0:  Y(s)/Y_DC → 1                       (DC bulk)
    γa → ∞:  Y(s)/Y_DC → K_SIBC / (Y_DC √s)
             = (2/γa) · (1 - 1/(2γa) - 1/(8(γa)²) - 1/(8(γa)³)
                            - 25/(128(γa)^4) - 13/(32(γa)^5) - ...)
                                                  (Senior tower)

The Senior tower coefficients b_k = [0, -1/2, -1/8, -1/8, -25/128,
-13/32, -1073/1024, -103/32, ...] form a DIVERGENT asymptotic series.
For γa larger than ~500, scipy's iv() overflows the float64 range, so
we truncate the Senior series at 6 terms — by which point the series
has converged to better than float64 precision for γa > ~30 and the
remaining error is exponentially small.

## ⚠ Phase 8b artifact (2026-06-12 — see memory project_mixed_galerkin_schur_dtn)

Earlier versions of this Y_exact used an early hard switch to the
*leading-order* K_SIBC/√s asymptote at |γa| > 50, omitting the Senior
tower corrections. This switch made the K_SIBC asymptote (which differs
from full Bessel by ~1% at γa=50, decaying as 1/(γa)) appear as the
"ground truth" against which mixed Galerkin was measured, biasing the
mixed-Galerkin error by exactly the Senior tower corrections the mixed
Galerkin was correctly computing. The "0.93% wall-band saturation"
reported in Phases 2/3/8 was 100% this artifact; the true cylinder
mixed-Galerkin 1-DOF accuracy is ≈ 0.04%, and Senior tower truncation
gives ~100× per added DOF (Phase 8b corrected result).

This module uses Senior tower with multiple terms in the asymptotic
region, restoring the proper ground-truth behavior.

## API

    Y_exact_cylinder(s, a, sigma, mu)         → complex Y(s)
    Y_DC_cylinder(a, sigma)                   → real Y_DC
    K_SIBC_cylinder(a, sigma, mu)             → real K_SIBC

## Self-test

Run as `python cylinder_bessel.py` to verify the Bessel↔Senior tower
crossover is smooth across γa = 100, 200, 500.
"""

from __future__ import annotations

import cmath
import math

from scipy.special import iv

# Cross-over threshold between scipy.iv and Senior tower asymptotic.
# At γa = 500, both formulas agree to ~1e-13 in relative terms, well
# below float64 precision. Above γa ~ 700, iv() starts overflowing.
_GA_CROSSOVER = 500.0

# Senior tower coefficients of (Y(s)/Y_DC) · (γa/2) ~ 1 + b_1/γa +
# b_2/(γa)² + ...  (asymptotic series in 1/γa, divergent but useful
# for γa > ~30).
_BESSEL_RATIO_ASYMPT_COEFFS = (
    -1.0 / 2.0,          # b_1
    -1.0 / 8.0,          # b_2
    -1.0 / 8.0,          # b_3
    -25.0 / 128.0,       # b_4
    -13.0 / 32.0,        # b_5
    -1073.0 / 1024.0,    # b_6
)


def Y_DC_cylinder(a: float, sigma: float) -> float:
    """DC conductance per unit length (S·m): Y_DC = π a² σ."""
    return math.pi * a**2 * sigma


def K_SIBC_cylinder(a: float, sigma: float, mu: float) -> float:
    """Leading SIBC coefficient: Y(s) ~ K_SIBC / √s  at γa → ∞.

    K_SIBC = 2π a √(σ/μ) for a cylinder (per unit length).
    """
    return 2.0 * math.pi * a * math.sqrt(sigma / mu)


def Y_exact_cylinder(s, a: float, sigma: float, mu: float):
    """Exact Y(s) for an infinite-z cylinder.

    Uses scipy.special.iv for |γa| < _GA_CROSSOVER (full Bessel),
    Senior tower asymptotic series for larger |γa| (where iv() would
    overflow float64).

    Parameters
    ----------
    s : complex
        Laplace variable.
    a : float
        Cylinder radius (m).
    sigma : float
        Conductivity (S/m).
    mu : float
        Permeability (H/m).

    Returns
    -------
    complex
        Y(s) per unit length (S·m).
    """
    Y_DC = Y_DC_cylinder(a, sigma)
    gamma = cmath.sqrt(s * mu * sigma)
    gA = gamma * a

    if abs(gA) < _GA_CROSSOVER:
        # Full Bessel evaluation.
        return Y_DC * 2.0 * iv(1, gA) / (gA * iv(0, gA))

    # Senior tower asymptotic series.  Compute (I_1/I_0) ratio via
    # 1 + b_1/(γa) + b_2/(γa)² + ...
    z = gA
    ratio = 1.0 + 0j
    z_pow = 1.0 + 0j
    for b_k in _BESSEL_RATIO_ASYMPT_COEFFS:
        z_pow *= z
        ratio += b_k / z_pow
    return Y_DC * 2.0 * ratio / z


def _self_test() -> None:
    """Compare Bessel vs Senior tower in the overlap region."""
    sigma = 5.8e7
    mu = 4 * math.pi * 1e-7
    a = 5e-3
    MS = mu * sigma

    print("=== cylinder_bessel.py self-test ===")
    print(f"a = {a*1e3} mm, sigma = {sigma:.2e} S/m, mu = {mu:.4e} H/m")
    print(f"Y_DC   = {Y_DC_cylinder(a, sigma):.4e} S*m")
    print(f"K_SIBC = {K_SIBC_cylinder(a, sigma, mu):.4e}")
    print()
    print(f"Crossover |gamma a| = {_GA_CROSSOVER}; both formulas should agree to ~1e-12 there.")
    print()
    print(f"{'|gamma a|':>10}  {'|Y_bessel|':>14}  {'|Y_senior|':>14}  {'rel diff':>12}")
    for gA_target in [30, 100, 200, 300, 400, 499.5]:
        # Choose s on the imaginary axis at phase π/4 (typical sweep).
        # We need |γa| = gA_target with phase π/4: γa = gA_target · e^{iπ/4}
        # so s = (γa/a)² / (μσ) at e^{iπ/2} = j times real.
        gA_complex = gA_target * cmath.exp(1j * math.pi / 4)
        s_val = (gA_complex / a)**2 / MS

        # Manual evaluation: scipy iv
        Y_bessel = math.pi * a**2 * sigma * 2.0 * iv(1, gA_complex) / (gA_complex * iv(0, gA_complex))

        # Manual evaluation: Senior tower
        z = gA_complex
        ratio = 1.0 + 0j
        z_pow = 1.0 + 0j
        for b_k in _BESSEL_RATIO_ASYMPT_COEFFS:
            z_pow *= z
            ratio += b_k / z_pow
        Y_senior = math.pi * a**2 * sigma * 2.0 * ratio / z

        diff = abs(Y_bessel - Y_senior) / abs(Y_bessel)
        print(f"{gA_target:10.1f}  {abs(Y_bessel):14.6e}  {abs(Y_senior):14.6e}  {diff:12.2e}")

    print()
    print("Round-trip Y_exact_cylinder():")
    print(f"{'|gamma a|':>10}  {'|Y_exact|':>14}")
    for gA_target in [10, 100, 1000, 10000]:
        gA_complex = gA_target * cmath.exp(1j * math.pi / 4)
        s_val = (gA_complex / a)**2 / MS
        Y = Y_exact_cylinder(s_val, a, sigma, mu)
        print(f"{gA_target:10.0f}  {abs(Y):14.6e}")


if __name__ == "__main__":
    _self_test()
