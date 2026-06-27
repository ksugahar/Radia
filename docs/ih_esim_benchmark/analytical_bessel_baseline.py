"""Linear-mu cylinder analytical baseline for ESIM cell verification.

For a constant-mu conductor (cylinder of radius R, conductivity sigma,
permeability mu_0 * mu_r), Maxwell's equations admit a closed-form
surface impedance solution involving modified Bessel functions:

    Z_s(linear) = rho * gamma * I_1(gamma R) / I_0(gamma R)

where rho = 1/sigma and gamma = (1 + j) / delta with skin depth
delta = sqrt(2 / (omega mu_0 mu_r sigma)).

This script:

1. Computes the analytical Z_s for a fixed (R, sigma, mu_r) at a sweep
   of frequencies.
2. Runs the radia ESIMFiniteSlabSolver in cylinder mode at the same
   inputs (no BH curve — pure linear-mu mode) and compares Z_s.
3. Reports relative error on |Z_s|, Re(Z_s), Im(Z_s).

A passing baseline is < 0.1 % relative error on |Z_s| across the
frequency sweep, which proves the radia cell solver's discretisation
(uniform mesh, cell-centred fluxes, L'Hopital BC at r=0) is correctly
implemented.

This is a regression-quality check, not a research demo.  For
nonlinear ESIM validation, use benchmark.py.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import iv as bessel_iv

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "radia"))
from esim_cell_problem import ESIMFiniteSlabSolver  # noqa: E402

MU_0 = 4e-7 * math.pi


def analytical_Zs_cylinder(R, sigma, mu_r, frequency):
    """Closed-form linear-mu cylinder surface impedance.

    Z_s = rho * gamma * I_1(gamma R) / I_0(gamma R)
    where rho = 1/sigma, gamma = (1 + j) / delta,
          delta = sqrt(2 / (omega mu_0 mu_r sigma)).
    """
    omega = 2 * math.pi * frequency
    rho = 1.0 / sigma
    delta = math.sqrt(2.0 / (omega * MU_0 * mu_r * sigma))
    gamma = complex(1.0, 1.0) / delta
    arg = gamma * R
    z = rho * gamma * bessel_iv(1, arg) / bessel_iv(0, arg)
    return complex(z), delta


def radia_Zs_cylinder(R, sigma, mu_r, frequency, n_nodes=2000):
    """Radia ESIMFiniteSlabSolver Z_s, cylinder mode, linear mu.

    Bypasses BH curve by passing mu_r (constant relative permeability).
    """
    esim = ESIMFiniteSlabSolver(
        half_thickness=R,
        bh_curve=None,
        sigma=sigma,
        frequency=frequency,
        mu_r=mu_r,
        n_nodes=n_nodes,
        geometry='cylinder',
    )
    res = esim.solve(H0=1.0)
    return complex(res['Z'])


def sweep():
    R = 5.0e-3          # workpiece radius [m]
    sigma = 2.0e6       # [S/m]
    mu_r = 100.0        # constant permeability
    freqs = [1.0e3, 1.0e4, 5.0e4, 1.0e5, 5.0e5, 1.0e6]

    print(f"Linear-mu cylinder cross-check")
    print(f"  R = {R*1e3:.1f} mm, sigma = {sigma:.2e} S/m, mu_r = {mu_r}")
    print()
    print(f"{'freq[Hz]':>10s}  {'delta[mm]':>9s}  {'xi':>6s}  "
          f"{'Re_anal':>11s}  {'Re_radia':>11s}  {'relErr':>9s}    "
          f"{'Im_anal':>11s}  {'Im_radia':>11s}  {'relErr':>9s}")
    print("-" * 130)

    max_rel = 0.0
    for f in freqs:
        Z_anal, delta = analytical_Zs_cylinder(R, sigma, mu_r, f)
        Z_radia = radia_Zs_cylinder(R, sigma, mu_r, f)
        xi = R / delta
        rel_re = abs(Z_radia.real - Z_anal.real) / max(abs(Z_anal.real), 1e-30)
        rel_im = abs(Z_radia.imag - Z_anal.imag) / max(abs(Z_anal.imag), 1e-30)
        rel_abs = abs(abs(Z_radia) - abs(Z_anal)) / max(abs(Z_anal), 1e-30)
        max_rel = max(max_rel, rel_abs)
        print(f"{f:>10.2e}  {delta*1e3:>9.4f}  {xi:>6.2f}  "
              f"{Z_anal.real:>11.4e}  {Z_radia.real:>11.4e}  {rel_re:>9.2e}    "
              f"{Z_anal.imag:>11.4e}  {Z_radia.imag:>11.4e}  {rel_im:>9.2e}")
    print()
    print(f"Max |Z_s| relative error across sweep: {max_rel:.3e}")
    print(f"Pass criterion (< 5e-3): {'PASS' if max_rel < 5e-3 else 'FAIL'}")
    return max_rel


if __name__ == "__main__":
    sweep()
