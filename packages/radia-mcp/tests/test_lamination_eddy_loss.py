r"""Lamination eddy-current loss -- core-loss eddy term (#66) -- test.

Classical P = pi^2 sigma d^2 f^2 B^2 / 6 [W/m^3] (thin limit) with the exact skin-effect reduction
F(xi)=(3/xi)(sinh xi - sin xi)/(cosh xi - cos xi), xi=d/delta. The closed-form F is checked against
the numerically-integrated exact 1-D slab loss (an independent reference); the thin limit recovers
the classical d^2 law and the heavy-skin limit -> 3/xi. Tool-independent (analytic + exact slab)."""
import cmath
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (classical_eddy_loss_density, lamination_eddy_skin_factor,
                                           lamination_eddy_loss_density, MU0)


def _F_numeric(xi, n=40001):
    """Exact slab loss / classical loss for a fixed average flux density (delta = 1 units)."""
    import numpy as np
    trap = getattr(np, "trapezoid", None) or np.trapz
    d, k = xi, (1 + 1j)                                  # delta=1: omega=2, k=(1+j)/delta
    Bs = k * d / (2 * cmath.tanh(k * d / 2))             # surface B so average B = 1
    x = np.linspace(0.0, d, n)
    Jz = -Bs * k * np.sinh(k * (x - d / 2)) / cmath.cosh(k * d / 2)
    return trap(np.abs(Jz) ** 2, x) / (2 * d) / (2.0 ** 2 * d ** 2 / 24.0)


def test_skin_factor_matches_exact_slab():
    np = pytest.importorskip("numpy")
    for xi in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
        F = lamination_eddy_skin_factor(xi, 1.0)        # skin_depth=1 so the first arg IS xi
        assert math.isclose(F, _F_numeric(xi), rel_tol=1e-4)


def test_skin_factor_limits():
    # thin limit -> 1 (classical d^2 law), monotone decreasing, heavy skin -> 3/xi
    assert math.isclose(lamination_eddy_skin_factor(1e-3, 1.0), 1.0, rel_tol=1e-9)
    assert lamination_eddy_skin_factor(0.5, 1.0) > lamination_eddy_skin_factor(2.0, 1.0) \
        > lamination_eddy_skin_factor(8.0, 1.0)
    assert lamination_eddy_skin_factor(0.5, 1.0) < 1.0
    assert math.isclose(lamination_eddy_skin_factor(50.0, 1.0), 3.0 / 50.0, rel_tol=1e-6)
    assert math.isclose(lamination_eddy_skin_factor(8.0, 1.0), 3.0 / 8.0, rel_tol=2e-2)


def test_classical_loss_law():
    sigma, d, f, B = 2.0e6, 0.35e-3, 50.0, 1.5
    P = classical_eddy_loss_density(sigma, d, f, B)
    # P = pi^2 sigma d^2 f^2 B^2 / 6
    assert math.isclose(P, math.pi ** 2 * sigma * d ** 2 * f ** 2 * B ** 2 / 6.0, rel_tol=1e-12)
    # d^2 law: doubling thickness quadruples loss; f^2 law: doubling f quadruples loss
    assert math.isclose(classical_eddy_loss_density(sigma, 2 * d, f, B), 4 * P, rel_tol=1e-12)
    assert math.isclose(classical_eddy_loss_density(sigma, d, 2 * f, B), 4 * P, rel_tol=1e-12)
    # sane magnitude: M19 0.35 mm @ 50 Hz, 1.5 T ~ 0.3 W/kg (rho 7650)
    assert 0.1 < P / 7650.0 < 1.0


def test_skin_corrected_density():
    sigma, d, f, B, mur = 2.0e6, 0.35e-3, 50.0, 1.5, 5000.0
    P_cl = classical_eddy_loss_density(sigma, d, f, B)
    # thin sheet at line frequency (delta still >> d) -> within a few % of classical
    P = lamination_eddy_loss_density(sigma, d, f, B, mu_r=mur)
    assert P <= P_cl and math.isclose(P, P_cl, rel_tol=5e-2)
    # the helper composes classical * F(d/delta) exactly; a thick sheet rolls off below classical
    d2 = 3e-3
    delta = math.sqrt(2.0 / (2 * math.pi * f * MU0 * mur * sigma))
    P2 = lamination_eddy_loss_density(sigma, d2, f, B, mu_r=mur)
    P2_cl = classical_eddy_loss_density(sigma, d2, f, B)
    assert math.isclose(P2, P2_cl * lamination_eddy_skin_factor(d2, delta), rel_tol=1e-12)
    assert P2 < 0.8 * P2_cl                                    # skin effect rolls the loss off


if __name__ == "__main__":
    test_skin_factor_matches_exact_slab()
    test_skin_factor_limits()
    test_classical_loss_law()
    test_skin_corrected_density()
    print("[OK] lamination eddy loss: F(xi)==exact slab, thin->classical pi^2 sigma d^2 f^2 B^2/6, heavy->3/xi.")
