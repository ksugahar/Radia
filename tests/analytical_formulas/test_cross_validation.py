"""Cross-validation of analytical_formulas against the rest of Radia.

These tests bridge the closed-form reference layer with the existing
analytical implementations elsewhere in the package, so that a sign
flip or factor mismatch cannot pass silently in either.
"""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    rect_magnet_2d_B,
    shielding_factor_sphere,
)
from radia.analytical_formulas.rect_magnet_2d import MU_0


# ---------------------------------------------------------------------------
# 2D rectangular bar  vs.  3D cuboid (long-z limit)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aspect_z", [10.0, 30.0, 100.0])
def test_rect_2d_matches_long_3d_cuboid(aspect_z):
    """The 2D bar B is the L_z -> infty limit of a 3D cuboid at z = 0.

    Both routines accept M in A/m and return B in Tesla, so a direct
    component-wise comparison is meaningful.
    """
    from radia.analytical_magnet import CuboidMagnet

    a, b = 0.02, 0.01
    Lz = aspect_z * max(a, b)
    Mx = 1.0e6

    # 3D cuboid covers [-a, +a] x [-b, +b] x [-Lz, +Lz].
    cuboid = CuboidMagnet(
        center=[0.0, 0.0, 0.0],
        dimensions=[2 * a, 2 * b, 2 * Lz],
        magnetization=[Mx, 0.0, 0.0],
    )

    # Sample external points at z = 0 (mid-plane) along a line at y = b/2.
    xs = np.array([0.05, 0.10, 0.15])
    for x in xs:
        B3d = np.array(cuboid.get_B([x, b * 0.5, 0.0]))
        B2d = rect_magnet_2d_B([x, b * 0.5], a, b, M=(Mx, 0.0))
        # B_z must be ~zero for both; B_x, B_y of the 3D field at z = 0
        # converges to the 2D result as Lz / a -> infty.
        assert abs(B3d[2]) < 1e-10
        np.testing.assert_allclose(
            B3d[:2], B2d, rtol=4.0 / aspect_z,
            atol=1e-3 * MU_0 * Mx,
        )


def test_rect_2d_far_field_matches_2d_dipole_law():
    """At large r, |B| of any 2D bar -> mu_0 m / (2 pi r**2),
    where the 2D moment per unit length is m = M_x * cross-section.
    """
    a, b = 0.005, 0.003
    Mx = 1.0e6
    m_2D = Mx * 4.0 * a * b
    rs = np.array([0.5, 1.0, 5.0, 20.0])
    for r in rs:
        # On the +x axis, the M-along-x bar's far-field B is along +x.
        B = rect_magnet_2d_B([r, 0.0], a, b, M=(Mx, 0.0))
        Bd = MU_0 * m_2D / (2.0 * math.pi * r * r)
        # Convergence rate goes like (a/r)**2, so tolerance scales likewise.
        rel_tol = 5.0 * (max(a, b) / r) ** 2
        assert B[0] == pytest.approx(Bd, rel=rel_tol)


# ---------------------------------------------------------------------------
# Spherical shielding via the high-mu_r asymptote
# ---------------------------------------------------------------------------


def test_sphere_shielding_high_mu_limit():
    """Closed-form high-mu_r limit. From the exact eq 24

        S_sph = 9 mu_r / [(mu_r + 2)(2 mu_r + 1) - 2 (a/b)**3 (mu_r - 1)**2]

    the leading order at mu_r >> 1 with fixed wall thickness is

        S_sph ~ S_asym / (1 + S_asym),
        S_asym = 9 / [2 mu_r (1 - (a/b)**3)],

    which for a thin wall ``t = b - a << b`` reduces to
    ``S_asym ~ 3 / (2 mu_r t / b)``, the textbook Thomson asymptote.
    """
    for mu_r in (1e4, 1e6, 1e8):
        for t_over_b in (1e-3, 1e-2, 0.05):
            a = 1.0 - t_over_b
            S = shielding_factor_sphere(a, 1.0, mu_r)
            S_asym = 9.0 / (2.0 * mu_r * (1.0 - a ** 3))
            S_pred = S_asym / (1.0 + S_asym)
            # Convergence rate is 1 / mu_r at fixed t/b.
            assert S == pytest.approx(S_pred, rel=5.0 / mu_r)


# ---------------------------------------------------------------------------
# Plate eddy: total dipole moment against simple geometric estimate
# ---------------------------------------------------------------------------


def test_plate_eddy_total_moment_against_quick_estimate():
    """The total in-plane current density integrated as a magnetic
    moment ``m_z = (1/2) integral (r x J)_z dx dy`` should be
    proportional to ``sigma * Bdot * (a*b)**2`` by dimensional
    analysis. The PDF gives the leading coefficient through

        m_z = sigma Bdot int_{-a}^{+a} int_{-b}^{+b} T_z dx dy

    after integration by parts (with T_z = 0 on the boundary).
    Here we just check the dimensional scaling and the negative sign
    expected from Lenz's law (m_z opposes B_dot > 0).
    """
    from radia.analytical_formulas import plate_eddy_T

    sigma, Bdot = 5.0e7, 1.0
    a, b = 0.10, 0.05
    nx, ny = 81, 81
    x = np.linspace(-a, a, nx)
    y = np.linspace(-b, b, ny)
    X, Y = np.meshgrid(x, y)
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    T = plate_eddy_T(pts, a, b, sigma, Bdot, n_terms=120).reshape(X.shape)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    integrated_T = np.sum(T) * dx * dy
    m_z = integrated_T  # m_z = integral T dx dy (per unit thickness)
    # Sign: for Bdot > 0, induced current opposes increase -> T < 0 mostly,
    # so m_z < 0.
    assert m_z < 0.0
    # Dimensional check: |m_z| / (sigma Bdot a**2 b**2) is O(1) (the dimensionless
    # geometric factor depends weakly on b/a).
    scale = sigma * Bdot * a ** 2 * b ** 2
    assert 0.001 < abs(m_z) / scale < 1.0
