"""Tests for Part 6 / 8 / 9 extensions to analytical_formulas.

Covers:

* ``plate_eddy_dissipation``                 (Part 6 §3)
* ``shielding_factor_*_thin_ac``             (Part 8 §2)
* ``spherical_shell_internal_field``         (Part 6 §2)
* ``conductor_impedance.*``                  (Part 6 §4-§5)
* ``adaptive_quadrature.*``                  (Part 9 §2)
* ``cuboid_average_field.average_B_in_box``  (Part 6 §7, numerical path)
"""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    adaptive_integrate,
    average_B_in_box,
    cylinder_ac_impedance,
    cylinder_dc_resistance,
    cylinder_internal_inductance,
    patterson_nodes_weights,
    planar_surface_impedance,
    plate_eddy_dissipation,
    shielding_factor_cylinder_thin_ac,
    shielding_factor_sphere_thin_ac,
    skin_depth,
    spherical_shell_internal_field,
)
from radia.analytical_formulas.shielding import (
    MU_0,
    shielding_factor_sphere,
)


# ---------------------------------------------------------------------------
# plate_eddy_dissipation
# ---------------------------------------------------------------------------


def test_plate_dissipation_positive_and_quadratic_in_Bdot():
    a, b, d, sigma = 0.05, 0.05, 1e-3, 5e7
    P1 = plate_eddy_dissipation(a, b, d, sigma, 1.0)
    P2 = plate_eddy_dissipation(a, b, d, sigma, 2.0)
    assert P1 > 0
    assert P2 / P1 == pytest.approx(4.0, rel=1e-12)


def test_plate_dissipation_linear_in_d():
    a, b, sigma = 0.05, 0.05, 5e7
    P1 = plate_eddy_dissipation(a, b, 1e-3, sigma, 1.0)
    P2 = plate_eddy_dissipation(a, b, 2e-3, sigma, 1.0)
    assert P2 / P1 == pytest.approx(2.0, rel=1e-12)


def test_plate_dissipation_symmetric_in_aspect_swap():
    """P(a, b) = P(b, a) by axis-swap symmetry of the plate."""
    P_ab = plate_eddy_dissipation(0.05, 0.10, 1e-3, 5e7, 1.0)
    P_ba = plate_eddy_dissipation(0.10, 0.05, 1e-3, 5e7, 1.0)
    assert P_ab == pytest.approx(P_ba, rel=1e-3)


# ---------------------------------------------------------------------------
# shielding (AC + shell internal)
# ---------------------------------------------------------------------------


def test_sphere_thin_ac_dc_limit():
    """At ``omega = 0`` the AC formula gives ``S = 1`` (no shielding)."""
    S = shielding_factor_sphere_thin_ac(0.05, 1e-3, 5e7, omega=0.0)
    assert S.real == pytest.approx(1.0, abs=1e-12)
    assert S.imag == pytest.approx(0.0, abs=1e-12)


def test_sphere_thin_ac_high_frequency_attenuates():
    """|S| decreases with frequency for a fixed shell."""
    omega1, omega2 = 2 * math.pi * 60, 2 * math.pi * 1e4
    S1 = shielding_factor_sphere_thin_ac(0.05, 1e-3, 5e7, omega1)
    S2 = shielding_factor_sphere_thin_ac(0.05, 1e-3, 5e7, omega2)
    assert abs(S2) < abs(S1)


def test_cylinder_thin_ac_consistent_with_sphere_at_dc():
    """At ``omega = 0`` both sphere and cylinder give ``S = 1``."""
    S_cyl = shielding_factor_cylinder_thin_ac(0.05, 1e-3, 5e7, omega=0.0)
    S_sph = shielding_factor_sphere_thin_ac(0.05, 1e-3, 5e7, omega=0.0)
    assert S_cyl.real == pytest.approx(1.0)
    assert S_sph.real == pytest.approx(1.0)


def test_spherical_shell_internal_matches_static_S():
    """``spherical_shell_internal_field`` returns ``S`` consistent with
    the static :func:`shielding_factor_sphere`.
    """
    a, b, mu_r = 0.05, 0.06, 1000.0
    out = spherical_shell_internal_field(1.0, a, b, mu_r)
    S_ref = shielding_factor_sphere(a, b, mu_r)
    assert out["S"] == pytest.approx(S_ref, rel=1e-2)


def test_spherical_shell_field_argument_validation():
    with pytest.raises(ValueError):
        spherical_shell_internal_field(1.0, 0.0, 0.06, 100.0)
    with pytest.raises(ValueError):
        spherical_shell_internal_field(1.0, 0.06, 0.05, 100.0)
    with pytest.raises(ValueError):
        spherical_shell_internal_field(1.0, 0.05, 0.06, 0.5)


# ---------------------------------------------------------------------------
# conductor_impedance
# ---------------------------------------------------------------------------


def test_skin_depth_dc_infinity():
    assert math.isinf(skin_depth(5e7, 0.0))


def test_skin_depth_textbook_copper_60Hz():
    """Copper at 60 Hz has skin depth around 8.4 mm (textbook value)."""
    d = skin_depth(5.96e7, 2 * math.pi * 60)
    assert 8.0e-3 < d < 9.0e-3


def test_planar_Z_phase_45deg():
    """``arg(Z_s) = pi / 4`` regardless of sigma, omega."""
    Z = planar_surface_impedance(5e7, 2 * math.pi * 1e3)
    phase = math.atan2(Z.imag, Z.real)
    assert phase == pytest.approx(math.pi / 4, rel=1e-10)


def test_cylinder_dc_limit():
    """At ``omega = 0`` ``Z = R_dc`` (real, no imaginary)."""
    Z = cylinder_ac_impedance(1e-3, 5e7, omega=0.0)
    assert Z.imag == 0.0
    assert Z.real == pytest.approx(cylinder_dc_resistance(1e-3, 5e7))


def test_cylinder_low_freq_picks_up_internal_inductance():
    """At low frequency ``Im(Z) = omega L_int = omega mu_0 / (8 pi)``."""
    omega = 2 * math.pi * 100  # 100 Hz, much below skin-effect onset
    Z = cylinder_ac_impedance(1e-3, 5e7, omega)
    L_pred = cylinder_internal_inductance()
    assert Z.imag == pytest.approx(omega * L_pred, rel=2e-2)


def test_cylinder_high_freq_R_increases():
    """Skin effect: R(omega) grows above R_dc at high frequency."""
    R_dc = cylinder_dc_resistance(1e-3, 5e7)
    R_high = cylinder_ac_impedance(1e-3, 5e7, 2 * math.pi * 1e9).real
    assert R_high > 5 * R_dc


# ---------------------------------------------------------------------------
# adaptive_quadrature (Patterson)
# ---------------------------------------------------------------------------


def test_patterson_n0_n3_node_counts():
    for n in (0, 1, 2, 3):
        x, w = patterson_nodes_weights(n)
        assert len(x) == 2 ** (n + 1) - 1


def test_patterson_polynomial_exactness():
    """Each Patterson rule of order n is exact for polynomials of
    degree at most a known bound (verified empirically n=0..3).
    """
    expected_max_p = {0: 0, 1: 4, 2: 10, 3: 22}
    for n, max_p in expected_max_p.items():
        x, w = patterson_nodes_weights(n)
        for p in range(0, max_p + 1, 2):
            approx = float(np.sum(w * x ** p))
            exact = 2.0 / (p + 1)
            assert approx == pytest.approx(exact, abs=1e-12)


def test_patterson_invalid_order():
    with pytest.raises(ValueError):
        patterson_nodes_weights(4)


def test_adaptive_integrate_sin():
    val, info = adaptive_integrate(np.sin, 0.0, math.pi, tol=1e-10)
    assert val == pytest.approx(2.0, abs=1e-10)
    assert info["converged"]


def test_adaptive_integrate_caches_evaluations():
    """All cumulative-rule abscissae should be evaluated only once."""
    counter = [0]

    def f(x):
        counter[0] += 1
        return math.exp(-x ** 2)

    val, info = adaptive_integrate(f, -1.0, 1.0, tol=1e-10)
    # The order-3 rule has 15 nodes; should not exceed that even after
    # multiple refinement rounds.
    assert counter[0] == info["n_evals"]
    assert info["n_evals"] <= 15


# ---------------------------------------------------------------------------
# cuboid_average_field (numerical path)
# ---------------------------------------------------------------------------


def test_average_B_far_field_dipole_limit():
    """Far from a magnetised cube, the volume-averaged B over a small
    target box equals the 3D dipole field at the target centre.
    """
    src_min = np.array([-0.01, -0.01, -0.005])
    src_max = +1 * np.array([0.01, 0.01, 0.005])
    src_size = src_max - src_min
    M = np.array([1e6, 0.0, 0.0])
    m = M[0] * np.prod(src_size)
    d = 1.0
    tgt_min = np.array([d - 1e-4, -1e-4, -1e-4])
    tgt_max = np.array([d + 1e-4, +1e-4, +1e-4])
    B_avg = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max, n_quad=4)
    # On +x axis, M along x: dipole B_x = (mu_0 / 4 pi) * 2 m / d**3
    B_dipole = (MU_0 / (4 * math.pi)) * 2 * m / d ** 3
    assert B_avg[0] == pytest.approx(B_dipole, rel=2e-3)
    assert abs(B_avg[1]) < 1e-12
    assert abs(B_avg[2]) < 1e-12


def test_average_B_y_z_zero_by_symmetry():
    """For ``M = (M_x, 0, 0)`` and target/source symmetric in y & z,
    the volume-average ``B_y`` and ``B_z`` vanish exactly.
    """
    src_min = np.array([-0.01, -0.01, -0.005])
    src_max = np.array([0.01, 0.01, 0.005])
    tgt_min = np.array([0.05, -0.001, -0.001])
    tgt_max = np.array([0.07, 0.001, 0.001])
    M = np.array([1e6, 0.0, 0.0])
    B_avg = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max, n_quad=6)
    assert abs(B_avg[1]) < 1e-15
    assert abs(B_avg[2]) < 1e-15


def test_average_B_linear_in_M():
    src_min = np.array([-0.01, -0.01, -0.005])
    src_max = np.array([0.01, 0.01, 0.005])
    tgt_min = np.array([0.04, -0.001, -0.001])
    tgt_max = np.array([0.06, 0.001, 0.001])
    B1 = average_B_in_box([1e6, 0., 0.], src_min, src_max, tgt_min, tgt_max,
                          n_quad=4)
    B2 = average_B_in_box([3e6, 0., 0.], src_min, src_max, tgt_min, tgt_max,
                          n_quad=4)
    # B_x scales linearly in M; B_y, B_z are zero by symmetry
    # (both at floating-point dust level), so an absolute tolerance
    # equal to the dust suppresses the relative-tolerance failure.
    np.testing.assert_allclose(B2, 3.0 * B1, rtol=1e-12, atol=1e-15)
