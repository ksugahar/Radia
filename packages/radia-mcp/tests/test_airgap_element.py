r"""Air-gap harmonic element (AGE) -- analytic gate (step 1).

Locks the AGE core before the NGSolve assembly: the annular harmonic transfer
reproduces the ring boundary data, and the closed-form harmonic torque equals the
direct Maxwell-stress quadrature at EVERY radius in the gap -- i.e. the transmitted
torque is RADIUS-/mesh-INDEPENDENT (the AGE gem).  Pure numpy/math, no FE solve.
"""
import cmath
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.airgap_element import (airgap_harmonic_torque,
                                                    annular_field,
                                                    annular_harmonic_coeffs,
                                                    annular_radial_trace,
                                                    maxwell_torque_numeric,
                                                    planar_dtn_matrix,
                                                    planar_field,
                                                    planar_harmonic_coeffs,
                                                    planar_harmonic_thrust,
                                                    planar_radial_trace,
                                                    maxwell_thrust_numeric)

RI, RO = 0.10, 0.11           # rotor / stator ring radii [m] (a 10 mm-bore, 1 mm gap)


def test_transfer_reproduces_ring_boundary_data():
    for n in (1, 2, 3, 6):
        A_in, A_out = 0.7 + 0.2j, -0.3 + 0.5j
        # A_n(ri) == Re[A_in], A_n(ro) == Re[A_out] at theta=0
        assert math.isclose(annular_field(n, RI, RO, A_in, A_out, RI, 0.0),
                            A_in.real, rel_tol=1e-9, abs_tol=1e-12)
        assert math.isclose(annular_field(n, RI, RO, A_in, A_out, RO, 0.0),
                            A_out.real, rel_tol=1e-9, abs_tol=1e-12)


def test_harmonic_torque_is_radius_independent():
    """The gem: closed-form torque == direct Maxwell-stress quadrature at every r."""
    n = 4
    A_in, A_out = 1.0 + 0.0j, 0.6 * cmath.exp(1j * 0.7)   # 0.7 rad rotor/stator phase
    T_formula = airgap_harmonic_torque(n, RI, RO, A_in, A_out, axial_length=0.05)
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        r = RI + frac * (RO - RI)
        T_num = maxwell_torque_numeric(n, RI, RO, A_in, A_out, r, axial_length=0.05)
        assert math.isclose(T_formula, T_num, rel_tol=1e-6, abs_tol=1e-12), (frac, T_formula, T_num)


def test_torque_zero_in_phase_max_at_quadrature():
    n = 2
    # rotor & stator harmonics IN PHASE (real ratio) -> no torque
    assert abs(airgap_harmonic_torque(n, RI, RO, 1.0 + 0j, 0.5 + 0j)) < 1e-12
    # 90 deg phase -> the extremum of T ~ Im(alpha conj(beta))
    T90 = abs(airgap_harmonic_torque(n, RI, RO, 1.0 + 0j, 0.5j))
    T45 = abs(airgap_harmonic_torque(n, RI, RO, 1.0 + 0j, 0.5 * cmath.exp(1j * math.pi / 4)))
    assert T90 > T45 > 0.0


def test_torque_antisymmetric_in_phase_sign():
    n = 3
    Tp = airgap_harmonic_torque(n, RI, RO, 1.0 + 0j, 0.5 * cmath.exp(+1j * 0.6))
    Tm = airgap_harmonic_torque(n, RI, RO, 1.0 + 0j, 0.5 * cmath.exp(-1j * 0.6))
    assert math.isclose(Tp, -Tm, rel_tol=1e-9)            # reversing the lead reverses the torque


def test_radial_trace_matches_finite_difference():
    n, h = 2, 1e-7
    A_in, A_out = 0.8 + 0.1j, 0.4 - 0.3j
    d_ri, d_ro = annular_radial_trace(n, RI, RO, A_in, A_out)
    alpha, beta = annular_harmonic_coeffs(n, RI, RO, A_in, A_out)
    fd_ri = ((alpha * (RI + h) ** n + beta * (RI + h) ** -n)
             - (alpha * (RI - h) ** n + beta * (RI - h) ** -n)) / (2 * h)
    assert math.isclose(d_ri.real, fd_ri.real, rel_tol=1e-5)
    assert math.isclose(d_ri.imag, fd_ri.imag, rel_tol=1e-5)


# ---------------------------------------------------------------------------
# Planar gap DtN tests
# ---------------------------------------------------------------------------

G = 0.001           # gap width [m] (1 mm)
K_WAVE = 100.0      # wavenumber [1/m] (pole pitch ~63 mm)
X_PERIOD = 2.0 * math.pi / K_WAVE


def test_planar_field_reproduces_boundary_data():
    """A(0) == A_bottom and A(g) == A_top (Dirichlet interpolation)."""
    for k in (50.0, K_WAVE, 300.0):
        for A_b, A_t in ((1.0 + 0j, 0.0 + 0j), (0.0 + 0j, 1.0 + 0j),
                          (0.7 + 0.2j, -0.3 + 0.5j)):
            assert math.isclose(planar_field(k, G, A_b, A_t, 0.0),
                                 A_b.real, abs_tol=1e-12)
            assert math.isclose(planar_field(k, G, A_b, A_t, G),
                                 A_t.real, abs_tol=1e-12)


def test_planar_radial_trace_matches_finite_difference():
    """planar_radial_trace matches centred FD of the field."""
    h = 1e-8
    A_b, A_t = 0.7 + 0.2j, -0.3 + 0.5j
    d0, dg = planar_radial_trace(K_WAVE, G, A_b, A_t)
    fd0 = (planar_field(K_WAVE, G, A_b, A_t, h) - planar_field(K_WAVE, G, A_b, A_t, -h)) / (2 * h)
    fdg = (planar_field(K_WAVE, G, A_b, A_t, G + h) - planar_field(K_WAVE, G, A_b, A_t, G - h)) / (2 * h)
    assert math.isclose(d0.real, fd0, rel_tol=1e-5)
    assert math.isclose(dg.real, fdg, rel_tol=1e-5)


def test_planar_dtn_matrix_explicit():
    """Verify the closed-form DtN matrix M against planar_radial_trace unit responses."""
    k, g = K_WAVE, G
    M = planar_dtn_matrix(k, g)
    # Column 1: A_bottom=1, A_top=0
    d0_in, dg_in = planar_radial_trace(k, g, 1.0, 0.0)
    assert abs(M[0][0] - d0_in) < 1e-12 * abs(d0_in)
    assert abs(M[1][0] - dg_in) < 1e-12 * abs(dg_in)
    # Column 2: A_bottom=0, A_top=1
    d0_out, dg_out = planar_radial_trace(k, g, 0.0, 1.0)
    assert abs(M[0][1] - d0_out) < 1e-12 * abs(d0_out)
    assert abs(M[1][1] - dg_out) < 1e-12 * abs(dg_out)


def test_planar_thrust_is_position_independent():
    """The gem: thrust from planar_harmonic_thrust == maxwell_thrust_numeric at every z."""
    A_b, A_t = 1.0 + 0j, 0.5 * cmath.exp(1j * 0.7)
    F_formula = planar_harmonic_thrust(K_WAVE, G, A_b, A_t, X_PERIOD)
    for frac in (0.0, 0.2, 0.5, 0.8, 1.0):
        z = frac * G
        F_num = maxwell_thrust_numeric(K_WAVE, G, A_b, A_t, z, X_PERIOD)
        assert math.isclose(F_formula, F_num, rel_tol=1e-9, abs_tol=1e-25), (frac, F_formula, F_num)


def test_planar_thrust_zero_in_phase():
    """No thrust when A_bottom and A_top are in phase (real ratio)."""
    F = planar_harmonic_thrust(K_WAVE, G, 1.0 + 0j, 0.5 + 0j, X_PERIOD)
    assert abs(F) < 1e-20


def test_planar_thrust_antisymmetric_phase_sign():
    """Reversing the lead angle reverses the thrust direction."""
    A_b = 1.0 + 0j
    Fp = planar_harmonic_thrust(K_WAVE, G, A_b, 0.5 * cmath.exp(+1j * 0.6), X_PERIOD)
    Fm = planar_harmonic_thrust(K_WAVE, G, A_b, 0.5 * cmath.exp(-1j * 0.6), X_PERIOD)
    assert math.isclose(Fp, -Fm, rel_tol=1e-9)


def test_planar_thin_gap_limit():
    """For k*g << 1: A ≈ linear interpolation and dA/dz ≈ (A_top - A_bottom)/g."""
    k, g_thin = 1.0, 1e-5     # k*g = 1e-5 << 1
    A_b, A_t = 1.0 + 0j, 0.3 + 0j
    d0, dg = planar_radial_trace(k, g_thin, A_b, A_t)
    expected = (A_t - A_b) / g_thin
    assert abs(d0.real - expected.real) / abs(expected.real) < 1e-6
    assert abs(dg.real - expected.real) / abs(expected.real) < 1e-6


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("[OK] AGE analytic gate: harmonic transfer + mesh-INDEPENDENT harmonic torque "
          "(== Maxwell-stress quadrature at every radius). Ready for the NGSolve Schur assembly.")
    print("[OK] Planar DtN gate: field, trace, DtN matrix, position-INDEPENDENT thrust.")


if __name__ == "__main__":
    main()
