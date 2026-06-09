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
                                                    maxwell_torque_numeric)

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


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("[OK] AGE analytic gate: harmonic transfer + mesh-INDEPENDENT harmonic torque "
          "(== Maxwell-stress quadrature at every radius). Ready for the NGSolve Schur assembly.")


if __name__ == "__main__":
    main()
