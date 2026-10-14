"""Surface-impedance (SIBC) plumbing and skin-effect physics, current API.

Rewrite of the 2026-02 script-style demo this file used to hold: that
version drove a `PEECBuilder.set_frequency()` API that no longer exists
(frequency moved into `PEECCircuitSolver.compute_port_impedance(freq,
Zs=...)`), printed instead of asserting, and wrote a PNG into the CWD.
The physics it demonstrated is asserted here against closed forms:

* the DC resistance of a straight bar is exact (rho L / A);
* the Zs plumbing is exact by construction -- for a single segment,
  Z_branch = diag(R + Zs) + jwL, so the port resistance must be
  R_dc + Re(Zs) to machine precision (this is the contract the SIBC
  workflow stands on);
* the frequency-dependent resistance itself comes from the SHIPPED
  closed form `radia.analytical_formulas.cylinder_ac_impedance` (full
  Bessel round-wire solution), pinned against its thin-skin asymptote
  R_ac/R_dc -> a/(2 delta) + 1/4 + O(delta/a).
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.analytical_formulas import (  # noqa: E402
    cylinder_ac_impedance,
    cylinder_dc_resistance,
    dowell_rectangular_ac_impedance,
    skin_depth,
)
from radia.peec_matrices import PEECBuilder  # noqa: E402
from radia.peec_topology import PEECCircuitSolver  # noqa: E402

SIGMA_CU = 5.8e7
RHO_CU = 1.0 / SIGMA_CU


def _bar_solver(length=0.2, w=2e-3, h=2e-3):
    b = PEECBuilder()
    n0 = b.add_node_at(0.0, 0.0, 0.0)
    n1 = b.add_node_at(length, 0.0, 0.0)
    b.add_connected_segment(n0, n1, w, h, sigma=SIGMA_CU)
    b.add_port(n0, n1)
    topo = b.build_topology()
    return PEECCircuitSolver(topo), topo


def test_dc_resistance_is_the_analytic_bar_value():
    length, w, h = 0.2, 2e-3, 2e-3
    solver, topo = _bar_solver(length, w, h)
    r_analytic = RHO_CU * length / (w * h)
    assert float(np.asarray(topo["R"]).sum()) == pytest.approx(r_analytic,
                                                               rel=1e-12)
    z = solver.compute_port_impedance(10.0)
    assert z.real == pytest.approx(r_analytic, rel=1e-9)


def test_zs_adds_to_the_port_resistance_exactly():
    """The SIBC contract: Z_branch = diag(R + Zs) + jwL, so for one
    segment the port resistance is R_dc + Re(Zs) to machine precision
    and the reactive part of Zs adds to the inductive rise."""
    solver, topo = _bar_solver()
    r_dc = float(np.asarray(topo["R"]).sum())
    freq = 1.0e5
    zs = np.array([3.7e-3 + 2.1e-3j])
    z0 = solver.compute_port_impedance(freq)
    z1 = solver.compute_port_impedance(freq, Zs=zs)
    assert z1.real - z0.real == pytest.approx(zs[0].real, rel=1e-9)
    assert z1.imag - z0.imag == pytest.approx(zs[0].imag, rel=1e-9)
    assert z0.real == pytest.approx(r_dc, rel=1e-9)


def test_bessel_round_wire_sibc_through_the_port():
    """Feed the shipped full-Bessel round-wire impedance in as Zs (the
    internal-impedance EXCESS over DC, which is what an SIBC adds on
    top of the builder's DC resistance) and check the port reproduces
    the closed form."""
    length, a = 0.2, 1.0e-3
    b = PEECBuilder()
    n0 = b.add_node_at(0.0, 0.0, 0.0)
    n1 = b.add_node_at(length, 0.0, 0.0)
    # square section with the SAME DC resistance as the round wire
    side = math.sqrt(math.pi) * a
    b.add_connected_segment(n0, n1, side, side, sigma=SIGMA_CU)
    b.add_port(n0, n1)
    solver = PEECCircuitSolver(b.build_topology())

    freq = 1.0e6
    omega = 2 * math.pi * freq
    z_wire = cylinder_ac_impedance(a, SIGMA_CU, omega) * length
    r_dc = cylinder_dc_resistance(a, SIGMA_CU) * length
    zs = np.array([z_wire - r_dc])

    z = solver.compute_port_impedance(freq, Zs=zs)
    assert z.real == pytest.approx(z_wire.real, rel=1e-9)


def test_skin_effect_resistance_matches_the_thin_skin_asymptote():
    """R_ac/R_dc of the shipped Bessel solution vs the classical
    a/(2 delta) + 1/4 asymptote, and monotone growth with frequency --
    the physics the old demo plotted, as assertions."""
    a = 1.0e-3
    r_dc = cylinder_dc_resistance(a, SIGMA_CU)
    prev = 1.0
    for freq in (1.0e4, 1.0e5, 1.0e6, 1.0e7):
        omega = 2 * math.pi * freq
        ratio = cylinder_ac_impedance(a, SIGMA_CU, omega).real / r_dc
        assert ratio > prev * 1.01, (freq, ratio, prev)
        prev = ratio
        delta = skin_depth(SIGMA_CU, omega)
        if a / delta > 6.0:            # asymptote valid for thin skin
            asym = a / (2.0 * delta) + 0.25
            assert ratio == pytest.approx(asym, rel=0.02), (freq, ratio,
                                                            asym)


def test_rectangular_dowell_preserves_dc_and_orientation():
    width, thickness = 10.0e-3, 0.2e-3
    r_dc = 1.0 / (SIGMA_CU * width * thickness)
    assert dowell_rectangular_ac_impedance(
        width, thickness, SIGMA_CU, 0.0) == pytest.approx(r_dc)
    assert dowell_rectangular_ac_impedance(
        thickness, width, SIGMA_CU, 0.0) == pytest.approx(r_dc)


def test_rectangular_dowell_matches_closed_form_resistance_factor():
    width, thickness = 8.0e-3, 1.0e-3
    omega = 2.0 * math.pi * 100.0e3
    delta = skin_depth(SIGMA_CU, omega)
    xi = thickness / (2.0 * delta)
    expected_factor = xi * (
        math.sinh(2.0 * xi) + math.sin(2.0 * xi)
    ) / (math.cosh(2.0 * xi) - math.cos(2.0 * xi))
    r_dc = 1.0 / (SIGMA_CU * width * thickness)
    z = dowell_rectangular_ac_impedance(
        width, thickness, SIGMA_CU, omega)
    assert z.real / r_dc == pytest.approx(expected_factor, rel=1.0e-12)


def test_rectangular_dowell_internal_inductance_and_skin_asymptote():
    width, thickness = 10.0e-3, 0.1e-3
    omega_low = 1.0
    z_low = dowell_rectangular_ac_impedance(
        width, thickness, SIGMA_CU, omega_low)
    expected_l_internal = 4.0e-7 * math.pi * thickness / (12.0 * width)
    assert z_low.imag / omega_low == pytest.approx(
        expected_l_internal, rel=1.0e-10)

    omega_high = 2.0 * math.pi * 1.0e9
    delta = skin_depth(SIGMA_CU, omega_high)
    expected_surface = (1.0 + 1.0j) / (
        2.0 * SIGMA_CU * width * delta)
    z_high = dowell_rectangular_ac_impedance(
        width, thickness, SIGMA_CU, omega_high)
    assert z_high == pytest.approx(expected_surface, rel=1.0e-12)


@pytest.mark.parametrize("width,thickness", [(0.0, 1.0), (1.0, -1.0)])
def test_rectangular_dowell_rejects_nonphysical_dimensions(width, thickness):
    with pytest.raises(ValueError):
        dowell_rectangular_ac_impedance(
            width, thickness, SIGMA_CU, 1.0)
