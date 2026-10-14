r"""TEM transmission-line per-unit-length parameters -- closed-form regression test.

Each test pins the exact textbook value at a concrete point, adds ONE
non-tautological independent cross-check (a limiting case that matches a
DIFFERENT known result), and exercises the ValueError path. The recurring
independent gate is the geometry-cancelling identity L*C == mu0*mu_r*eps0*eps_r
(so vp == c/sqrt(eps_r*mu_r)), which the C and L closed forms must satisfy even
though each carries the geometry factor separately."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.transmission_line import (
    EPS0,
    MU0,
    C0,
    coaxial_line_parameters,
    tem_lc_identity_summary,
    two_wire_line_parameters,
    microstrip_line_parameters,
    wire_capacitance_per_length,
)


def test_coaxial_line_parameters():
    a, b, eps_r, mu_r = 1.0e-3, 3.5e-3, 2.1, 1.0
    r = coaxial_line_parameters(a, b, eps_r, mu_r)

    ln_ratio = math.log(b / a)
    assert math.isclose(r["ln_ratio"], ln_ratio, rel_tol=1e-12)
    assert math.isclose(r["C_per_m"], 2 * math.pi * EPS0 * eps_r / ln_ratio, rel_tol=1e-12)
    assert math.isclose(r["L_per_m"], MU0 * mu_r * ln_ratio / (2 * math.pi), rel_tol=1e-12)
    # exact value: Z0 = (eta0/(2*pi)) sqrt(mu_r/eps_r) ln(b/a) ~ 51.8 ohm
    z0_exact = (MU0 * C0 / (2 * math.pi)) * math.sqrt(mu_r / eps_r) * ln_ratio
    assert math.isclose(r["Z0"], z0_exact, rel_tol=1e-9)
    assert math.isclose(r["Z0"], 51.83, rel_tol=1e-3)

    # non-tautological gate: geometry cancels -> L*C == mu0*mu_r*eps0*eps_r exactly,
    # so vp == c/sqrt(eps_r*mu_r) regardless of b/a
    assert math.isclose(r["L_per_m"] * r["C_per_m"], MU0 * mu_r * EPS0 * eps_r, rel_tol=1e-12)
    assert math.isclose(r["vp"], C0 / math.sqrt(eps_r * mu_r), rel_tol=1e-9)

    with pytest.raises(ValueError):
        coaxial_line_parameters(3.0e-3, 1.0e-3)  # outer <= inner


def test_tem_lc_identity_summary_matches_geometry_helpers():
    coax = coaxial_line_parameters(0.8e-3, 3.2e-3, eps_r=2.2, mu_r=1.0)
    ident = tem_lc_identity_summary(coax["C_per_m"], coax["L_per_m"], eps_r=2.2, mu_r=1.0)
    assert math.isclose(ident["Z0"], coax["Z0"], rel_tol=1e-12)
    assert math.isclose(ident["vp"], coax["vp"], rel_tol=1e-12)
    assert abs(ident["LC_relative_error"]) < 1e-12
    assert abs(ident["vp_relative_error"]) < 1e-9

    two = two_wire_line_parameters(0.4e-3, 5.0e-3, eps_r=1.7, mu_r=1.2)
    ident_two = tem_lc_identity_summary(two["C_per_m"], two["L_per_m"], eps_r=1.7, mu_r=1.2)
    assert math.isclose(ident_two["LC_product"], MU0 * 1.2 * EPS0 * 1.7, rel_tol=1e-12)
    assert math.isclose(ident_two["vp_expected"], C0 / math.sqrt(1.7 * 1.2), rel_tol=1e-12)

    with pytest.raises(ValueError):
        tem_lc_identity_summary(0.0, coax["L_per_m"])
    with pytest.raises(ValueError):
        tem_lc_identity_summary(coax["C_per_m"], coax["L_per_m"], eps_r=0.0)


def test_two_wire_line_parameters():
    a, D, eps_r, mu_r = 0.5e-3, 6.0e-3, 1.0, 1.0
    r = two_wire_line_parameters(a, D, eps_r, mu_r)

    arg = D / (2 * a)
    assert math.isclose(r["acosh_arg"], arg, rel_tol=1e-12)
    assert math.isclose(r["C_per_m"], math.pi * EPS0 * eps_r / math.acosh(arg), rel_tol=1e-12)
    assert math.isclose(r["L_per_m"], MU0 * mu_r * math.acosh(arg) / math.pi, rel_tol=1e-12)
    # exact value: Z0 = (eta0/pi)*acosh(D/2a) in air (sqrt(L/C) vs direct
    # product agree to floating-point round-off)
    assert math.isclose(r["Z0"], (MU0 * C0 / math.pi) * math.acosh(arg), rel_tol=1e-9)

    # non-tautological gate 1: geometry cancels -> L*C == mu0*eps0 exactly (air line)
    assert math.isclose(r["L_per_m"] * r["C_per_m"], MU0 * mu_r * EPS0 * eps_r, rel_tol=1e-12)
    assert math.isclose(r["vp"], C0, rel_tol=1e-9)

    # non-tautological gate 2: wide-spacing limit acosh(x) -> ln(2x),
    # so a far-separated two-wire line approaches the logarithmic form
    aw, Dw = 1.0e-4, 1.0  # D/(2a) = 5000, well into the asymptote
    rw = two_wire_line_parameters(aw, Dw)
    x = Dw / (2 * aw)
    z0_log = (MU0 * C0 / math.pi) * math.log(2 * x)
    assert math.isclose(rw["Z0"], z0_log, rel_tol=1e-6)

    with pytest.raises(ValueError):
        two_wire_line_parameters(0.5e-3, 1.0e-3)  # D <= 2a (overlap)


def test_microstrip_line_parameters():
    # wide strip on eps_r=4 substrate, u = w/h = 2
    eps_r = 4.0
    r = microstrip_line_parameters(2.0e-3, 1.0e-3, eps_r)
    u = 2.0
    eps_eff_exact = (eps_r + 1) / 2 + (eps_r - 1) / 2 * (1 + 12 / u) ** -0.5
    assert math.isclose(r["u"], u, rel_tol=1e-12)
    assert math.isclose(r["eps_eff"], eps_eff_exact, rel_tol=1e-12)
    z0_exact = (120 * math.pi / math.sqrt(eps_eff_exact)) / (
        u + 1.393 + 0.667 * math.log(u + 1.444)
    )
    assert math.isclose(r["Z0"], z0_exact, rel_tol=1e-12)
    assert math.isclose(r["vp"], C0 / math.sqrt(eps_eff_exact), rel_tol=1e-12)
    # effective permittivity is strictly between air (1) and the bulk substrate
    assert 1.0 < r["eps_eff"] < eps_r

    # non-tautological gate 1: an air "substrate" (eps_r=1) is dispersionless,
    # eps_eff -> 1 and the line carries a pure TEM wave at c
    r_air = microstrip_line_parameters(2.0e-3, 1.0e-3, 1.0)
    assert math.isclose(r_air["eps_eff"], 1.0, rel_tol=1e-12)
    assert math.isclose(r_air["vp"], C0, rel_tol=1e-12)

    # non-tautological gate 2: a very wide strip (u -> inf) fully fills the
    # substrate, so eps_eff -> eps_r (parallel-plate limit). Convergence is
    # ~1/sqrt(u), so probe deep into the asymptote.
    r_wide = microstrip_line_parameters(1.0e6, 1.0, eps_r)  # u = 1e6
    assert math.isclose(r_wide["eps_eff"], eps_r, rel_tol=1e-2)
    assert r_wide["eps_eff"] < eps_r  # approaches from below

    with pytest.raises(ValueError):
        microstrip_line_parameters(1.0e-3, 1.0e-3, 0.5)  # eps_r < 1


def test_wire_capacitance_per_length():
    # exact value: coax form matches 2*pi*eps0/ln(b/a)
    r_coax = wire_capacitance_per_length("coax", 1.0e-3, 3.0e-3, eps_r=1.0)
    assert math.isclose(
        r_coax["C_per_m"], 2 * math.pi * EPS0 / math.log(3.0), rel_tol=1e-12
    )
    assert r_coax["geometry"] == "coax"

    # consistency: 'coax' geometry == coaxial_line_parameters C
    assert math.isclose(
        r_coax["C_per_m"],
        coaxial_line_parameters(1.0e-3, 3.0e-3)["C_per_m"],
        rel_tol=1e-12,
    )

    # non-tautological gate: image symmetry. A wire at height h above a ground
    # plane equals HALF a two-wire line of separation D = 2h, hence carries
    # TWICE its capacitance: C_plane(h) == 2 * C_two_wire(D = 2h).
    radius, h = 0.5e-3, 5.0e-3
    c_plane = wire_capacitance_per_length("wire_plane", radius, h)["C_per_m"]
    c_two = wire_capacitance_per_length("two_wire", radius, 2.0 * h)["C_per_m"]
    assert math.isclose(c_plane, 2.0 * c_two, rel_tol=1e-12)

    with pytest.raises(ValueError):
        wire_capacitance_per_length("stripline", 1.0e-3, 3.0e-3)  # unknown geometry


if __name__ == "__main__":
    test_coaxial_line_parameters()
    test_two_wire_line_parameters()
    test_microstrip_line_parameters()
    test_wire_capacitance_per_length()
    print("[OK] TEM transmission-line per-unit-length parameters: closed forms verified.")
