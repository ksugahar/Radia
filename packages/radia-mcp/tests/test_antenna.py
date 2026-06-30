r"""Antenna closed-form helpers -- regression test.

Pure textbook far-field figures of merit (radiation resistance, directivity,
array factor, aperture pattern). Each test pins the exact closed-form value at a
concrete point, adds one non-tautological independent cross-check (a limiting
case matching a DIFFERENT known result / special value / scaling law), and
exercises the ValueError path."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.antenna import (
    CIN_2PI,
    ETA0,
    circular_aperture_pattern,
    farfield_gain_directivity_efficiency_gate,
    finite_dipole_directivity,
    half_wave_dipole_impedance,
    short_dipole_radiation_resistance,
    small_loop_radiation_resistance,
    uniform_linear_array_factor,
)


def test_short_dipole_radiation_resistance():
    # exact closed form: R = 20 pi^2 (L/lam)^2; L/lam = 0.1 -> ~1.974 Ohm
    r = short_dipole_radiation_resistance(0.1)
    assert math.isclose(r["R_rad"], 20.0 * math.pi ** 2 * 0.01, rel_tol=1e-12)
    assert math.isclose(r["R_rad"], 1.9739208802, rel_tol=1e-8)
    # directivity of the elementary electric dipole is exactly 3/2 (sin^2 pattern)
    assert r["directivity"] == 1.5
    # independent cross-check: a 4x longer short dipole radiates 16x the power
    # (R scales as length^2) -- a scaling law, not a restatement of the formula
    r2 = short_dipole_radiation_resistance(0.4)
    assert math.isclose(r2["R_rad"] / r["R_rad"], 16.0, rel_tol=1e-12)
    with pytest.raises(ValueError):
        short_dipole_radiation_resistance(0.0)


def test_small_loop_radiation_resistance():
    # exact closed form: R = 20 pi^2 (C/lam)^4 n^2; C/lam = 0.1, n=1
    r = small_loop_radiation_resistance(0.1)
    assert math.isclose(r["R_rad"], 20.0 * math.pi ** 2 * 0.1 ** 4, rel_tol=1e-12)
    assert math.isclose(r["R_rad"], 0.0197392088, rel_tol=1e-7)
    # the small loop is the magnetic dual of the electric dipole: same D0 = 3/2
    assert r["directivity"] == 1.5
    # independent cross-check: n turns multiply R by n^2 (current-moment scaling)
    r3 = small_loop_radiation_resistance(0.1, n_turns=3)
    assert math.isclose(r3["R_rad"] / r["R_rad"], 9.0, rel_tol=1e-12)
    # independent cross-check vs a DIFFERENT formula: loop R is the 4th-power
    # radiator, dipole R the 2nd-power one; at C/lam = L/lam the ratio is (C/lam)^2
    x = 0.05
    loop = small_loop_radiation_resistance(x)["R_rad"]
    dip = short_dipole_radiation_resistance(x)["R_rad"]
    assert math.isclose(loop / dip, x ** 2, rel_tol=1e-12)
    with pytest.raises(ValueError):
        small_loop_radiation_resistance(0.1, n_turns=0)


def test_half_wave_dipole_impedance():
    z = half_wave_dipole_impedance()
    # exact closed form: R = eta0/(4 pi) * Cin(2 pi) ~= 73.08 Ohm
    assert math.isclose(z["R_rad"], (ETA0 / (4.0 * math.pi)) * CIN_2PI, rel_tol=1e-12)
    assert math.isclose(z["R_rad"], 73.08, rel_tol=2e-4)
    assert math.isclose(z["X"], 42.5, rel_tol=1e-12)
    assert z["Z"] == complex(z["R_rad"], z["X"])
    # independent cross-check: directivity D0 = 4/Cin(2 pi) = 1.6409 (a value
    # derived separately from the cosine integral, not from R)
    assert math.isclose(z["directivity"], 4.0 / CIN_2PI, rel_tol=1e-12)
    assert math.isclose(z["directivity"], 1.640922, rel_tol=1e-5)


def test_finite_dipole_directivity():
    # half-wave dipole: numerically integrated D0 must reproduce 1.6409
    d = finite_dipole_directivity(0.5)
    assert math.isclose(d["directivity"], 1.6409, abs_tol=1e-3)
    # independent cross-check: the numerical-integration directivity at L/lam=0.5
    # matches the *closed-form* half-wave value 4/Cin(2 pi) to < 1e-3
    hw = half_wave_dipole_impedance()["directivity"]
    assert abs(d["directivity"] - hw) < 1e-3
    # dBi consistency (a separate quantity): D0 = 10^(dBi/10)
    assert math.isclose(d["directivity"], 10.0 ** (d["directivity_dBi"] / 10.0),
                        rel_tol=1e-12)
    with pytest.raises(ValueError):
        finite_dipole_directivity(-0.5)


def test_uniform_linear_array_factor():
    N = 8
    d = 0.5
    # main beam (broadside, theta = 90 deg, scan = 90 deg): psi = 0 -> L'Hopital
    main = uniform_linear_array_factor(N, d, 90.0, scan_deg=90.0)
    assert math.isclose(main["psi"], 0.0, abs_tol=1e-12)
    assert math.isclose(main["AF"], float(N), rel_tol=1e-12)       # |AF| = N
    assert math.isclose(main["AF_normalized"], 1.0, rel_tol=1e-12)  # AF/N = 1
    # independent cross-check vs a DIFFERENT known fact: the array factor has an
    # EXACT null where N psi/2 = pi, i.e. psi = 2 pi/N. For broadside, psi =
    # 2 pi d (cos theta) = 2 pi/N -> cos theta = 1/(N d). Construct that angle.
    cos_null = 1.0 / (N * d)
    theta_null = math.degrees(math.acos(cos_null))
    nul = uniform_linear_array_factor(N, d, theta_null, scan_deg=90.0)
    assert math.isclose(nul["psi"], 2.0 * math.pi / N, rel_tol=1e-9)
    assert math.isclose(nul["AF"], 0.0, abs_tol=1e-7)
    assert math.isclose(nul["AF_normalized"], 0.0, abs_tol=1e-8)
    with pytest.raises(ValueError):
        uniform_linear_array_factor(0, d, 90.0)


def test_circular_aperture_pattern():
    Dlam = 10.0
    # boresight: E(0) = 1 exactly (2 J1(u)/u -> 1 as u -> 0)
    bore = circular_aperture_pattern(Dlam, 0.0)
    assert math.isclose(bore["E_normalized"], 1.0, rel_tol=1e-12)
    assert math.isclose(bore["u"], 0.0, abs_tol=1e-12)
    # uniform-aperture directivity D0 = (pi D/lam)^2
    assert math.isclose(bore["directivity_uniform"], (math.pi * Dlam) ** 2, rel_tol=1e-12)
    # independent cross-check: the FIRST NULL is the first zero of J1 (u=3.8317),
    # i.e. sin(theta_null) = 3.8317 / (pi D/lam). Evaluate E there -> ~0.
    U1 = 3.8317059702
    sin_null = U1 / (math.pi * Dlam)
    theta_null_deg = math.degrees(math.asin(sin_null))
    assert math.isclose(bore["first_null_deg"], theta_null_deg, rel_tol=1e-9)
    nul = circular_aperture_pattern(Dlam, theta_null_deg)
    assert math.isclose(nul["u"], U1, rel_tol=1e-9)
    assert abs(nul["E_normalized"]) < 1e-6
    with pytest.raises(ValueError):
        circular_aperture_pattern(0.0, 0.0)


def test_farfield_gain_directivity_efficiency_gate():
    directivity_dbi = 6.0
    eta = 0.72
    gain_dbi = 10.0 * math.log10((10.0 ** (directivity_dbi / 10.0)) * eta)
    got = farfield_gain_directivity_efficiency_gate(
        gain_dbi=gain_dbi,
        directivity_dbi=directivity_dbi,
        radiated_power_w=7.2,
        accepted_power_w=10.0,
    )

    assert got["status"] == "ok"
    assert got["radiation_efficiency"] == pytest.approx(0.72)
    assert got["efficiency_from_gain_over_directivity"] == pytest.approx(0.72)
    assert got["gain_relative_error"] < 1.0e-12
    assert all(got["checks"].values())

    too_high_gain = farfield_gain_directivity_efficiency_gate(
        gain_dbi=directivity_dbi + 0.1,
        directivity_dbi=directivity_dbi,
        radiated_power_w=7.2,
        accepted_power_w=10.0,
    )
    assert too_high_gain["status"] == "needs_attention"
    assert too_high_gain["checks"]["gain_not_above_directivity"] is False

    inconsistent_power = farfield_gain_directivity_efficiency_gate(
        gain_dbi=gain_dbi,
        directivity_dbi=directivity_dbi,
        radiated_power_w=5.0,
        accepted_power_w=10.0,
    )
    assert inconsistent_power["status"] == "needs_attention"
    assert inconsistent_power["checks"]["gain_matches_directivity_times_efficiency"] is False

    with pytest.raises(ValueError):
        farfield_gain_directivity_efficiency_gate(gain_dbi, directivity_dbi, 1.0, 0.0)


if __name__ == "__main__":
    test_short_dipole_radiation_resistance()
    test_small_loop_radiation_resistance()
    test_half_wave_dipole_impedance()
    test_finite_dipole_directivity()
    test_uniform_linear_array_factor()
    test_circular_aperture_pattern()
    test_farfield_gain_directivity_efficiency_gate()
    print("[OK] antenna closed forms: dipole/loop R_rad, D0, array factor, aperture.")
