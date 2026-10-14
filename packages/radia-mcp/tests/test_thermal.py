"""Closed-form thermal post-processors -- regression test.

One focused test per function: an exact closed-form value, ONE non-tautological
independent gate (a limiting case that reproduces a DIFFERENT known result or a
special number / scaling), and the ValueError guard. Pure analytic, no field
solve, no external solver."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.thermal import (
    SIGMA_SB,
    forced_convection_dittus_boelter,
    induction_heating_surface_power,
    lumped_thermal_time_constant,
    natural_convection_vertical_plate,
    radiation_exchange_two_surface,
    thermal_constriction_resistance,
)


def test_lumped_thermal_time_constant():
    k, A, L, h = 200.0, 1e-4, 0.02, 50.0
    rho, c_p, P = 2700.0, 900.0, 5.0
    r = lumped_thermal_time_constant(k, A, L, h, rho, c_p, P)
    r_cond = L / (k * A)
    r_film = 1.0 / (h * A)
    c_th = rho * c_p * A * L
    assert math.isclose(r["R_th"], r_cond + r_film, rel_tol=1e-12)
    assert math.isclose(r["C_th"], c_th, rel_tol=1e-12)
    # gate: steady rise is power * series resistance, and tau = R C exactly
    assert math.isclose(r["dT_steady"], P * (r_cond + r_film), rel_tol=1e-12)
    assert math.isclose(r["tau"], (r_cond + r_film) * c_th, rel_tol=1e-12)
    # gate: a first-order system reaches 1 - 1/e at t = tau (independent of inputs)
    assert math.isclose(r["frac_at_tau"], 1.0 - 1.0 / math.e, rel_tol=1e-12)
    assert 0.632 < r["frac_at_tau"] < 0.633
    with pytest.raises(ValueError):
        lumped_thermal_time_constant(-k, A, L, h, rho, c_p, P)


def test_thermal_constriction_resistance():
    k, a = 400.0, 1e-3
    r = thermal_constriction_resistance(k, a)
    assert math.isclose(r["R_isothermal"], 1.0 / (4.0 * k * a), rel_tol=1e-12)
    assert math.isclose(r["R_isoflux"], 8.0 / (3.0 * math.pi ** 2 * k * a), rel_tol=1e-12)
    # gate: the isoflux/isothermal ratio is the pure number 32/(3 pi^2), no k or a
    assert math.isclose(r["constriction_factor"], 32.0 / (3.0 * math.pi ** 2), rel_tol=1e-12)
    assert "R_two_body" not in r
    # gate: two equal half-spaces -> twice the single-body isothermal resistance
    k2 = 250.0
    r2 = thermal_constriction_resistance(k, a, k2=k2)
    assert math.isclose(r2["R_two_body"], (1.0 / (4.0 * a)) * (1.0 / k + 1.0 / k2), rel_tol=1e-12)
    same = thermal_constriction_resistance(k, a, k2=k)
    assert math.isclose(same["R_two_body"], 2.0 * same["R_isothermal"], rel_tol=1e-12)
    with pytest.raises(ValueError):
        thermal_constriction_resistance(k, a, k2=-1.0)


def test_induction_heating_surface_power():
    sigma, mu_r, f, Hs = 1.0e7, 1.0, 50_000.0, 2000.0
    r = induction_heating_surface_power(sigma, mu_r, f, Hs)
    mu0 = 4.0e-7 * math.pi
    delta = math.sqrt(2.0 / (2.0 * math.pi * f * mu0 * mu_r * sigma))
    assert math.isclose(r["reference_depth"], delta, rel_tol=1e-12)
    assert math.isclose(r["Rs"], 1.0 / (sigma * delta), rel_tol=1e-12)
    # gate: the Rs*Hs^2 form equals the explicit sqrt(pi f mu0 mu_r/sigma) form
    explicit = Hs ** 2 * math.sqrt(math.pi * f * mu0 * mu_r / sigma)
    assert math.isclose(r["power_per_area"], explicit, rel_tol=1e-12)
    # gate: P/A scales as sqrt(f) -- quadruple f -> double the power density
    r4 = induction_heating_surface_power(sigma, mu_r, 4.0 * f, Hs)
    assert math.isclose(r4["power_per_area"], 2.0 * r["power_per_area"], rel_tol=1e-12)
    with pytest.raises(ValueError):
        induction_heating_surface_power(sigma, mu_r, -f, Hs)


def test_forced_convection_dittus_boelter():
    Re, Pr, k_f, D = 1.0e5, 4.0, 0.6, 0.02
    hot = forced_convection_dittus_boelter(Re, Pr, k_f, D, heating=True)
    assert math.isclose(hot["Nu"], 0.023 * Re ** 0.8 * Pr ** 0.4, rel_tol=1e-12)
    assert math.isclose(hot["h"], hot["Nu"] * k_f / D, rel_tol=1e-12)
    assert hot["Nu"] > 0
    # gate: heating exponent 0.4 vs cooling 0.3 -> ratio is exactly Pr^0.1
    cold = forced_convection_dittus_boelter(Re, Pr, k_f, D, heating=False)
    assert math.isclose(cold["Nu"], 0.023 * Re ** 0.8 * Pr ** 0.3, rel_tol=1e-12)
    assert math.isclose(hot["Nu"] / cold["Nu"], Pr ** 0.1, rel_tol=1e-12)
    assert hot["Nu"] > cold["Nu"]              # Pr>1 so heating gives larger Nu
    with pytest.raises(ValueError):
        forced_convection_dittus_boelter(-Re, Pr, k_f, D)


def test_natural_convection_vertical_plate():
    Ra, Pr, k_f, Lh = 1.0e9, 0.7, 0.026, 0.3
    r = natural_convection_vertical_plate(Ra, Pr, k_f, Lh)
    denom = (1.0 + (0.492 / Pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    nu = (0.825 + 0.387 * Ra ** (1.0 / 6.0) / denom) ** 2
    assert math.isclose(r["Nu"], nu, rel_tol=1e-12)
    assert math.isclose(r["h"], nu * k_f / Lh, rel_tol=1e-12)
    # gate: as Ra -> 0 the bracket -> 0.825 so Nu -> 0.825^2 (conduction floor)
    r0 = natural_convection_vertical_plate(0.0, Pr, k_f, Lh)
    assert math.isclose(r0["Nu"], 0.825 ** 2, rel_tol=1e-12)
    tiny = natural_convection_vertical_plate(1e-24, Pr, k_f, Lh)  # Ra^(1/6)=1e-4
    assert math.isclose(tiny["Nu"], 0.825 ** 2, rel_tol=1e-3)
    with pytest.raises(ValueError):
        natural_convection_vertical_plate(Ra, -Pr, k_f, Lh)


def test_radiation_exchange_two_surface():
    e1, A1, T1, T2 = 0.8, 0.05, 400.0, 300.0
    r = radiation_exchange_two_surface(e1, A1, T1, T2)
    # small-object-in-enclosure closed form
    assert math.isclose(r["Q"], e1 * SIGMA_SB * A1 * (T1 ** 4 - T2 ** 4), rel_tol=1e-12)
    # gate: linearised film coefficient is 4 e sigma T_m^3 about the mean T
    t_m = 0.5 * (T1 + T2)
    assert math.isclose(r["h_rad_linearized"], 4.0 * e1 * SIGMA_SB * t_m ** 3, rel_tol=1e-12)
    # gate: net flux vanishes at thermal equilibrium T1 == T2
    eq = radiation_exchange_two_surface(e1, A1, 350.0, 350.0)
    assert eq["Q"] == 0.0
    # gate: scales linearly with sigma_SB (Q / sigma_SB == e1 A1 (T1^4 - T2^4))
    assert math.isclose(r["Q"] / SIGMA_SB, e1 * A1 * (T1 ** 4 - T2 ** 4), rel_tol=1e-12)
    # gate: general two-surface gray form with a finite enclosure reduces toward
    # the small-object limit as the enclosure grows (A2 >> A1, e2 -> 1)
    big = radiation_exchange_two_surface(e1, A1, T1, T2, emissivity_2=1.0,
                                         area_2=1.0e6, view_factor=1.0)
    assert math.isclose(big["Q"], r["Q"], rel_tol=1e-3)
    with pytest.raises(ValueError):
        radiation_exchange_two_surface(1.5, A1, T1, T2)


if __name__ == "__main__":
    test_lumped_thermal_time_constant()
    test_thermal_constriction_resistance()
    test_induction_heating_surface_power()
    test_forced_convection_dittus_boelter()
    test_natural_convection_vertical_plate()
    test_radiation_exchange_two_surface()
    print("[OK] thermal closed forms verified.")
