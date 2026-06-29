r"""Electrostatics closed-form regression tests (capacitance / conduction /
polarizability / field stress).  Each test checks the exact textbook value plus
one independent, non-tautological gate -- a limiting case that must reproduce a
DIFFERENT known result (isolated sphere, RC duality, conducting-sphere
saturation, sign flips, etc.) -- and the ValueError guard.  Pure analytic
reference (no external solver)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.electrostatics import (
    EPS0,
    isolated_disk_capacitance,
    spreading_resistance_disk,
    coaxial_shell_resistance,
    coaxial_capacitor_energy_force,
    layered_parallel_plate_capacitance,
    parallel_plate_capacitor_energy_force,
    capacitance_gradient_force_summary,
    dielectric_sphere_polarizability,
    dielectric_sphere_uniform_field,
    uniformly_polarized_sphere_field,
    charged_conducting_sphere_surface_stress,
    sphere_above_plane_capacitance,
    dielectrophoresis_force,
)
from radia_mcp.radia_ngsolve.force import electrostatic_traction_summary


def test_isolated_disk_capacitance():
    a = 0.03
    C = isolated_disk_capacitance(a)["C"]
    assert math.isclose(C, 8.0 * EPS0 * a, rel_tol=1e-12)
    # GATE: ratio to the isolated SPHERE of equal radius (C_sphere = 4 pi eps0 a) is exactly 2/pi
    C_sphere = 4.0 * math.pi * EPS0 * a
    assert math.isclose(C / C_sphere, 2.0 / math.pi, rel_tol=1e-12)
    # dielectric scaling
    assert math.isclose(isolated_disk_capacitance(a, eps_r=3.0)["C"], 3.0 * C, rel_tol=1e-12)
    with pytest.raises(ValueError):
        isolated_disk_capacitance(0.0)


def test_spreading_resistance_disk():
    rho, a = 1.7e-8, 1e-4
    out = spreading_resistance_disk(rho, a)
    assert math.isclose(out["R"], rho / (4.0 * a), rel_tol=1e-12)
    assert math.isclose(out["hemispherical_R"], rho / (2.0 * math.pi * a), rel_tol=1e-12)
    # GATE: the flat-spot/hemisphere ratio is the pure geometric pi/2 (independent of rho, a)
    assert math.isclose(out["R"] / out["hemispherical_R"], math.pi / 2.0, rel_tol=1e-12)
    with pytest.raises(ValueError):
        spreading_resistance_disk(rho, -1.0)


def test_coaxial_shell_resistance():
    rho, a, b, L = 2.0e-3, 0.01, 0.02, 0.5
    R = coaxial_shell_resistance(rho, a, b, L)["R"]
    assert math.isclose(R, rho * math.log(b / a) / (2.0 * math.pi * L), rel_tol=1e-12)
    # GATE: RC duality -- coax capacitance C = 2 pi eps0 L / ln(b/a) (computed independently);
    # the ln(b/a) cancels so R*C must equal rho*eps0 exactly.
    C_coax = 2.0 * math.pi * EPS0 * L / math.log(b / a)
    assert math.isclose(R * C_coax, rho * EPS0, rel_tol=1e-12)
    with pytest.raises(ValueError):
        coaxial_shell_resistance(rho, 0.02, 0.01, L)   # r_outer < r_inner


def test_coaxial_capacitor_energy_force():
    er, a, b, length, voltage = 2.5, 0.01, 0.03, 0.2, 120.0
    out = coaxial_capacitor_energy_force(er, a, b, length, voltage)
    eps = EPS0 * er
    log_ratio = math.log(b / a)
    capacitance = 2.0 * math.pi * eps * length / log_ratio
    e_inner = voltage / (a * log_ratio)
    e_outer = voltage / (b * log_ratio)

    assert out["C"] == pytest.approx(capacitance)
    assert out["energy"] == pytest.approx(0.5 * capacitance * voltage * voltage)
    assert out["electric_field_inner_V_per_m"] == pytest.approx(e_inner)
    assert out["electric_field_outer_V_per_m"] == pytest.approx(e_outer)
    assert out["pressure_inner_Pa"] == pytest.approx(0.5 * eps * e_inner * e_inner)
    assert out["pressure_outer_Pa"] == pytest.approx(0.5 * eps * e_outer * e_outer)
    assert out["inner_radius_force_N"] == pytest.approx(out["inner_pressure_area_force_N"])
    assert out["outer_radius_force_N"] == pytest.approx(out["outer_pressure_area_force_N"])

    inner_gradient = capacitance_gradient_force_summary(
        out["C"],
        out["dCdr_inner_F_per_m"],
        voltage_V=voltage,
    )
    outer_gradient = capacitance_gradient_force_summary(
        out["C"],
        out["dCdr_outer_F_per_m"],
        voltage_V=voltage,
    )
    assert inner_gradient["fixed_voltage_force_N"] == pytest.approx(out["inner_radius_force_N"])
    assert outer_gradient["fixed_voltage_force_N"] == pytest.approx(out["outer_radius_force_N"])
    assert out["inner_radius_force_N"] > 0.0
    assert out["outer_radius_force_N"] < 0.0
    assert out["pressure_inner_Pa"] / out["pressure_outer_Pa"] == pytest.approx((b / a) ** 2)

    with pytest.raises(ValueError):
        coaxial_capacitor_energy_force(er, b, a, length, voltage)


def test_layered_parallel_plate_capacitance():
    area, d = 1e-4, 1e-3
    # single layer must reduce to the plain slab C = eps0 eps_r A / d
    one = layered_parallel_plate_capacitance(area, [d], [4.0])
    assert math.isclose(one["C"], EPS0 * 4.0 * area / d, rel_tol=1e-12)
    assert math.isclose(one["eps_eff"], 4.0, rel_tol=1e-12)
    # GATE: two IDENTICAL layers -> eps_eff is unchanged (== eps_r), and total C halves vs one layer
    two = layered_parallel_plate_capacitance(area, [d, d], [4.0, 4.0])
    assert math.isclose(two["eps_eff"], 4.0, rel_tol=1e-12)
    assert math.isclose(two["C"], 0.5 * one["C"], rel_tol=1e-12)
    # series stack with contrasting layers: eps_eff lies between the two
    mix = layered_parallel_plate_capacitance(area, [d, d], [2.0, 8.0])
    assert 2.0 < mix["eps_eff"] < 8.0
    assert math.isclose(mix["eps_eff"], 2.0 / (1.0 / 2.0 + 1.0 / 8.0), rel_tol=1e-12)
    with pytest.raises(ValueError):
        layered_parallel_plate_capacitance(area, [d, d], [4.0])   # length mismatch


def test_parallel_plate_capacitor_energy_force():
    eps_r, area, gap, V = 1.0, 1e-3, 1e-3, 100.0
    out = parallel_plate_capacitor_energy_force(eps_r, area, gap, V)
    C = EPS0 * eps_r * area / gap
    E = V / gap
    pressure = 0.5 * EPS0 * eps_r * E * E
    assert math.isclose(out["C"], C, rel_tol=1e-12)
    assert math.isclose(out["energy"], 0.5 * C * V * V, rel_tol=1e-12)
    assert math.isclose(out["electric_field_V_per_m"], E, rel_tol=1e-12)
    assert math.isclose(out["pressure_Pa"], pressure, rel_tol=1e-12)
    assert math.isclose(out["energy_density_J_per_m3"], pressure, rel_tol=1e-12)
    # GATE: identity |F| == energy/gap (force = energy density * area)
    assert math.isclose(out["force"], out["energy"] / gap, rel_tol=1e-12)
    assert math.isclose(out["force"], pressure * area, rel_tol=1e-12)
    # GATE: same pressure from the electrostatic Maxwell-stress traction helper
    traction = electrostatic_traction_summary(
        (0.0, 0.0, E),
        (0.0, 0.0, 1.0),
        area_m2=area,
        eps=EPS0 * eps_r,
    )
    assert math.isclose(traction["normal_traction_Pa"], out["pressure_Pa"], rel_tol=1e-12)
    assert math.isclose(traction["force_N"][2], out["force"], rel_tol=1e-12)
    # GATE: 1/gap^2 scaling -- halving the gap quadruples the force
    out2 = parallel_plate_capacitor_energy_force(eps_r, area, gap / 2.0, V)
    assert math.isclose(out2["force"], 4.0 * out["force"], rel_tol=1e-12)
    with pytest.raises(ValueError):
        parallel_plate_capacitor_energy_force(eps_r, area, 0.0, V)


def test_parallel_plate_gap_sweep_scaling():
    eps_r, width_per_m, voltage = 1.0, 0.005, 100.0
    gaps = [0.0005, 0.001, 0.002]
    rows = [
        parallel_plate_capacitor_energy_force(eps_r, width_per_m, gap, voltage)
        for gap in gaps
    ]

    fields = [row["electric_field_V_per_m"] for row in rows]
    caps = [row["C"] for row in rows]
    forces = [row["force"] for row in rows]

    assert fields[0] == pytest.approx(2.0 * fields[1])
    assert fields[1] == pytest.approx(2.0 * fields[2])
    assert caps[0] == pytest.approx(2.0 * caps[1])
    assert caps[1] == pytest.approx(2.0 * caps[2])
    assert forces[0] == pytest.approx(4.0 * forces[1])
    assert forces[1] == pytest.approx(4.0 * forces[2])
    for gap, row in zip(gaps, rows):
        assert row["force"] == pytest.approx(row["energy"] / gap)


def test_capacitance_gradient_force_matches_parallel_plate_gap_force():
    eps_r, area, gap, voltage = 2.5, 2.0e-4, 5.0e-4, 120.0
    plate = parallel_plate_capacitor_energy_force(eps_r, area, gap, voltage)
    dcdgap = -plate["C"] / gap
    charge = plate["C"] * voltage
    summary = capacitance_gradient_force_summary(
        plate["C"],
        dcdgap,
        voltage_V=voltage,
        charge_C=charge,
    )

    assert summary["fixed_voltage_force_N"] == pytest.approx(-plate["force"])
    assert summary["fixed_charge_force_N"] == pytest.approx(-plate["force"])
    assert summary["fixed_voltage_force_N"] == pytest.approx(-plate["energy"] / gap)
    assert summary["charge_for_voltage_C"] == pytest.approx(charge)
    assert summary["charge_consistency_error_C"] == pytest.approx(0.0)


def test_capacitance_gradient_force_rejects_bad_inputs():
    with pytest.raises(ValueError):
        capacitance_gradient_force_summary(0.0, 1.0, voltage_V=1.0)
    with pytest.raises(ValueError):
        capacitance_gradient_force_summary(1.0, 1.0)


def test_dielectric_sphere_polarizability():
    a = 0.02
    out = dielectric_sphere_polarizability(a, 5.0)
    cm = (5.0 - 1.0) / (5.0 + 2.0)
    assert math.isclose(out["alpha"], 4.0 * math.pi * EPS0 * a ** 3 * cm, rel_tol=1e-12)
    assert math.isclose(out["M_avg_per_E0"], 3.0 * cm, rel_tol=1e-12)
    # GATE: no contrast (eps_r = 1) -> zero polarizability
    assert dielectric_sphere_polarizability(a, 1.0)["alpha"] == 0.0
    # GATE: conducting limit eps_r -> infinity saturates to alpha = 4 pi eps0 a^3
    big = dielectric_sphere_polarizability(a, 1e9)["alpha"]
    assert math.isclose(big, 4.0 * math.pi * EPS0 * a ** 3, rel_tol=1e-6)
    with pytest.raises(ValueError):
        dielectric_sphere_polarizability(a, 0.0)


def test_dielectric_sphere_uniform_field():
    a, er, E0 = 0.01, 5.0, 1.0
    out = dielectric_sphere_uniform_field(a, er, E0, sample_radius=3.0 * a)
    cm = (er - 1.0) / (er + 2.0)
    assert math.isclose(out["clausius_mossotti"], cm, rel_tol=1e-12)
    assert math.isclose(out["interior_field"], 3.0 * E0 / (er + 2.0), rel_tol=1e-12)
    assert math.isclose(out["polarizability"],
                        dielectric_sphere_polarizability(a, er)["alpha"], rel_tol=1e-12)
    assert math.isclose(out["dipole_moment"], out["polarizability"] * E0, rel_tol=1e-12)
    assert math.isclose(out["axial_field"], E0 * (1.0 + 2.0 * cm / 27.0), rel_tol=1e-12)
    assert math.isclose(out["equatorial_field"], E0 * (1.0 - cm / 27.0), rel_tol=1e-12)

    vacuum = dielectric_sphere_uniform_field(a, 1.0, E0, sample_radius=2.0 * a)
    assert vacuum["interior_field"] == E0
    assert vacuum["axial_field"] == E0
    assert vacuum["equatorial_field"] == E0
    assert dielectric_sphere_uniform_field(a, 1e9)["interior_field"] < 1e-8 * E0

    with pytest.raises(ValueError):
        dielectric_sphere_uniform_field(a, er, E0, sample_radius=a)


def test_uniformly_polarized_sphere_field():
    P, a = 1.3, 0.01
    out = uniformly_polarized_sphere_field(P, a)
    # GATE: interior depolarizing field is exactly -P/(3 eps0) (depolarization factor 1/3)
    assert math.isclose(out["E_in"], -P / (3.0 * EPS0), rel_tol=1e-12)
    # consistency: D_in = eps0 E_in + P = (2/3) P
    assert math.isclose(out["D_in"], (2.0 / 3.0) * P, rel_tol=1e-12)
    assert math.isclose(out["D_in"], EPS0 * out["E_in"] + P, rel_tol=1e-12)
    # exterior dipole moment = P * sphere volume
    assert math.isclose(out["dipole_moment"], P * (4.0 / 3.0) * math.pi * a ** 3, rel_tol=1e-12)
    with pytest.raises(ValueError):
        uniformly_polarized_sphere_field(P, 0.0)


def test_charged_conducting_sphere_surface_stress():
    Q, a = 1e-9, 0.05
    out = charged_conducting_sphere_surface_stress(Q, a)
    E = Q / (4.0 * math.pi * EPS0 * a ** 2)
    assert math.isclose(out["surface_field"], E, rel_tol=1e-12)
    # GATE: Maxwell stress p = eps0 E^2 / 2 equals sigma^2/(2 eps0) with sigma = Q/(4 pi a^2)
    sigma = Q / (4.0 * math.pi * a ** 2)
    assert math.isclose(out["pressure"], 0.5 * EPS0 * E * E, rel_tol=1e-12)
    assert math.isclose(out["pressure"], sigma * sigma / (2.0 * EPS0), rel_tol=1e-12)
    # pressure is outward (positive) regardless of charge sign
    assert charged_conducting_sphere_surface_stress(-Q, a)["pressure"] > 0
    with pytest.raises(ValueError):
        charged_conducting_sphere_surface_stress(Q, 0.0)


def test_sphere_above_plane_capacitance():
    a = 0.01
    C_iso = 4.0 * math.pi * EPS0 * a
    # GATE: far from the plane (h >> a) the image series collapses to the isolated sphere 4 pi eps0 a.
    # The leading proximity correction is 1/(2 cosh(alpha)) ~ a/(2h), so at h = 1e4 a it is ~5e-5.
    far = sphere_above_plane_capacitance(a, 1e4 * a)["C"]
    assert far > C_iso                                    # proximity always raises C
    assert math.isclose(far, C_iso, rel_tol=1e-3)         # but converges to the isolated value
    # closer in, proximity to the ground plane raises the capacitance further above the isolated value
    near = sphere_above_plane_capacitance(a, 1.5 * a)["C"]
    assert near > far > C_iso
    # explicit two-term check of the series at a chosen height
    out = sphere_above_plane_capacitance(a, 2.0 * a, n_terms=2)
    alpha = math.acosh(2.0)
    expect = 4.0 * math.pi * EPS0 * a * (1.0 + math.sinh(alpha) / math.sinh(2.0 * alpha))
    assert math.isclose(out["C"], expect, rel_tol=1e-12)
    assert out["n_terms"] == 2
    with pytest.raises(ValueError):
        sphere_above_plane_capacitance(a, 0.5 * a)   # sphere intersects the plane (h <= a)


def test_dielectrophoresis_force():
    a, em, ep, g = 5e-6, 2.0, 8.0, 1e10
    out = dielectrophoresis_force(a, em, ep, g)
    cm = (ep - em) / (ep + 2.0 * em)
    assert math.isclose(out["clausius_mossotti"], cm, rel_tol=1e-12)
    assert math.isclose(out["force"], 2.0 * math.pi * EPS0 * em * a ** 3 * cm * g, rel_tol=1e-12)
    # GATE: Clausius-Mossotti sign flips across eps_p = eps_m (positive vs negative DEP)
    assert dielectrophoresis_force(a, em, ep, g)["clausius_mossotti"] > 0          # eps_p > eps_m
    assert dielectrophoresis_force(a, em, 1.0, g)["clausius_mossotti"] < 0          # eps_p < eps_m
    assert dielectrophoresis_force(a, em, em, g)["clausius_mossotti"] == 0.0        # no contrast
    with pytest.raises(ValueError):
        dielectrophoresis_force(a, em, 0.0, g)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("[OK] electrostatics closed forms: capacitance / conduction / polarizability / stress.")
