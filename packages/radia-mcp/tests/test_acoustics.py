"""Closed-form acoustic radiation helpers -- fast regression checks."""

import builtins
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.acoustics import (
    _baffled_piston_resistance_reactance_ratios,
    acoustic_boundary_power_summary,
    acoustic_dtn_from_impedance,
    acoustic_impedance_reflection_summary,
    acoustic_impedance_reflection_sweep_summary,
    acoustic_impedance_radiation_pressure_summary,
    acoustic_impedance_from_dtn,
    baffled_circular_piston_radiation,
    helmholtz_green_3d,
    helmholtz_green_low_frequency_series,
    helmholtz_green_low_frequency_teaching_report,
    low_frequency_helmholtz_kernel_manifest_gate,
    planar_helmholtz_dtn_symbol,
    planar_mode_radiation_impedance,
    pulsating_sphere_radiation,
    spherical_hankel2,
    spherical_helmholtz_dtn_eigenvalue,
    spherical_mode_radiation_impedance,
)


def test_pulsating_sphere_impedance_and_power_conservation():
    a = 0.1
    c = 343.0
    rho = 1.2041
    ka = 1.25
    f = ka * c / (2.0 * math.pi * a)
    v0 = 0.02
    out = pulsating_sphere_radiation(a, f, v0, rho=rho, c=c, sample_radius=12.0 * a)

    denom = 1.0 + ka * ka
    assert out["ka"] == pytest.approx(ka)
    assert out["specific_resistance"] / (rho * c) == pytest.approx(ka * ka / denom)
    assert out["specific_reactance"] / (rho * c) == pytest.approx(ka / denom)

    expected_power = 0.5 * 4.0 * math.pi * a * a * out["specific_resistance"] * v0 * v0
    assert out["radiated_power"] == pytest.approx(expected_power)
    assert out["sample_power"] == pytest.approx(out["radiated_power"], rel=1e-14)


def test_pulsating_sphere_scaling_and_far_pressure():
    a = 0.05
    c = 343.0
    v0 = 0.01
    low1 = pulsating_sphere_radiation(a, 0.05 * c / (2.0 * math.pi * a), v0)
    low2 = pulsating_sphere_radiation(a, 0.10 * c / (2.0 * math.pi * a), v0)
    assert low2["specific_resistance"] / low1["specific_resistance"] == pytest.approx(4.0, rel=0.01)
    assert low2["specific_reactance"] / low1["specific_reactance"] == pytest.approx(2.0, rel=0.01)

    r10 = pulsating_sphere_radiation(a, 500.0, v0, c=c, sample_radius=10.0 * a)
    r20 = pulsating_sphere_radiation(a, 500.0, v0, c=c, sample_radius=20.0 * a)
    assert abs(r20["sample_pressure"]) / abs(r10["sample_pressure"]) == pytest.approx(0.5)

    high = pulsating_sphere_radiation(a, 10.0 * c / (2.0 * math.pi * a), v0)
    assert high["radiation_efficiency"] == pytest.approx(100.0 / 101.0)


def test_pulsating_sphere_validation():
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.0, 100.0, 1.0)
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.1, 0.0, 1.0)
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.1, 100.0, 1.0, sample_radius=0.05)


def test_helmholtz_green_low_frequency_series_terms():
    r = 2.0
    k = 1.0e-3
    out = helmholtz_green_low_frequency_series(r, k, order=4)

    assert out["laplace_term"].real == pytest.approx(1.0 / (4.0 * math.pi * r))
    assert out["laplace_term"].imag == pytest.approx(0.0)
    assert out["terms"][1].real == pytest.approx(0.0)
    assert out["terms"][1].imag == pytest.approx(-k / (4.0 * math.pi))
    assert out["terms"][2].real == pytest.approx(-(k * k) * r / (8.0 * math.pi))
    assert out["terms"][2].imag == pytest.approx(0.0)
    assert out["abs_error"] < 2.0e-17
    assert out["approx"] == pytest.approx(helmholtz_green_3d(r, k))


def test_helmholtz_green_series_convergence_and_validation():
    r = 0.3
    for kr in (1.0e-4, 1.0e-2, 0.1, 0.5):
        k = kr / r
        err2 = helmholtz_green_low_frequency_series(r, k, order=2)["abs_error"]
        err6 = helmholtz_green_low_frequency_series(r, k, order=6)["abs_error"]
        assert err6 < err2
        assert err6 < 1.0e-6

    with pytest.raises(ValueError):
        helmholtz_green_3d(0.0, 1.0)
    with pytest.raises(ValueError):
        helmholtz_green_low_frequency_series(1.0, 1.0, order=-1)


def test_low_frequency_helmholtz_teaching_report_exposes_cancellation_scale():
    r = 0.75
    k = 1.0e-9
    report = helmholtz_green_low_frequency_teaching_report(r, k, order=6)

    assert report["kind"] == "low_frequency_helmholtz_teaching_report"
    assert report["policy"] == "readable_bem_kernel_split_not_production_quadrature"
    assert report["time_convention"] == "exp(+i omega t), outgoing exp(-i k r)"
    assert report["kr_abs"] == pytest.approx(0.75e-9)
    assert report["cancellation_ratio"] > 1.0e8
    assert report["stable_error"] < 2.0e-17
    assert report["correction_agreement"] < 1.0e-15
    assert report["stable_correction"].imag == pytest.approx(-k / (4.0 * math.pi))
    assert report["stable_correction"].real == pytest.approx(-(k * k) * r / (8.0 * math.pi))


def test_low_frequency_helmholtz_kernel_manifest_gate_requires_split_identity():
    report = helmholtz_green_low_frequency_teaching_report(0.75, 1.0e-9, order=6)
    report["kernel_family"] = "helmholtz_single_layer"
    report["low_frequency_strategy"] = "laplace_plus_taylor_regular_part"

    gate = low_frequency_helmholtz_kernel_manifest_gate(
        report,
        expected_kernel_family="helmholtz_single_layer",
        expected_low_frequency_strategy="laplace_plus_taylor_regular_part",
        expected_time_convention="exp(+i omega t), outgoing exp(-i k r)",
        max_kr_abs=1.0e-6,
        min_cancellation_ratio=1.0e8,
    )

    assert gate["status"] == "ok"
    assert gate["checks"]["kernel_family_recorded"] is True
    assert gate["checks"]["low_frequency_strategy_recorded"] is True
    assert gate["checks"]["expected_kernel_family_matches"] is True
    assert gate["checks"]["expected_low_frequency_strategy_matches"] is True
    assert gate["checks"]["expected_time_convention_matches"] is True
    assert gate["checks"]["kr_abs_within_low_frequency_limit"] is True
    assert gate["checks"]["cancellation_ratio_large_enough"] is True

    wrong_strategy = dict(report, low_frequency_strategy="direct_exp_minus_laplace")
    wrong_gate = low_frequency_helmholtz_kernel_manifest_gate(
        wrong_strategy,
        expected_low_frequency_strategy="laplace_plus_taylor_regular_part",
    )
    assert wrong_gate["status"] == "needs_attention"
    assert wrong_gate["checks"]["expected_low_frequency_strategy_matches"] is False

    high_kr = dict(report, kr_abs=0.2)
    high_kr_gate = low_frequency_helmholtz_kernel_manifest_gate(high_kr, max_kr_abs=1.0e-3)
    assert high_kr_gate["status"] == "needs_attention"
    assert high_kr_gate["checks"]["kr_abs_within_low_frequency_limit"] is False


def test_spherical_dtn_monopole_matches_closed_form_and_impedance():
    a = 0.2
    c = 343.0
    rho = 1.2041
    ka = 1.25
    k = ka / a
    f = k * c / (2.0 * math.pi)

    h0 = spherical_hankel2(0, ka)
    assert h0 == pytest.approx(1j * complex(math.cos(ka), -math.sin(ka)) / ka)

    dtn0 = spherical_helmholtz_dtn_eigenvalue(a, k, 0)
    assert dtn0 == pytest.approx(-1.0 / a - 1j * k)

    mode = spherical_mode_radiation_impedance(a, f, 0, rho=rho, c=c)
    sphere = pulsating_sphere_radiation(a, f, 1.0, rho=rho, c=c)
    assert mode["specific_impedance"] == pytest.approx(sphere["specific_impedance"])
    assert mode["radiation_efficiency"] == pytest.approx(sphere["radiation_efficiency"])
    assert mode["reactance_ratio"] == pytest.approx(sphere["reactance_ratio"])


def test_spherical_dtn_matches_radial_finite_difference_for_higher_modes():
    a = 0.4
    k = 5.0
    delta = 1.0e-6 * a
    for degree in range(1, 5):
        h_boundary = spherical_hankel2(degree, k * a)

        def normalized_outgoing(r):
            return spherical_hankel2(degree, k * r) / h_boundary

        finite_difference = (
            normalized_outgoing(a + delta) - normalized_outgoing(a - delta)
        ) / (2.0 * delta)
        dtn = spherical_helmholtz_dtn_eigenvalue(a, k, degree)
        assert dtn == pytest.approx(finite_difference, rel=1.0e-9, abs=1.0e-9)


def test_spherical_mode_radiation_impedance_low_frequency_ordering():
    a = 0.1
    c = 343.0
    low_ka = 0.05
    high_ka = 0.10
    previous_ratio = 0.0
    for degree in range(4):
        low = spherical_mode_radiation_impedance(a, low_ka * c / (2.0 * math.pi * a), degree)
        high = spherical_mode_radiation_impedance(a, high_ka * c / (2.0 * math.pi * a), degree)
        ratio = high["radiation_efficiency"] / low["radiation_efficiency"]
        assert high["radiation_efficiency"] > low["radiation_efficiency"] > 0.0
        assert high["reactance_ratio"] > low["reactance_ratio"] > 0.0
        assert ratio > previous_ratio
        previous_ratio = ratio


def test_spherical_acoustic_helpers_validate_inputs():
    with pytest.raises(ValueError):
        spherical_hankel2(-1, 1.0)
    with pytest.raises(ValueError):
        spherical_hankel2(0, 0.0)
    with pytest.raises(ValueError):
        spherical_helmholtz_dtn_eigenvalue(0.0, 1.0, 0)
    with pytest.raises(ValueError):
        spherical_helmholtz_dtn_eigenvalue(1.0, 0.0, 0)
    with pytest.raises(ValueError):
        spherical_mode_radiation_impedance(1.0, 0.0, 0)


def test_planar_dtn_symbol_matches_normal_oblique_and_evanescent_limits():
    k = 4.0
    normal = planar_helmholtz_dtn_symbol(k, 0.0)
    assert normal["regime"] == "propagating"
    assert normal["normal_wavenumber"] == pytest.approx(k)
    assert normal["dtn_eigenvalue"] == pytest.approx(-1j * k)
    assert abs(normal["symbol_identity_residual"]) < 1.0e-14

    theta = math.radians(60.0)
    kt = k * math.sin(theta)
    oblique = planar_helmholtz_dtn_symbol(k, kt)
    assert oblique["normal_wavenumber"] == pytest.approx(k * math.cos(theta))
    assert oblique["dtn_eigenvalue"] == pytest.approx(-1j * k * math.cos(theta))
    assert abs(oblique["symbol_identity_residual"]) < 1.0e-14

    evanescent = planar_helmholtz_dtn_symbol(k, 2.0 * k)
    assert evanescent["regime"] == "evanescent"
    assert evanescent["normal_wavenumber"].real == pytest.approx(0.0)
    assert evanescent["normal_wavenumber"].imag == pytest.approx(-math.sqrt(3.0) * k)
    assert evanescent["dtn_eigenvalue"].real == pytest.approx(-math.sqrt(3.0) * k)
    assert evanescent["dtn_eigenvalue"].imag == pytest.approx(0.0)
    assert abs(evanescent["symbol_identity_residual"]) < 1.0e-14


def test_planar_mode_radiation_impedance_angle_and_evanescent_modes():
    rho = 1.2041
    c = 343.0
    f = 1000.0
    k = 2.0 * math.pi * f / c

    normal = planar_mode_radiation_impedance(f, incidence_angle_rad=0.0, rho=rho, c=c)
    assert normal["specific_impedance"] == pytest.approx(rho * c)
    assert normal["normalized_impedance"] == pytest.approx(1.0)
    assert normal["radiation_efficiency"] == pytest.approx(1.0)

    theta = math.radians(60.0)
    oblique = planar_mode_radiation_impedance(f, incidence_angle_rad=theta, rho=rho, c=c)
    assert oblique["normalized_impedance"] == pytest.approx(1.0 / math.cos(theta))
    assert oblique["dtn_eigenvalue"] == pytest.approx(-1j * k * math.cos(theta))

    evanescent = planar_mode_radiation_impedance(f, tangential_wavenumber=1.5 * k, rho=rho, c=c)
    assert evanescent["regime"] == "evanescent"
    assert evanescent["radiation_efficiency"] == pytest.approx(0.0)
    assert evanescent["reactance_ratio"] == pytest.approx(1.0 / math.sqrt(1.5 * 1.5 - 1.0))
    assert evanescent["specific_impedance"].imag > 0.0


def test_planar_acoustic_helpers_validate_inputs():
    with pytest.raises(ValueError):
        planar_helmholtz_dtn_symbol(0.0, 0.0)
    with pytest.raises(ValueError):
        planar_helmholtz_dtn_symbol(1.0, -1.0)
    with pytest.raises(ValueError):
        planar_mode_radiation_impedance(0.0, incidence_angle_rad=0.0)
    with pytest.raises(ValueError):
        planar_mode_radiation_impedance(100.0)
    with pytest.raises(ValueError):
        planar_mode_radiation_impedance(100.0, tangential_wavenumber=1.0, incidence_angle_rad=0.0)
    with pytest.raises(ValueError):
        planar_mode_radiation_impedance(100.0, incidence_angle_rad=0.5 * math.pi)


def test_acoustic_impedance_dtn_conversion_matches_planar_and_spherical_modes():
    rho, c, f = 1.2, 340.0, 1000.0
    omega = 2.0 * math.pi * f
    k = omega / c

    normal = acoustic_dtn_from_impedance(f, specific_impedance=rho * c, rho=rho)
    assert normal["dtn_eigenvalue"] == pytest.approx(-1j * k)
    roundtrip_normal = acoustic_impedance_from_dtn(f, normal["dtn_eigenvalue"], rho=rho)
    assert roundtrip_normal["specific_impedance"] == pytest.approx(rho * c)

    oblique = planar_mode_radiation_impedance(f, incidence_angle_rad=math.radians(50.0), rho=rho, c=c)
    from_imp = acoustic_dtn_from_impedance(f, specific_impedance=oblique["specific_impedance"], rho=rho)
    from_adm = acoustic_dtn_from_impedance(f, specific_admittance=1.0 / oblique["specific_impedance"], rho=rho)
    assert from_imp["dtn_eigenvalue"] == pytest.approx(oblique["dtn_eigenvalue"])
    assert from_adm["dtn_eigenvalue"] == pytest.approx(oblique["dtn_eigenvalue"])

    sphere = spherical_mode_radiation_impedance(0.17, f, 2, rho=rho, c=c)
    roundtrip = acoustic_impedance_from_dtn(f, sphere["dtn_eigenvalue"], rho=rho)
    assert roundtrip["specific_impedance"] == pytest.approx(sphere["specific_impedance"])


def test_acoustic_impedance_dtn_conversion_validation():
    with pytest.raises(ValueError):
        acoustic_dtn_from_impedance(0.0, specific_impedance=1.0)
    with pytest.raises(ValueError):
        acoustic_dtn_from_impedance(100.0, specific_impedance=0.0)
    with pytest.raises(ValueError):
        acoustic_dtn_from_impedance(100.0)
    with pytest.raises(ValueError):
        acoustic_dtn_from_impedance(100.0, specific_impedance=1.0, specific_admittance=1.0)
    with pytest.raises(ValueError):
        acoustic_impedance_from_dtn(100.0, 0.0)
    with pytest.raises(ValueError):
        acoustic_impedance_from_dtn(100.0, 1.0, rho=0.0)


def test_acoustic_boundary_power_summary_peak_and_rms_conventions():
    rho, c = 1.2, 340.0
    z0 = rho * c
    pressure = 2.0
    velocity = pressure / z0
    area = 0.25

    peak = acoustic_boundary_power_summary(pressure, velocity, area=area, amplitude="peak")
    rms = acoustic_boundary_power_summary(pressure, velocity, area=area, amplitude="rms")

    assert peak["specific_impedance"] == pytest.approx(z0)
    assert peak["active_intensity"] == pytest.approx(pressure * pressure / (2.0 * z0))
    assert peak["reactive_intensity"] == pytest.approx(0.0)
    assert peak["active_power"] == pytest.approx(area * pressure * pressure / (2.0 * z0))
    assert rms["active_power"] == pytest.approx(2.0 * peak["active_power"])
    assert rms["phasor_average_factor"] == pytest.approx(1.0)


def test_acoustic_boundary_power_summary_reactive_trace_and_validation():
    z_reactive = 1.0j * 120.0
    velocity = 0.03 - 0.01j
    pressure = z_reactive * velocity
    area = 0.4
    out = acoustic_boundary_power_summary(pressure, velocity, area=area)

    assert out["active_intensity"] == pytest.approx(0.0)
    assert out["reactive_intensity"] == pytest.approx(0.5 * 120.0 * abs(velocity) ** 2)
    assert out["reactive_power"] == pytest.approx(area * out["reactive_intensity"])
    assert out["specific_impedance"] == pytest.approx(z_reactive)

    with pytest.raises(ValueError):
        acoustic_boundary_power_summary(1.0, 1.0, area=-1.0)
    with pytest.raises(ValueError):
        acoustic_boundary_power_summary(1.0, 1.0, amplitude="phasor")
    with pytest.raises(ValueError):
        acoustic_boundary_power_summary(complex(float("nan"), 0.0), 1.0)


def test_acoustic_impedance_reflection_matched_and_mismatched_loads():
    rho, c = 1.2041, 343.0
    z0 = rho * c
    matched = acoustic_impedance_reflection_summary(z0, incident_pressure=2.0, rho=rho, c=c)
    assert matched["pressure_reflection_coefficient"] == pytest.approx(0.0)
    assert matched["absorption_coefficient"] == pytest.approx(1.0)
    assert matched["incident_intensity"] == pytest.approx(2.0 * 2.0 / (2.0 * z0))
    assert matched["absorbed_intensity"] == pytest.approx(matched["incident_intensity"])
    assert abs(matched["power_balance_residual"]) < 1.0e-15

    twice = acoustic_impedance_reflection_summary(2.0 * z0, rho=rho, c=c)
    assert twice["pressure_reflection_coefficient"] == pytest.approx(1.0 / 3.0)
    assert twice["power_reflection_coefficient"] == pytest.approx(1.0 / 9.0)
    assert twice["absorption_coefficient"] == pytest.approx(8.0 / 9.0)
    assert twice["boundary_active_intensity_into_load"] == pytest.approx(twice["absorbed_intensity"])


def test_acoustic_impedance_reflection_reactive_and_oblique_limits():
    rho, c = 1.2, 340.0
    z0 = rho * c
    reactive = acoustic_impedance_reflection_summary(1j * z0, rho=rho, c=c)
    assert abs(reactive["pressure_reflection_coefficient"]) == pytest.approx(1.0)
    assert reactive["absorption_coefficient"] == pytest.approx(0.0)
    assert reactive["absorbed_intensity"] == pytest.approx(0.0)
    assert reactive["boundary_reactive_intensity_into_load"] != pytest.approx(0.0)

    pressure_release = acoustic_impedance_reflection_summary(0.0, rho=rho, c=c)
    assert pressure_release["pressure_reflection_coefficient"] == pytest.approx(-1.0)
    assert pressure_release["total_boundary_pressure"] == pytest.approx(0.0)
    assert pressure_release["absorption_coefficient"] == pytest.approx(0.0)

    theta = math.radians(60.0)
    z_normal = z0 / math.cos(theta)
    oblique_matched = acoustic_impedance_reflection_summary(z_normal, incidence_angle_rad=theta, rho=rho, c=c)
    assert oblique_matched["characteristic_normal_impedance"] == pytest.approx(z_normal)
    assert oblique_matched["pressure_reflection_coefficient"] == pytest.approx(0.0)
    assert oblique_matched["absorption_coefficient"] == pytest.approx(1.0)

    with pytest.raises(ValueError):
        acoustic_impedance_reflection_summary(z0, incidence_angle_rad=0.5 * math.pi)
    with pytest.raises(ValueError):
        acoustic_impedance_reflection_summary(complex(float("inf"), 0.0))
    with pytest.raises(ValueError):
        acoustic_impedance_reflection_summary(z0, amplitude="complex")


def test_acoustic_impedance_slot_power_identity_for_passive_loads():
    z0 = 1.0
    cases = [
        0.0 + 0.25j,
        0.0 + 1.0j,
        0.1 + 0.0j,
        0.1 + 0.5j,
        1.0 + 0.0j,
        2.0 + 0.5j,
        10.0 + 0.0j,
    ]
    max_abs_residual = 0.0
    max_rel_residual = 0.0
    reactive_absorption = []

    for z_norm in cases:
        out = acoustic_impedance_reflection_summary(z_norm, rho=1.0, c=1.0)
        reflection = complex(out["pressure_reflection_coefficient"])
        p_total = 1.0 + reflection
        v_total = 1.0 - reflection
        boundary_power = 0.5 * (p_total * v_total.conjugate()).real
        expected_power = 0.5 * (1.0 - abs(reflection) ** 2)
        abs_residual = abs(boundary_power - expected_power)
        max_abs_residual = max(max_abs_residual, abs_residual)
        if abs(expected_power) > 1.0e-14:
            max_rel_residual = max(max_rel_residual, abs_residual / abs(expected_power))
        if abs(z_norm.real) < 1.0e-14:
            reactive_absorption.append(abs(out["absorption_coefficient"]))

        assert out["absorption_coefficient"] >= -1.0e-14
        assert out["boundary_active_intensity_into_load"] == pytest.approx(expected_power)

    assert max_abs_residual < 1.0e-12
    assert max_rel_residual < 1.0e-12
    assert max(reactive_absorption) < 1.0e-12


def test_acoustic_impedance_negative_real_part_is_active_not_absorbing():
    out = acoustic_impedance_reflection_summary(-2.0, rho=1.0, c=1.0)
    reflection = complex(out["pressure_reflection_coefficient"])
    expected_power = 0.5 * (1.0 - abs(reflection) ** 2)

    assert abs(reflection) > 1.0
    assert out["absorption_coefficient"] < 0.0
    assert out["boundary_active_intensity_into_load"] == pytest.approx(expected_power)
    assert expected_power < 0.0


def test_acoustic_impedance_radiation_pressure_absorber_reflector_limits():
    rho, c = 1.2041, 343.0
    z0 = rho * c
    incident_pressure = 2.0
    area = 0.5
    incident_intensity = incident_pressure * incident_pressure / (2.0 * z0)

    matched = acoustic_impedance_radiation_pressure_summary(
        z0,
        area=area,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
    )
    twice = acoustic_impedance_radiation_pressure_summary(
        2.0 * z0,
        area=area,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
    )
    reactive = acoustic_impedance_radiation_pressure_summary(
        1j * z0,
        area=area,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
    )

    assert matched["power_reflection_coefficient"] == pytest.approx(0.0)
    assert matched["normal_momentum_pressure_Pa"] == pytest.approx(incident_intensity / c)
    assert matched["normal_force_N"] == pytest.approx(area * incident_intensity / c)
    assert abs(matched["force_balance_residual_N"]) < 1.0e-18

    assert twice["power_reflection_coefficient"] == pytest.approx(1.0 / 9.0)
    assert twice["absorption_coefficient"] == pytest.approx(8.0 / 9.0)
    assert twice["normal_momentum_pressure_Pa"] == pytest.approx((10.0 / 9.0) * incident_intensity / c)
    assert twice["normal_force_N"] == pytest.approx(twice["force_from_absorptance_reflectance_N"])

    assert reactive["power_reflection_coefficient"] == pytest.approx(1.0)
    assert reactive["absorption_coefficient"] == pytest.approx(0.0)
    assert reactive["normal_momentum_pressure_Pa"] == pytest.approx(2.0 * incident_intensity / c)

    with pytest.raises(ValueError):
        acoustic_impedance_radiation_pressure_summary(z0, area=-1.0)


def test_acoustic_impedance_reflection_sweep_tracks_absorption_force_and_passivity():
    rho = 1.2
    c = 340.0
    z0 = rho * c
    area = 0.25
    incident_pressure = 2.0
    frequencies = [100.0, 200.0, 300.0]
    sweep = acoustic_impedance_reflection_sweep_summary(
        frequencies,
        [z0, 2.0 * z0, 1j * z0],
        area=area,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
    )

    incident_intensity = 0.5 * incident_pressure * incident_pressure / z0
    assert sweep["n_points"] == 3
    assert sweep["frequency_monotonic_increasing"] is True
    assert sweep["status"] == "ok"
    assert sweep["max_absorption_frequency_Hz"] == pytest.approx(100.0)
    assert sweep["min_absorption_frequency_Hz"] == pytest.approx(300.0)
    assert sweep["max_force_frequency_Hz"] == pytest.approx(300.0)
    assert sweep["min_force_frequency_Hz"] == pytest.approx(100.0)
    assert sweep["rows"][0]["absorption_coefficient"] == pytest.approx(1.0)
    assert sweep["rows"][1]["power_reflection_coefficient"] == pytest.approx(1.0 / 9.0)
    assert sweep["rows"][1]["absorption_coefficient"] == pytest.approx(8.0 / 9.0)
    assert sweep["rows"][2]["power_reflection_coefficient"] == pytest.approx(1.0)
    assert sweep["max_normal_force_N"] == pytest.approx(area * 2.0 * incident_intensity / c)

    active = acoustic_impedance_reflection_sweep_summary(
        [400.0],
        [-2.0 * z0],
        rho=rho,
        c=c,
        passivity_tolerance=1.0e-6,
    )
    assert active["status"] == "needs_attention"
    assert active["passivity_violation_count"] == 1
    assert active["max_passivity_excess_absorption"] == pytest.approx(8.0)
    assert active["passivity_violation_rows"][0]["absorption_coefficient"] == pytest.approx(-8.0)


def test_baffled_circular_piston_impedance_scaling_and_power():
    pytest.importorskip("scipy")
    a = 0.08
    c = 343.0
    rho = 1.2041
    v0 = 0.02
    low1 = baffled_circular_piston_radiation(a, 0.025 * c / (2.0 * math.pi * a), v0, rho=rho, c=c)
    low2 = baffled_circular_piston_radiation(a, 0.050 * c / (2.0 * math.pi * a), v0, rho=rho, c=c)
    mid = baffled_circular_piston_radiation(a, 1.0 * c / (2.0 * math.pi * a), v0, rho=rho, c=c)
    high = baffled_circular_piston_radiation(a, 20.0 * c / (2.0 * math.pi * a), v0, rho=rho, c=c)

    assert low2["radiation_efficiency"] / low1["radiation_efficiency"] == pytest.approx(4.0, rel=1.0e-3)
    assert low2["reactance_ratio"] / low1["reactance_ratio"] == pytest.approx(2.0, rel=1.0e-3)
    assert low1["radiation_efficiency"] == pytest.approx(low1["low_ka_resistance_asymptote"], rel=5.0e-4)
    assert low1["reactance_ratio"] == pytest.approx(low1["low_ka_reactance_asymptote"], rel=5.0e-4)
    assert high["radiation_efficiency"] == pytest.approx(1.0, abs=0.05)
    assert mid["radiated_power"] == pytest.approx(
        0.5 * mid["surface_area"] * mid["specific_resistance"] * v0 * v0
    )
    assert mid["volume_velocity_impedance"] == pytest.approx(
        mid["specific_impedance"] / mid["surface_area"]
    )

    with pytest.raises(ValueError):
        baffled_circular_piston_radiation(0.0, 100.0)
    with pytest.raises(ValueError):
        baffled_circular_piston_radiation(a, 0.0)


def test_baffled_circular_piston_fallback_without_scipy():
    low1 = _baffled_piston_resistance_reactance_ratios(0.025, prefer_scipy=False)
    low2 = _baffled_piston_resistance_reactance_ratios(0.050, prefer_scipy=False)
    high = _baffled_piston_resistance_reactance_ratios(20.0, prefer_scipy=False)

    assert low1[2] == "fallback"
    assert low2[0] / low1[0] == pytest.approx(4.0, rel=1.0e-3)
    assert low2[1] / low1[1] == pytest.approx(2.0, rel=1.0e-3)
    assert low1[0] == pytest.approx(0.5 * 0.025 * 0.025, rel=5.0e-4)
    assert low1[1] == pytest.approx(8.0 * 0.025 / (3.0 * math.pi), rel=5.0e-4)
    assert high[0] == pytest.approx(1.0, abs=0.05)

    with pytest.raises(ValueError):
        _baffled_piston_resistance_reactance_ratios(0.0, prefer_scipy=False)


def test_baffled_circular_piston_public_fallback_when_scipy_missing(monkeypatch):
    original_import = builtins.__import__

    def _blocked_import(name, *args, **kwargs):
        if name == "scipy" or name.startswith("scipy."):
            raise ModuleNotFoundError("No module named 'scipy'")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _blocked_import)

    out = baffled_circular_piston_radiation(0.08, 0.025 * 343.0 / (2.0 * math.pi * 0.08))
    assert out["special_function_source"] == "fallback"
    assert out["radiation_efficiency"] == pytest.approx(out["low_ka_resistance_asymptote"], rel=5.0e-4)
