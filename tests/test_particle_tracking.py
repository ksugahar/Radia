from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy.integrate")

from radia.particle_tracking import (
    FIVE_MOMENTUM_OFFSETS,
    ParticleSpecies,
    SPEED_OF_LIGHT_M_S,
    TrackingBox,
    TrackingPlane,
    speed_from_kinetic_voltage,
    track_lorentz_ivp,
    track_two_momentum_exit_dispersion,
    fit_five_momentum_exit_optics,
    track_five_momentum_exit_optics,
    uniform_magnetic_trajectory,
    velocity_from_kinetic_voltage,
)


ELECTRON = ParticleSpecies(-1.60217733e-19, 9.1093897e-31)


def _zero(x, y, z):
    return 0.0, 0.0, 0.0


def _uniform_bz(x, y, z):
    return 0.0, 0.0, 1.0e-3


def test_uniform_magnetic_orbit_matches_closed_form_and_preserves_energy():
    times = np.linspace(0.0, 5.0e-8, 101)
    velocity0 = velocity_from_kinetic_voltage(
        ELECTRON, 100.0, (1.0, 0.0, 0.0), relativistic=False
    )
    result = track_lorentz_ivp(
        ELECTRON,
        (0.0, 0.0, 0.0),
        velocity0,
        times,
        magnetic_flux_density_t=_uniform_bz,
        relativistic=False,
    )
    exact = uniform_magnetic_trajectory(
        ELECTRON, 100.0, times, magnetic_flux_density_t=1.0e-3
    )
    np.testing.assert_allclose(result["position_m"], exact["position_m"], atol=2.0e-11)
    np.testing.assert_allclose(result["velocity_m_s"], exact["velocity_m_s"], atol=4.0e-3)
    assert result["maximum_relative_kinetic_energy_drift"] < 2.0e-10


def test_charge_sign_reverses_magnetic_deflection():
    positron = ParticleSpecies(abs(ELECTRON.charge_c), ELECTRON.mass_kg)
    times = [0.0, 5.0e-9]
    electron = uniform_magnetic_trajectory(
        ELECTRON, 10.0, times, magnetic_flux_density_t=1.0e-3
    )
    positive = uniform_magnetic_trajectory(
        positron, 10.0, times, magnetic_flux_density_t=1.0e-3
    )
    assert electron["position_m"][-1, 1] > 0.0
    assert positive["position_m"][-1, 1] < 0.0


def test_relativistic_voltage_conversion_stays_below_c():
    speed = speed_from_kinetic_voltage(ELECTRON, 5.0e6, relativistic=True)
    classical = speed_from_kinetic_voltage(ELECTRON, 5.0e6, relativistic=False)
    assert speed < SPEED_OF_LIGHT_M_S
    assert classical > SPEED_OF_LIGHT_M_S


def test_electric_field_changes_kinetic_energy():
    def electric_x(x, y, z):
        return -1000.0, 0.0, 0.0

    result = track_lorentz_ivp(
        ELECTRON,
        (0.0, 0.0, 0.0),
        (1.0e5, 0.0, 0.0),
        np.linspace(0.0, 1.0e-8, 21),
        electric_field_v_m=electric_x,
    )
    assert result["kinetic_energy_j"][-1] > result["kinetic_energy_j"][0]


def test_first_stop_box_face_is_reported():
    result = track_lorentz_ivp(
        ELECTRON,
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        np.linspace(0.0, 1.0, 21),
        electric_field_v_m=_zero,
        stop_box=TrackingBox((-1.0, -1.0, -1.0), (0.5, 1.0, 1.0)),
    )
    assert result["stop_event"]["face"] == "x_maximum"
    assert result["stop_event"]["time_s"] == pytest.approx(0.25)
    assert result["stop_event"]["position_m"][0] == pytest.approx(0.5)


def test_exit_plane_reports_exact_position_and_velocity():
    plane = TrackingPlane((0.25, 0.0, 0.0), (2.0, 0.0, 0.0), direction=1)
    result = track_lorentz_ivp(
        ELECTRON,
        (0.0, 0.1, 0.0),
        (2.0, 0.0, 0.0),
        np.linspace(0.0, 1.0, 11),
        magnetic_flux_density_t=_zero,
        stop_plane=plane,
    )
    assert result["stop_event"]["face"] == "plane"
    assert result["stop_event"]["time_s"] == pytest.approx(0.125)
    np.testing.assert_allclose(
        result["stop_event"]["position_m"], (0.25, 0.1, 0.0), atol=1e-14
    )
    np.testing.assert_allclose(
        result["stop_event"]["velocity_m_s"], (2.0, 0.0, 0.0), atol=1e-14
    )


def test_two_momenta_have_zero_exit_dispersion_in_zero_field():
    result = track_two_momentum_exit_dispersion(
        ELECTRON,
        (0.0, 0.1, 0.0),
        (2.0, 0.0, 0.0),
        np.linspace(0.0, 1.0, 21),
        TrackingPlane((0.5, 0.0, 0.0), (1.0, 0.0, 0.0), direction=1),
        relative_momentum_offset=1.0e-3,
        transverse_direction=(0.0, 1.0, 0.0),
        magnetic_flux_density_t=_zero,
    )
    assert result["eta_m"] == pytest.approx(0.0, abs=1e-12)
    assert result["coincident_exit_error_m"] == pytest.approx(0.0, abs=1e-12)


def test_two_momentum_exit_dispersion_detects_magnetic_separation():
    result = track_two_momentum_exit_dispersion(
        ELECTRON,
        (0.0, 0.0, 0.0),
        (2.0e6, 0.0, 0.0),
        np.linspace(0.0, 2.0e-7, 101),
        TrackingPlane((0.005, 0.0, 0.0), (1.0, 0.0, 0.0), direction=1),
        relative_momentum_offset=1.0e-3,
        transverse_direction=(0.0, 1.0, 0.0),
        magnetic_flux_density_t=_uniform_bz,
    )
    assert abs(result["eta_m"]) > 1.0e-4
    assert result["coincident_exit_error_m"] > 1.0e-7


def test_five_momentum_fit_recovers_quadratic_exit_optics():
    offsets = np.asarray(FIVE_MOMENTUM_OFFSETS)
    positions = 2.0e-4 + 3.0e-3 * offsets + 4.0 * offsets**2
    angles = -1.0e-4 - 5.0e-4 * offsets + 0.2 * offsets**2

    result = fit_five_momentum_exit_optics(offsets, positions, angles)

    np.testing.assert_allclose(
        result["linear_regression_weights"],
        (-400.0, -200.0, 0.0, 200.0, 400.0),
        atol=1.0e-12,
    )
    assert result["x0_m"] == pytest.approx(2.0e-4)
    assert result["psi0_rad"] == pytest.approx(-1.0e-4)
    assert result["eta_m"] == pytest.approx(3.0e-3)
    assert result["eta_prime_rad"] == pytest.approx(-5.0e-4)
    assert result["x_quadratic_m"] == pytest.approx(4.0)
    assert result["psi_quadratic_rad"] == pytest.approx(0.2)
    assert result["max_x_residual_m"] < 1.0e-14
    assert result["max_psi_residual_rad"] < 1.0e-14
    assert result["pass_all"]


def test_five_momenta_have_zero_exit_optics_in_zero_field():
    result = track_five_momentum_exit_optics(
        ELECTRON,
        (0.0, 0.1, 0.0),
        (2.0, 0.0, 0.0),
        np.linspace(0.0, 1.0, 21),
        TrackingPlane((0.5, 0.0, 0.0), (1.0, 0.0, 0.0), direction=1),
        reference_exit_point_m=(0.5, 0.1, 0.0),
        transverse_direction=(0.0, 1.0, 0.0),
        longitudinal_direction=(1.0, 0.0, 0.0),
        magnetic_flux_density_t=_zero,
    )

    assert result["x0_m"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["psi0_rad"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["eta_m"] == pytest.approx(0.0, abs=1.0e-12)
    assert result["eta_prime_rad"] == pytest.approx(0.0, abs=1.0e-12)
    assert len(result["tracks"]) == 5
    assert result["pass_all"]


def test_five_momentum_fit_rejects_a_noncanonical_sample_count():
    with pytest.raises(ValueError, match="exactly five"):
        fit_five_momentum_exit_optics(
            (-1.0e-3, 0.0, 1.0e-3),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        )


@pytest.mark.parametrize(
    "offsets, message",
    [
        ((-1.0e-3, 0.0, -5.0e-4, 5.0e-4, 1.0e-3), "strictly increasing"),
        ((-1.0e-3, -5.0e-4, 2.0e-4, 5.0e-4, 1.0e-3), "exactly one zero"),
    ],
)
def test_five_momentum_track_rejects_invalid_offsets_before_tracking(
        offsets, message):
    with pytest.raises(ValueError, match=message):
        track_five_momentum_exit_optics(
            ELECTRON,
            (0.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            np.linspace(0.0, 1.0, 21),
            TrackingPlane((0.5, 0.0, 0.0), (1.0, 0.0, 0.0), direction=1),
            reference_exit_point_m=(0.5, 0.0, 0.0),
            transverse_direction=(0.0, 1.0, 0.0),
            relative_momentum_offsets=offsets,
            magnetic_flux_density_t=_zero,
        )


@pytest.mark.parametrize(
    "factory",
    [
        lambda: ParticleSpecies(0.0, 1.0),
        lambda: ParticleSpecies(1.0, 0.0),
        lambda: TrackingBox((0.0, 0.0, 0.0), (0.0, 1.0, 1.0)),
        lambda: TrackingPlane((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
    ],
)
def test_invalid_identity_is_rejected(factory):
    with pytest.raises(ValueError):
        factory()
