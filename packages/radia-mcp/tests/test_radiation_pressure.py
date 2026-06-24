"""RF radiation-pressure helpers from Poynting flux."""

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    ETA0,
    oblique_radiation_pressure_summary,
    plane_wave_intensity_from_electric_field,
    poynting_patch_force_summary,
    radiation_force_from_normal_scattering,
    radiation_force_from_power,
    radiation_pressure_from_intensity,
    radiation_pressure_summary,
    radiation_scattering_force_summary,
    two_port_scattering_momentum_force_summary,
)


def test_plane_wave_intensity_accepts_rms_and_peak_conventions():
    e_rms = math.sqrt(ETA0)
    e_peak = math.sqrt(2.0 * ETA0)

    assert plane_wave_intensity_from_electric_field(e_rms, amplitude="rms") == pytest.approx(1.0)
    assert plane_wave_intensity_from_electric_field(e_peak, amplitude="peak") == pytest.approx(1.0)


def test_radiation_pressure_absorber_and_reflector_limits():
    intensity = 3.0

    absorber = radiation_pressure_from_intensity(intensity, absorptance=1.0, reflectance=0.0)
    reflector = radiation_pressure_from_intensity(intensity, absorptance=0.0, reflectance=1.0)
    half_absorbing_mirror = radiation_pressure_from_intensity(
        intensity, absorptance=0.25, reflectance=0.5
    )

    assert absorber == pytest.approx(intensity / C0)
    assert reflector == pytest.approx(2.0 * intensity / C0)
    assert half_absorbing_mirror == pytest.approx(1.25 * intensity / C0)


def test_radiation_force_from_power_is_area_integrated_pressure():
    power = 5.0
    area = 0.2
    intensity = power / area
    summary = radiation_pressure_summary(
        intensity,
        area_m2=area,
        absorptance=0.0,
        reflectance=1.0,
    )

    assert radiation_force_from_power(power, absorptance=0.0, reflectance=1.0) == pytest.approx(
        2.0 * power / C0
    )
    assert summary["incident_power_W"] == pytest.approx(power)
    assert summary["pressure_Pa"] == pytest.approx(2.0 * intensity / C0)
    assert summary["force_N"] == pytest.approx(2.0 * power / C0)
    assert summary["momentum_transfer_factor"] == pytest.approx(2.0)


def test_radiation_force_from_normal_scattering_matches_momentum_balance():
    power = 5.0

    absorber = radiation_scattering_force_summary(power, reflectance=0.0, transmittance=0.0)
    reflector = radiation_scattering_force_summary(power, reflectance=1.0, transmittance=0.0)
    transparent = radiation_scattering_force_summary(power, reflectance=0.0, transmittance=1.0)
    lossy_partial = radiation_scattering_force_summary(power, reflectance=0.5, transmittance=0.25)

    assert radiation_force_from_normal_scattering(power, 0.0, 0.0) == pytest.approx(power / C0)
    assert absorber["force_N"] == pytest.approx(power / C0)
    assert reflector["force_N"] == pytest.approx(2.0 * power / C0)
    assert transparent["force_N"] == pytest.approx(0.0)
    assert lossy_partial["absorptance"] == pytest.approx(0.25)
    assert lossy_partial["momentum_transfer_factor"] == pytest.approx(1.25)
    assert lossy_partial["force_N"] == pytest.approx(
        radiation_force_from_power(power, absorptance=0.25, reflectance=0.5)
    )
    assert lossy_partial["force_from_absorptance_reflectance_N"] == pytest.approx(
        lossy_partial["force_N"]
    )


def test_two_port_scattering_momentum_force_vector():
    power = 3.0
    kin = (1.0, 0.0, 0.0)

    straight = two_port_scattering_momentum_force_summary(
        power,
        kin,
        reflectance=0.0,
        transmittance=1.0,
        transmitted_direction=kin,
    )
    mirror = two_port_scattering_momentum_force_summary(
        power,
        kin,
        reflectance=1.0,
        transmittance=0.0,
    )
    bend = two_port_scattering_momentum_force_summary(
        power,
        kin,
        reflectance=0.0,
        transmittance=1.0,
        transmitted_direction=(0.0, 1.0, 0.0),
    )
    lossy_straight = two_port_scattering_momentum_force_summary(
        power,
        kin,
        reflectance=0.1,
        transmittance=0.6,
        transmitted_direction=kin,
    )

    assert straight["force_N"] == pytest.approx([0.0, 0.0, 0.0])
    assert mirror["force_N"] == pytest.approx([2.0 * power / C0, 0.0, 0.0])
    assert bend["force_N"] == pytest.approx([power / C0, -power / C0, 0.0])
    assert bend["force_magnitude_N"] == pytest.approx(math.sqrt(2.0) * power / C0)
    assert lossy_straight["force_N"] == pytest.approx([0.5 * power / C0, 0.0, 0.0])
    assert lossy_straight["axial_force_along_incident_direction_N"] == pytest.approx(
        lossy_straight["straight_through_equivalent_force_N"]
    )


def test_oblique_radiation_pressure_reduces_to_normal_incidence():
    intensity = 7.0
    area = 0.25
    normal = radiation_pressure_summary(intensity, area_m2=area)
    oblique = oblique_radiation_pressure_summary(intensity, 0.0, area_m2=area)

    assert oblique["normal_pressure_Pa"] == pytest.approx(normal["pressure_Pa"])
    assert oblique["normal_force_N"] == pytest.approx(normal["force_N"])
    assert oblique["tangential_force_N"] == pytest.approx(0.0)
    assert oblique["incident_power_on_patch_W"] == pytest.approx(intensity * area)


def test_oblique_radiation_pressure_splits_absorbed_tangential_momentum():
    intensity = 9.0
    area = 0.5
    angle = math.radians(60.0)
    c = math.cos(angle)
    s = math.sin(angle)

    absorber = oblique_radiation_pressure_summary(
        intensity, angle, area_m2=area, absorptance=1.0, reflectance=0.0
    )
    reflector = oblique_radiation_pressure_summary(
        intensity, angle, area_m2=area, absorptance=0.0, reflectance=1.0
    )

    assert absorber["normal_force_N"] == pytest.approx(intensity * area * c * c / C0)
    assert absorber["tangential_force_N"] == pytest.approx(intensity * area * s * c / C0)
    assert reflector["normal_force_N"] == pytest.approx(2.0 * intensity * area * c * c / C0)
    assert reflector["tangential_force_N"] == pytest.approx(0.0)
    assert reflector["incident_power_on_patch_W"] == pytest.approx(intensity * area * c)


def test_poynting_patch_force_vector_matches_oblique_summary():
    intensity = 9.0
    area = 0.5
    angle = math.radians(60.0)
    c = math.cos(angle)
    s = math.sin(angle)
    poynting = (intensity * s, 0.0, -intensity * c)
    normal = (0.0, 0.0, 1.0)
    scalar_absorber = oblique_radiation_pressure_summary(
        intensity, angle, area_m2=area, absorptance=1.0, reflectance=0.0
    )
    vector_absorber = poynting_patch_force_summary(
        poynting, normal, area_m2=area, absorptance=1.0, reflectance=0.0
    )
    scalar_reflector = oblique_radiation_pressure_summary(
        intensity, angle, area_m2=area, absorptance=0.0, reflectance=1.0
    )
    vector_reflector = poynting_patch_force_summary(
        poynting, normal, area_m2=area, absorptance=0.0, reflectance=1.0
    )

    assert vector_absorber["incident_power_on_patch_W"] == pytest.approx(intensity * area * c)
    assert vector_absorber["normal_force_into_surface_N"] == pytest.approx(
        scalar_absorber["normal_force_N"]
    )
    assert vector_absorber["tangential_force_magnitude_N"] == pytest.approx(
        scalar_absorber["tangential_force_N"]
    )
    assert vector_reflector["normal_force_into_surface_N"] == pytest.approx(
        scalar_reflector["normal_force_N"]
    )
    assert vector_reflector["tangential_force_magnitude_N"] == pytest.approx(0.0)


def test_radiation_pressure_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        plane_wave_intensity_from_electric_field(1.0, impedance_ohm=0.0)
    with pytest.raises(ValueError):
        plane_wave_intensity_from_electric_field(1.0, amplitude="phasor")
    with pytest.raises(ValueError):
        radiation_pressure_from_intensity(-1.0)
    with pytest.raises(ValueError):
        radiation_force_from_power(-1.0)
    with pytest.raises(ValueError):
        radiation_force_from_normal_scattering(-1.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        radiation_force_from_normal_scattering(1.0, -0.1, 0.0)
    with pytest.raises(ValueError):
        radiation_force_from_normal_scattering(1.0, 0.5, 0.6)
    with pytest.raises(ValueError):
        radiation_scattering_force_summary(1.0, 0.0, -0.1)
    with pytest.raises(ValueError):
        two_port_scattering_momentum_force_summary(1.0, (0.0, 0.0, 0.0), 0.0, 0.0)
    with pytest.raises(ValueError):
        two_port_scattering_momentum_force_summary(
            1.0,
            (1.0, 0.0, 0.0),
            0.0,
            0.5,
            transmitted_direction=(1.0, 0.0),
        )
    with pytest.raises(ValueError):
        radiation_pressure_summary(1.0, area_m2=-1.0)
    with pytest.raises(ValueError):
        radiation_pressure_from_intensity(1.0, absorptance=0.6, reflectance=0.5)
    with pytest.raises(ValueError):
        oblique_radiation_pressure_summary(1.0, -0.1)
    with pytest.raises(ValueError):
        oblique_radiation_pressure_summary(1.0, math.pi)
    with pytest.raises(ValueError):
        oblique_radiation_pressure_summary(1.0, 0.0, area_m2=-1.0)
    with pytest.raises(ValueError):
        poynting_patch_force_summary((1.0, 0.0), (0.0, 0.0, 1.0))
    with pytest.raises(ValueError):
        poynting_patch_force_summary((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))


if __name__ == "__main__":
    test_plane_wave_intensity_accepts_rms_and_peak_conventions()
    test_radiation_pressure_absorber_and_reflector_limits()
    test_radiation_force_from_power_is_area_integrated_pressure()
    test_radiation_force_from_normal_scattering_matches_momentum_balance()
    test_two_port_scattering_momentum_force_vector()
    test_oblique_radiation_pressure_reduces_to_normal_incidence()
    test_oblique_radiation_pressure_splits_absorbed_tangential_momentum()
    test_poynting_patch_force_vector_matches_oblique_summary()
    test_radiation_pressure_rejects_invalid_inputs()
    print("[OK] radiation-pressure helpers validated.")
