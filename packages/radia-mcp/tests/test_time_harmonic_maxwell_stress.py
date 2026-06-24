"""Time-harmonic Maxwell stress tensor helpers."""

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
    plane_wave_intensity_from_electric_field,
    radiation_pressure_from_intensity,
    time_average_maxwell_stress_tensor,
    time_average_maxwell_traction,
    time_average_maxwell_traction_summary,
)


def _assert_matrix_close(a, b, abs_tol=1.0e-14):
    assert len(a) == len(b)
    for row_a, row_b in zip(a, b):
        assert row_a == pytest.approx(row_b, abs=abs_tol)


def test_peak_plane_wave_stress_matches_poynting_momentum_flux():
    e0 = 1000.0
    E = (e0, 0.0, 0.0)
    H = (0.0, e0 / ETA0, 0.0)
    intensity = plane_wave_intensity_from_electric_field(e0, amplitude="peak")
    pressure = radiation_pressure_from_intensity(intensity)
    tensor = time_average_maxwell_stress_tensor(E, H, amplitude="peak")

    assert pressure == pytest.approx(intensity / C0)
    assert tensor[0][0] == pytest.approx(0.0, abs=1.0e-14)
    assert tensor[1][1] == pytest.approx(0.0, abs=1.0e-14)
    assert tensor[2][2] == pytest.approx(-pressure)
    assert tensor[0][1] == pytest.approx(0.0, abs=1.0e-18)


def test_rms_and_peak_phasors_give_same_average_stress():
    e_peak = 1000.0
    e_rms = e_peak / math.sqrt(2.0)
    peak = time_average_maxwell_stress_tensor(
        (e_peak, 0.0, 0.0),
        (0.0, e_peak / ETA0, 0.0),
        amplitude="peak",
    )
    rms = time_average_maxwell_stress_tensor(
        (e_rms, 0.0, 0.0),
        (0.0, e_rms / ETA0, 0.0),
        amplitude="rms",
    )

    _assert_matrix_close(rms, peak)


def test_time_average_traction_summary_decomposes_patch_force():
    e0 = 1000.0
    area = 0.25
    intensity = plane_wave_intensity_from_electric_field(e0, amplitude="peak")
    pressure = radiation_pressure_from_intensity(intensity)
    summary = time_average_maxwell_traction_summary(
        (e0, 0.0, 0.0),
        (0.0, e0 / ETA0, 0.0),
        (0.0, 0.0, 1.0),
        area_m2=area,
    )
    traction = time_average_maxwell_traction(
        (e0, 0.0, 0.0),
        (0.0, e0 / ETA0, 0.0),
        (0.0, 0.0, 1.0),
    )

    assert traction == pytest.approx([0.0, 0.0, -pressure])
    assert summary["normal_traction_Pa"] == pytest.approx(-pressure)
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(0.0)
    assert summary["force_N"] == pytest.approx([0.0, 0.0, -pressure * area])
    assert summary["average_energy_density_J_per_m3"] == pytest.approx(pressure)


def test_circular_polarized_plane_wave_keeps_same_momentum_flux():
    e0 = 1000.0
    scale = e0 / math.sqrt(2.0)
    E = (scale, 1j * scale, 0.0)
    H = (-1j * scale / ETA0, scale / ETA0, 0.0)
    tensor = time_average_maxwell_stress_tensor(E, H, amplitude="peak")
    intensity = plane_wave_intensity_from_electric_field(e0, amplitude="peak")
    pressure = radiation_pressure_from_intensity(intensity)

    assert tensor[0][0] == pytest.approx(0.0, abs=1.0e-14)
    assert tensor[1][1] == pytest.approx(0.0, abs=1.0e-14)
    assert tensor[2][2] == pytest.approx(-pressure)


def test_time_harmonic_stress_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        time_average_maxwell_stress_tensor((1.0, 0.0), (0.0, 1.0, 0.0))
    with pytest.raises(ValueError):
        time_average_maxwell_stress_tensor((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), eps=0.0)
    with pytest.raises(ValueError):
        time_average_maxwell_stress_tensor((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), mu=0.0)
    with pytest.raises(ValueError):
        time_average_maxwell_stress_tensor(
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            amplitude="phasor",
        )
    with pytest.raises(ValueError):
        time_average_maxwell_traction_summary(
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            area_m2=-1.0,
        )


if __name__ == "__main__":
    test_peak_plane_wave_stress_matches_poynting_momentum_flux()
    test_rms_and_peak_phasors_give_same_average_stress()
    test_time_average_traction_summary_decomposes_patch_force()
    test_circular_polarized_plane_wave_keeps_same_momentum_flux()
    test_time_harmonic_stress_rejects_invalid_inputs()
    print("[OK] time-harmonic Maxwell stress helpers validated.")
