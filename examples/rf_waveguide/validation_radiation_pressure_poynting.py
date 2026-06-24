"""Validation-class RF radiation pressure from Poynting power.

This example links waveguide port normalization to electromagnetic force:
integrated time-average power carries momentum ``P/c``.  A matched absorber
takes ``P/c`` normal force, while a short/ideal mirror takes ``2P/c``.  The
calculation is closed-form and public-safe; no commercial RF solver data is
needed.
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    ETA0,
    plane_wave_intensity_from_electric_field,
    radiation_force_from_power,
    radiation_pressure_summary,
)
from radia_mcp.radia_ngsolve.waveguide import rectangular_waveguide_te10_port_normalization  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_radiation_pressure_poynting_summary.json")

WR90_A = 0.02286
WR90_B = 0.01016
FREQUENCY = 10.0e9


def _assert_close(value, expected, rel=1e-12, abs_tol=1e-18):
    assert abs(value - expected) <= max(abs_tol, rel * max(abs(expected), 1.0))


def main():
    area = WR90_A * WR90_B
    cases = []
    for power in (0.25, 1.0, 4.0):
        port = rectangular_waveguide_te10_port_normalization(
            FREQUENCY, WR90_A, WR90_B, power_w=power
        )
        intensity = power / area
        absorber = radiation_pressure_summary(
            intensity, area_m2=area, absorptance=1.0, reflectance=0.0
        )
        reflector = radiation_pressure_summary(
            intensity, area_m2=area, absorptance=0.0, reflectance=1.0
        )
        cases.append(
            {
                "power_W": power,
                "port": port,
                "aperture_average_intensity_W_per_m2": intensity,
                "absorber": absorber,
                "reflector": reflector,
                "absorber_force_power_identity_abs_error": abs(
                    absorber["force_N"] - radiation_force_from_power(power)
                ),
                "reflector_force_power_identity_abs_error": abs(
                    reflector["force_N"]
                    - radiation_force_from_power(power, absorptance=0.0, reflectance=1.0)
                ),
            }
        )

    one_w = cases[1]
    four_w = cases[2]
    intensity = one_w["aperture_average_intensity_W_per_m2"]
    equivalent_e_rms = math.sqrt(ETA0 * intensity)
    equivalent_e_peak = math.sqrt(2.0 * ETA0 * intensity)
    intensity_from_rms = plane_wave_intensity_from_electric_field(
        equivalent_e_rms, amplitude="rms"
    )
    intensity_from_peak = plane_wave_intensity_from_electric_field(
        equivalent_e_peak, amplitude="peak"
    )
    checks = {
        "aperture_area_m2": area,
        "one_w_average_intensity_W_per_m2": intensity,
        "one_w_absorber_force_N": one_w["absorber"]["force_N"],
        "one_w_reflector_force_N": one_w["reflector"]["force_N"],
        "one_w_absorber_pressure_Pa": one_w["absorber"]["pressure_Pa"],
        "one_w_reflector_pressure_Pa": one_w["reflector"]["pressure_Pa"],
        "force_ratio_reflector_over_absorber": (
            one_w["reflector"]["force_N"] / one_w["absorber"]["force_N"]
        ),
        "force_ratio_4w_over_1w_absorber": (
            four_w["absorber"]["force_N"] / one_w["absorber"]["force_N"]
        ),
        "equivalent_plane_wave_E_rms_V_per_m": equivalent_e_rms,
        "equivalent_plane_wave_E_peak_V_per_m": equivalent_e_peak,
        "intensity_from_rms_abs_error": abs(intensity_from_rms - intensity),
        "intensity_from_peak_abs_error": abs(intensity_from_peak - intensity),
        "max_port_power_abs_error_W": max(
            row["port"]["poynting_abs_error_W"] for row in cases
        ),
        "max_force_power_identity_abs_error": max(
            max(
                row["absorber_force_power_identity_abs_error"],
                row["reflector_force_power_identity_abs_error"],
            )
            for row in cases
        ),
    }

    _assert_close(checks["one_w_absorber_force_N"], 1.0 / C0)
    _assert_close(checks["one_w_reflector_force_N"], 2.0 / C0)
    _assert_close(checks["force_ratio_reflector_over_absorber"], 2.0)
    _assert_close(checks["force_ratio_4w_over_1w_absorber"], 4.0)
    assert checks["intensity_from_rms_abs_error"] < 1e-12
    assert checks["intensity_from_peak_abs_error"] < 1e-12
    assert checks["max_port_power_abs_error_W"] < 1e-12
    assert checks["max_force_power_identity_abs_error"] < 1e-18

    summary = {
        "kind": "rf_radiation_pressure_poynting_validation",
        "validation_class": True,
        "guide": {
            "name": "WR-90",
            "width_a_m": WR90_A,
            "height_b_m": WR90_B,
            "frequency_Hz": FREQUENCY,
        },
        "rows": cases,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[RF radiation pressure from Poynting power]")
    for row in cases:
        print(
            f"  P={row['power_W']:.3g} W  "
            f"Iavg={row['aperture_average_intensity_W_per_m2']:.6g} W/m2  "
            f"F_abs={row['absorber']['force_N']:.6g} N  "
            f"F_ref={row['reflector']['force_N']:.6g} N"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
