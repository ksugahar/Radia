"""Validation-class time-harmonic Maxwell stress identities.

Frequency-domain FEM packages return complex phasors.  This example pins the
local post-processing convention: peak phasors need a ``1/2`` time-average
factor, RMS phasors do not, and a plane wave carries momentum flux ``I/c``.
The signs are those of the field stress ``<T> n``; a receiving surface takes
the opposite normal force from the incident field.
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from result_metadata import add_result_metadata  # noqa: E402

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    ETA0,
    plane_wave_intensity_from_electric_field,
    radiation_pressure_from_intensity,
    time_average_maxwell_stress_tensor,
    time_average_maxwell_traction_summary,
)


OUT_JSON = Path(__file__).with_name("validation_time_harmonic_maxwell_stress_summary.json")


def _max_abs_matrix_delta(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(len(a)) for j in range(len(a[i])))


def main():
    e_peak = 1000.0
    e_rms = e_peak / math.sqrt(2.0)
    area = 0.25
    intensity = plane_wave_intensity_from_electric_field(e_peak, amplitude="peak")
    pressure = radiation_pressure_from_intensity(intensity)

    peak = time_average_maxwell_traction_summary(
        (e_peak, 0.0, 0.0),
        (0.0, e_peak / ETA0, 0.0),
        (0.0, 0.0, 1.0),
        area_m2=area,
        amplitude="peak",
    )
    rms = time_average_maxwell_traction_summary(
        (e_rms, 0.0, 0.0),
        (0.0, e_rms / ETA0, 0.0),
        (0.0, 0.0, 1.0),
        area_m2=area,
        amplitude="rms",
    )
    circular_e = (e_peak / math.sqrt(2.0), 1j * e_peak / math.sqrt(2.0), 0.0)
    circular_h = (
        -1j * e_peak / (math.sqrt(2.0) * ETA0),
        e_peak / (math.sqrt(2.0) * ETA0),
        0.0,
    )
    circular_tensor = time_average_maxwell_stress_tensor(
        circular_e, circular_h, amplitude="peak"
    )

    checks = {
        "E_peak_V_per_m": e_peak,
        "E_rms_V_per_m": e_rms,
        "intensity_W_per_m2": intensity,
        "radiation_pressure_Pa": pressure,
        "momentum_flux_I_over_c_Pa": intensity / C0,
        "peak_Tzz_Pa": peak["stress_tensor_Pa"][2][2],
        "rms_Tzz_Pa": rms["stress_tensor_Pa"][2][2],
        "circular_Tzz_Pa": circular_tensor[2][2],
        "peak_normal_traction_Pa": peak["normal_traction_Pa"],
        "receiving_surface_force_N": -peak["force_N"][2],
        "peak_rms_tensor_max_abs_delta": _max_abs_matrix_delta(
            peak["stress_tensor_Pa"], rms["stress_tensor_Pa"]
        ),
        "peak_pressure_abs_error": abs(-peak["normal_traction_Pa"] - pressure),
        "rms_pressure_abs_error": abs(-rms["normal_traction_Pa"] - pressure),
        "circular_pressure_abs_error": abs(-circular_tensor[2][2] - pressure),
        "force_pressure_identity_abs_error": abs(-peak["force_N"][2] - pressure * area),
    }

    assert checks["peak_rms_tensor_max_abs_delta"] < 1.0e-14
    assert checks["peak_pressure_abs_error"] < 1.0e-14
    assert checks["rms_pressure_abs_error"] < 1.0e-14
    assert checks["circular_pressure_abs_error"] < 1.0e-14
    assert checks["force_pressure_identity_abs_error"] < 1.0e-14

    summary = {
        "kind": "time_harmonic_maxwell_stress_validation",
        "validation_class": True,
        "convention": "peak phasors use 1/2 time-average factor; RMS phasors use factor 1",
        "area_m2": area,
        "peak_summary": peak,
        "rms_summary": rms,
        "circular_stress_tensor_Pa": circular_tensor,
        "checks": checks,
    }
    summary = add_result_metadata(summary, __file__)
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[time-harmonic Maxwell stress]")
    print(
        f"  E_peak={e_peak:.6g} V/m I={intensity:.12g} W/m2 "
        f"pressure={pressure:.12g} Pa"
    )
    print(
        f"  Tzz_peak={checks['peak_Tzz_Pa']:.12g} Pa "
        f"Tzz_rms={checks['rms_Tzz_Pa']:.12g} Pa "
        f"Tzz_circular={checks['circular_Tzz_Pa']:.12g} Pa"
    )
    print(
        f"  receiving force on {area:.3g} m2 = "
        f"{checks['receiving_surface_force_N']:.12g} N"
    )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
