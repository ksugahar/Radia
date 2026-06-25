"""Validation-class acoustic impedance radiation pressure gate.

Run:

    python examples/acoustic_bem/validation_acoustic_impedance_radiation_pressure.py

For a local acoustic impedance load, the pressure reflection coefficient gives
both energy absorption and normal momentum transfer:

    p_rad = (1 + R) I_inc / c = (A + 2 R) I_inc / c

where ``R=|Gamma|^2`` and ``A=1-R``.  The matched absorber gives ``I/c``; a
lossless reflector gives ``2I/c``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import acoustic_impedance_radiation_pressure_summary  # noqa: E402


OUT_JSON = HERE / "validation_acoustic_impedance_radiation_pressure_summary.json"
RHO = 1.2041
C = 343.0
Z0 = RHO * C
INCIDENT_PRESSURE = 2.0
AREA = 0.5


def _complex_record(value):
    z = complex(value)
    return {
        "real": 0.0 if z.real == 0.0 else z.real,
        "imag": 0.0 if z.imag == 0.0 else z.imag,
        "abs": abs(z),
    }


def _json_clean(value):
    if isinstance(value, complex):
        return _complex_record(value)
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-18) -> float:
    error = abs(actual - expected)
    if error > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    incident_intensity = INCIDENT_PRESSURE * INCIDENT_PRESSURE / (2.0 * Z0)
    matched = acoustic_impedance_radiation_pressure_summary(
        Z0,
        area=AREA,
        incident_pressure=INCIDENT_PRESSURE,
        rho=RHO,
        c=C,
    )
    double_resistance = acoustic_impedance_radiation_pressure_summary(
        2.0 * Z0,
        area=AREA,
        incident_pressure=INCIDENT_PRESSURE,
        rho=RHO,
        c=C,
    )
    reactive = acoustic_impedance_radiation_pressure_summary(
        1.0j * Z0,
        area=AREA,
        incident_pressure=INCIDENT_PRESSURE,
        rho=RHO,
        c=C,
    )

    errors = {
        "matched_pressure_error_Pa": _assert_close(
            matched["normal_momentum_pressure_Pa"],
            incident_intensity / C,
        ),
        "double_resistance_pressure_error_Pa": _assert_close(
            double_resistance["normal_momentum_pressure_Pa"],
            (10.0 / 9.0) * incident_intensity / C,
        ),
        "reactive_pressure_error_Pa": _assert_close(
            reactive["normal_momentum_pressure_Pa"],
            2.0 * incident_intensity / C,
        ),
        "double_resistance_force_balance_error_N": abs(double_resistance["force_balance_residual_N"]),
    }

    checks = {
        "incident_intensity_W_per_m2": incident_intensity,
        "matched_pressure_Pa": matched["normal_momentum_pressure_Pa"],
        "matched_expected_pressure_Pa": incident_intensity / C,
        "double_resistance_reflectance": double_resistance["power_reflection_coefficient"],
        "double_resistance_absorption": double_resistance["absorption_coefficient"],
        "double_resistance_pressure_Pa": double_resistance["normal_momentum_pressure_Pa"],
        "double_resistance_expected_pressure_Pa": (10.0 / 9.0) * incident_intensity / C,
        "reactive_reflectance": reactive["power_reflection_coefficient"],
        "reactive_pressure_Pa": reactive["normal_momentum_pressure_Pa"],
        "reactive_expected_pressure_Pa": 2.0 * incident_intensity / C,
        "max_abs_error": max(errors.values()),
    }

    return {
        "kind": "acoustic_impedance_radiation_pressure_validation",
        "validation_class": True,
        "force_learning": (
            "acoustic impedance reflection maps to normal momentum pressure by "
            "(1+R)I/c=(A+2R)I/c"
        ),
        "rho": RHO,
        "c": C,
        "area": AREA,
        "incident_pressure": INCIDENT_PRESSURE,
        "checks": checks,
        "errors": errors,
        "rows": {
            "matched_absorber": matched,
            "double_resistance": double_resistance,
            "pure_reactive_reflector": reactive,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(_json_clean(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[acoustic impedance radiation pressure]")
    print(f"  incident_intensity={checks['incident_intensity_W_per_m2']:.12g} W/m^2")
    print(f"  matched_pressure={checks['matched_pressure_Pa']:.12e} Pa")
    print(
        "  Z=2Z0: "
        f"R={checks['double_resistance_reflectance']:.12g}, "
        f"A={checks['double_resistance_absorption']:.12g}, "
        f"pressure={checks['double_resistance_pressure_Pa']:.12e} Pa"
    )
    print(f"  reactive_pressure={checks['reactive_pressure_Pa']:.12e} Pa")
    print(f"  max_abs_error={checks['max_abs_error']:.3e}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
