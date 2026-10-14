"""Validation-class oblique-incidence radiation pressure.

For a plane wave incident on a flat patch at angle ``theta`` from the surface
normal, the normal force scales as ``cos(theta)^2``.  Absorption also transfers
tangential momentum; specular reflection changes only the normal momentum.

Run:

    python validation_test/rf_waveguide/validation_oblique_radiation_pressure.py
"""

from __future__ import annotations

import argparse
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
    oblique_radiation_pressure_summary,
    radiation_pressure_summary,
)


OUT_JSON = Path(__file__).with_name("validation_oblique_radiation_pressure_summary.json")

INTENSITY = 12.0
AREA = 0.75
ANGLES_DEG = (0.0, 30.0, 60.0)


def _assert_close(value: float, expected: float, rel: float = 1.0e-12, abs_tol: float = 1.0e-18) -> None:
    if abs(value - expected) > max(abs_tol, rel * max(abs(value), abs(expected), 1.0)):
        raise AssertionError(f"{value!r} != {expected!r}")


def build_rows() -> list[dict]:
    rows = []
    for angle_deg in ANGLES_DEG:
        angle = math.radians(angle_deg)
        c = math.cos(angle)
        s = math.sin(angle)
        absorber = oblique_radiation_pressure_summary(
            INTENSITY,
            angle,
            area_m2=AREA,
            absorptance=1.0,
            reflectance=0.0,
        )
        reflector = oblique_radiation_pressure_summary(
            INTENSITY,
            angle,
            area_m2=AREA,
            absorptance=0.0,
            reflectance=1.0,
        )
        rows.append({
            "angle_deg": angle_deg,
            "cos_angle": c,
            "sin_angle": s,
            "absorber": absorber,
            "reflector": reflector,
            "expected_absorber_normal_force_N": INTENSITY * AREA * c * c / C0,
            "expected_absorber_tangential_force_N": INTENSITY * AREA * s * c / C0,
            "expected_reflector_normal_force_N": 2.0 * INTENSITY * AREA * c * c / C0,
            "normal_force_ratio_reflector_over_absorber": (
                reflector["normal_force_N"] / absorber["normal_force_N"]
                if absorber["normal_force_N"] else math.inf
            ),
        })
    return rows


def validate(rows: list[dict]) -> dict:
    normal = radiation_pressure_summary(INTENSITY, area_m2=AREA)
    errors = []
    for row in rows:
        absorber = row["absorber"]
        reflector = row["reflector"]
        checks = [
            absorber["normal_force_N"] - row["expected_absorber_normal_force_N"],
            absorber["tangential_force_N"] - row["expected_absorber_tangential_force_N"],
            reflector["normal_force_N"] - row["expected_reflector_normal_force_N"],
            reflector["tangential_force_N"],
        ]
        errors.extend(abs(value) for value in checks)
    zero_deg = rows[0]
    _assert_close(zero_deg["absorber"]["normal_force_N"], normal["force_N"])
    _assert_close(zero_deg["absorber"]["tangential_force_N"], 0.0)
    _assert_close(zero_deg["normal_force_ratio_reflector_over_absorber"], 2.0)
    checks = {
        "normal_incidence_force_N": normal["force_N"],
        "max_force_component_abs_error_N": max(errors),
        "angle_60_absorber_normal_over_angle_0": (
            rows[2]["absorber"]["normal_force_N"] / rows[0]["absorber"]["normal_force_N"]
        ),
        "angle_60_absorber_tangential_force_N": rows[2]["absorber"]["tangential_force_N"],
        "reflector_tangential_force_max_abs_N": max(
            abs(row["reflector"]["tangential_force_N"]) for row in rows
        ),
    }
    _assert_close(checks["angle_60_absorber_normal_over_angle_0"], 0.25)
    _assert_close(checks["max_force_component_abs_error_N"], 0.0)
    _assert_close(checks["reflector_tangential_force_max_abs_N"], 0.0)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "oblique_radiation_pressure_validation",
        "validation_class": True,
        "intensity_W_per_m2": INTENSITY,
        "area_m2": AREA,
        "speed_m_per_s": C0,
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[Oblique radiation pressure]")
    for row in rows:
        print(
            f"  theta={row['angle_deg']:>4.0f} deg  "
            f"F_abs_n={row['absorber']['normal_force_N']:.6e} N  "
            f"F_abs_t={row['absorber']['tangential_force_N']:.6e} N  "
            f"F_ref_n={row['reflector']['normal_force_N']:.6e} N"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
