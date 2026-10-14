"""Validation-class vector radiation force from a Poynting vector.

This example turns a 3D time-average Poynting vector into the force on a flat
absorbing or reflecting patch.  It cross-checks the vector result against the
existing scalar oblique-incidence radiation-pressure identities.

Run:

    python validation_test/rf_waveguide/validation_poynting_patch_force_vector.py
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
    poynting_patch_force_summary,
)


OUT_JSON = Path(__file__).with_name("validation_poynting_patch_force_vector_summary.json")

INTENSITY = 9.0
AREA = 0.5
ANGLE_DEG = 60.0
NORMAL = (0.0, 0.0, 1.0)


def _assert_close(value: float, expected: float, rel: float = 1.0e-12, abs_tol: float = 1.0e-18) -> None:
    if abs(value - expected) > max(abs_tol, rel * max(abs(value), abs(expected), 1.0)):
        raise AssertionError(f"{value!r} != {expected!r}")


def build_summary() -> dict:
    angle = math.radians(ANGLE_DEG)
    c = math.cos(angle)
    s = math.sin(angle)
    poynting = (INTENSITY * s, 0.0, -INTENSITY * c)
    vector_absorber = poynting_patch_force_summary(
        poynting,
        NORMAL,
        area_m2=AREA,
        absorptance=1.0,
        reflectance=0.0,
    )
    vector_reflector = poynting_patch_force_summary(
        poynting,
        NORMAL,
        area_m2=AREA,
        absorptance=0.0,
        reflectance=1.0,
    )
    scalar_absorber = oblique_radiation_pressure_summary(
        INTENSITY,
        angle,
        area_m2=AREA,
        absorptance=1.0,
        reflectance=0.0,
    )
    scalar_reflector = oblique_radiation_pressure_summary(
        INTENSITY,
        angle,
        area_m2=AREA,
        absorptance=0.0,
        reflectance=1.0,
    )
    checks = {
        "incident_power_on_patch_W": vector_absorber["incident_power_on_patch_W"],
        "expected_incident_power_on_patch_W": INTENSITY * AREA * c,
        "absorber_normal_force_abs_error_N": abs(
            vector_absorber["normal_force_into_surface_N"] - scalar_absorber["normal_force_N"]
        ),
        "absorber_tangential_force_abs_error_N": abs(
            vector_absorber["tangential_force_magnitude_N"] - scalar_absorber["tangential_force_N"]
        ),
        "reflector_normal_force_abs_error_N": abs(
            vector_reflector["normal_force_into_surface_N"] - scalar_reflector["normal_force_N"]
        ),
        "reflector_tangential_force_N": vector_reflector["tangential_force_magnitude_N"],
        "absorber_force_vector_N": vector_absorber["force_N"],
        "reflector_force_vector_N": vector_reflector["force_N"],
        "speed_m_per_s": C0,
    }

    _assert_close(
        checks["incident_power_on_patch_W"],
        checks["expected_incident_power_on_patch_W"],
    )
    _assert_close(checks["absorber_normal_force_abs_error_N"], 0.0)
    _assert_close(checks["absorber_tangential_force_abs_error_N"], 0.0)
    _assert_close(checks["reflector_normal_force_abs_error_N"], 0.0)
    _assert_close(checks["reflector_tangential_force_N"], 0.0)

    return {
        "kind": "poynting_patch_force_vector",
        "validation_class": True,
        "force_learning": "3D Poynting vectors reduce to normal and tangential radiation-force components",
        "parameters": {
            "intensity_W_per_m2": INTENSITY,
            "area_m2": AREA,
            "incidence_angle_deg": ANGLE_DEG,
            "surface_normal": list(NORMAL),
            "poynting_W_per_m2": list(poynting),
        },
        "cases": {
            "absorber": vector_absorber,
            "reflector": vector_reflector,
            "scalar_absorber_reference": scalar_absorber,
            "scalar_reflector_reference": scalar_reflector,
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Poynting patch force vector]")
    print(f"  incident_power_on_patch_W: {checks['incident_power_on_patch_W']:.12g}")
    print(f"  absorber_force_vector_N: {checks['absorber_force_vector_N']}")
    print(f"  reflector_force_vector_N: {checks['reflector_force_vector_N']}")
    print(f"  absorber_normal_force_abs_error_N: {checks['absorber_normal_force_abs_error_N']:.3e}")
    print(f"  absorber_tangential_force_abs_error_N: {checks['absorber_tangential_force_abs_error_N']:.3e}")
    print(f"  reflector_normal_force_abs_error_N: {checks['reflector_normal_force_abs_error_N']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
