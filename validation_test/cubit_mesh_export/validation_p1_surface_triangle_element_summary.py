"""Validation-class P1 surface-triangle element teaching block.

Run:

    python validation_test/cubit_mesh_export/validation_p1_surface_triangle_element_summary.py

This example records the local geometry, surface gradients, P1 stiffness,
consistent mass, constant-load vector, and one-based sparse triplets for one
boundary triangle.  It is the small readable block used by FEM/BEM trace and
boundary assembly examples.
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

from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_surface_triangle_element_summary,
)


OUT_JSON = HERE / "validation_p1_surface_triangle_element_summary.json"
TRIANGLE = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0)]


def build_summary() -> dict[str, object]:
    summary = p1_surface_triangle_element_summary(
        TRIANGLE,
        density=3.0,
        coeff=2.0,
        source=5.0,
    )
    checks = {
        "area": summary["area"],
        "unit_normal": summary["unit_normal"],
        "mass_row_sums": summary["mass_row_sums"],
        "mass_integral_of_one": summary["mass_integral_of_one"],
        "constant_load_integral": summary["constant_load_integral"],
        "stiffness_nullspace_residual": summary["stiffness_nullspace_residual"],
        "gradient_partition_residual": summary["gradient_partition_residual"],
        "mass_triplet_count": len(summary["mass_triplets_1based"]),
        "stiffness_triplet_count": len(summary["stiffness_triplets_1based"]),
    }

    assert abs(float(checks["area"]) - 1.0) < 1.0e-12
    assert checks["unit_normal"] == (0.0, 0.0, 1.0)
    assert max(abs(value - 1.0) for value in checks["mass_row_sums"]) < 1.0e-12
    assert abs(float(checks["mass_integral_of_one"]) - 3.0) < 1.0e-12
    assert abs(float(checks["constant_load_integral"]) - 5.0) < 1.0e-12
    assert checks["stiffness_nullspace_residual"] < 1.0e-15
    assert checks["gradient_partition_residual"] < 1.0e-15
    assert checks["mass_triplet_count"] == 9
    assert checks["stiffness_triplet_count"] == 9

    return {
        "kind": "p1_surface_triangle_element_summary_validation",
        "validation_class": True,
        "learning_theme": (
            "P1 boundary triangles should expose geometry, mass, stiffness, "
            "load, and one-based sparse triplets as one readable assembly block"
        ),
        "checks": checks,
        "element_summary": summary,
    }


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = _json_clean(build_summary())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[P1 surface triangle element summary]")
    print(
        f"  area={checks['area']:.12g} mass_integral={checks['mass_integral_of_one']:.12g} "
        f"load_integral={checks['constant_load_integral']:.12g}"
    )
    print(
        f"  stiffness_nullspace_residual={checks['stiffness_nullspace_residual']:.3e} "
        f"gradient_partition_residual={checks['gradient_partition_residual']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
