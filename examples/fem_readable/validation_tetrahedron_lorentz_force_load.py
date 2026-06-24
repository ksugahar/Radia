"""Validation-class P1 tetrahedron Lorentz body-force load.

This is the volume-force companion to the surface Maxwell-traction example.
For constant fields in one P1 tetrahedron,

    f = J x B,      F_e = volume * f,

and the consistent P1 load gives ``F_e / 4`` to each node.  The calculation is
small enough to mirror directly in MATLAB/Gypsilab teaching scripts.

Run:

    python examples/fem_readable/validation_tetrahedron_lorentz_force_load.py
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

from radia_mcp.radia_ngsolve.force import tetrahedron_lorentz_force_summary  # noqa: E402


OUT_JSON = HERE / "validation_tetrahedron_lorentz_force_load_summary.json"

VERTICES = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]
CURRENT_DENSITY = [2.0, 0.0, 0.0]
B_FIELD = [0.0, 3.0, 0.0]


def _assert_close_vector(actual: list[float], expected: list[float], tol: float = 1.0e-14) -> None:
    if len(actual) != len(expected):
        raise AssertionError("vector length mismatch")
    for a, e in zip(actual, expected):
        if abs(a - e) > tol:
            raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    row = tetrahedron_lorentz_force_summary(VERTICES, CURRENT_DENSITY, B_FIELD)
    checks = {
        "volume_m3": row["volume_m3"],
        "expected_volume_m3": 1.0 / 6.0,
        "force_density_N_per_m3": row["force_density_N_per_m3"],
        "expected_force_density_N_per_m3": [0.0, 0.0, 6.0],
        "integrated_force_N": row["integrated_force_N"],
        "expected_integrated_force_N": [0.0, 0.0, 1.0],
        "nodal_force_sum_N": [
            sum(node[axis] for node in row["nodal_force_loads_N"])
            for axis in range(3)
        ],
        "expected_node_force_N": [0.0, 0.0, 0.25],
    }
    if abs(checks["volume_m3"] - checks["expected_volume_m3"]) > 1.0e-14:
        raise AssertionError("tet volume mismatch")
    _assert_close_vector(checks["force_density_N_per_m3"], checks["expected_force_density_N_per_m3"])
    _assert_close_vector(checks["integrated_force_N"], checks["expected_integrated_force_N"])
    _assert_close_vector(checks["nodal_force_sum_N"], checks["expected_integrated_force_N"])
    for node_force in row["nodal_force_loads_N"]:
        _assert_close_vector(node_force, checks["expected_node_force_N"])
    return {
        "kind": "p1_tetrahedron_lorentz_force_load",
        "validation_class": True,
        "vertices": VERTICES,
        "current_density_A_per_m2": CURRENT_DENSITY,
        "B_T": B_FIELD,
        "row": row,
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[P1 tetrahedron Lorentz force load]")
    print(f"  volume_m3: {checks['volume_m3']}")
    print(f"  force_density_N_per_m3: {checks['force_density_N_per_m3']}")
    print(f"  integrated_force_N: {checks['integrated_force_N']}")
    print(f"  expected_node_force_N: {checks['expected_node_force_N']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
