"""Validation-class build123d analytic box pressure force.

The CAD-side analytic reference is ``pressure * vector_area`` on the six faces
of an axis-aligned box.  If ``--vol`` is provided, the same pressure table is
applied to the named Netgen `.vol` boundaries and compared by face name.

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_force.py
    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_force.py --vol C:\\temp\\box.vol
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.build123d.modeling import box_face_pressure_force_rows  # noqa: E402
from radia_mcp.radia_ngsolve.netgen_vol import read_netgen_tri_tet_vol  # noqa: E402


SUMMARY_JSON = HERE / "validation_build123d_cubit_pressure_force_summary.json"
SIZE = (2.0, 3.0, 5.0)
PRESSURE = {"zmax": 2.0}


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _sum_forces(rows):
    return [sum(row["force_N"][axis] for row in rows) for axis in range(3)]


def _compare_force_rows(reference_rows, measured_rows):
    measured_by_name = {row["name"]: row for row in measured_rows}
    rows = []
    for ref in reference_rows:
        measured = measured_by_name.get(ref["name"])
        if measured is None:
            rows.append({"name": ref["name"], "passed": False, "reason": "missing measured row"})
            continue
        error = _norm([
            measured["force_N"][axis] - ref["force_N"][axis]
            for axis in range(3)
        ])
        rows.append({
            "name": ref["name"],
            "reference_force_N": ref["force_N"],
            "measured_force_N": measured["force_N"],
            "force_abs_error_N": error,
            "passed": error <= 1.0e-8,
            "reason": "ok" if error <= 1.0e-8 else "outside tolerance",
        })
    return rows


def build_summary(vol_path: Path | None = None) -> dict:
    reference_rows = box_face_pressure_force_rows(SIZE, PRESSURE, default_pressure=0.0)
    checks = {
        "reference_total_force_N": _sum_forces(reference_rows),
        "reference_zmax_force_N": {row["name"]: row for row in reference_rows}["zmax"]["force_N"],
    }
    if checks["reference_zmax_force_N"] != (0.0, 0.0, 12.0):
        raise AssertionError("analytic zmax force drifted")

    vol = {"provided": False, "path": None}
    if vol_path is not None:
        mesh = read_netgen_tri_tet_vol(vol_path)
        measured_rows = list(mesh.boundary_pressure_force_rows(PRESSURE, default_pressure=0.0))
        comparison_rows = _compare_force_rows(reference_rows, measured_rows)
        checks.update({
            "vol_max_force_abs_error_N": max(row.get("force_abs_error_N", math.inf) for row in comparison_rows),
            "vol_all_passed": all(row["passed"] for row in comparison_rows),
        })
        if not checks["vol_all_passed"]:
            raise AssertionError("external .vol pressure force comparison failed")
        vol = {
            "provided": True,
            "path": str(vol_path),
            "mesh_summary": mesh.summary(),
            "pressure_force_rows": measured_rows,
            "comparison_rows": comparison_rows,
        }

    return {
        "kind": "build123d_box_pressure_force_validation",
        "validation_class": True,
        "size": SIZE,
        "pressure_by_face_Pa": PRESSURE,
        "reference_rows": reference_rows,
        "checks": checks,
        "vol": vol,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=SUMMARY_JSON)
    args = parser.parse_args()

    summary = build_summary(args.vol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("[build123d pressure force]")
    print(f"  reference total force: {summary['checks']['reference_total_force_N']} N")
    print(f"  reference zmax force: {summary['checks']['reference_zmax_force_N']} N")
    if summary["vol"]["provided"]:
        print(f"  vol max force error: {summary['checks']['vol_max_force_abs_error_N']:.3e} N")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
