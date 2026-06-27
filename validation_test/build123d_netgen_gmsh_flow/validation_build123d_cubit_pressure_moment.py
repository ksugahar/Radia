"""Validation-class build123d analytic box pressure force/moment.

The CAD-side analytic reference is a constant pressure on each planar box face:

    F = p n A,
    M = (face_center - pivot) x F.

If ``--vol`` is provided, the same pressure table and pivot are applied to named
Netgen `.vol` boundaries and compared by face name.

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_moment.py
    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_moment.py --vol C:\\temp\\box.vol
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

from radia_mcp.build123d.modeling import box_face_pressure_moment_rows  # noqa: E402
from radia_mcp.radia_ngsolve.netgen_vol import read_netgen_tri_tet_vol  # noqa: E402


SUMMARY_JSON = HERE / "validation_build123d_cubit_pressure_moment_summary.json"
SIZE = (2.0, 3.0, 5.0)
PIVOT = (-1.0, -1.5, 0.0)
PRESSURE = {"zmax": 2.0}
UNIFORM_PRESSURE = {
    "xmin": 2.0,
    "xmax": 2.0,
    "ymin": 2.0,
    "ymax": 2.0,
    "zmin": 2.0,
    "zmax": 2.0,
}


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _sum_rows(rows, key):
    return [sum(row[key][axis] for row in rows) for axis in range(3)]


def _compare_moment_rows(reference_rows, measured_rows):
    measured_by_name = {row["name"]: row for row in measured_rows}
    rows = []
    for ref in reference_rows:
        measured = measured_by_name.get(ref["name"])
        if measured is None:
            rows.append({"name": ref["name"], "passed": False, "reason": "missing measured row"})
            continue
        force_error = _norm([
            measured["force_N"][axis] - ref["force_N"][axis]
            for axis in range(3)
        ])
        moment_error = _norm([
            measured["moment_about_pivot_Nm"][axis] - ref["moment_about_pivot_Nm"][axis]
            for axis in range(3)
        ])
        passed = force_error <= 1.0e-8 and moment_error <= 1.0e-8
        rows.append({
            "name": ref["name"],
            "reference_force_N": ref["force_N"],
            "measured_force_N": measured["force_N"],
            "force_abs_error_N": force_error,
            "reference_moment_Nm": ref["moment_about_pivot_Nm"],
            "measured_moment_Nm": measured["moment_about_pivot_Nm"],
            "moment_abs_error_Nm": moment_error,
            "passed": passed,
            "reason": "ok" if passed else "outside tolerance",
        })
    return rows


def build_summary(vol_path: Path | None = None) -> dict:
    reference_rows = box_face_pressure_moment_rows(
        SIZE,
        PRESSURE,
        default_pressure=0.0,
        pivot_m=PIVOT,
    )
    uniform_rows = box_face_pressure_moment_rows(SIZE, UNIFORM_PRESSURE, pivot_m=PIVOT)
    zmax = {row["name"]: row for row in reference_rows}["zmax"]
    checks = {
        "reference_uniform_total_force_N": _sum_rows(uniform_rows, "force_N"),
        "reference_uniform_total_moment_Nm": _sum_rows(uniform_rows, "moment_about_pivot_Nm"),
        "reference_zmax_force_N": zmax["force_N"],
        "reference_zmax_moment_Nm": zmax["moment_about_pivot_Nm"],
    }
    if _norm(checks["reference_uniform_total_force_N"]) > 1.0e-14:
        raise AssertionError("analytic uniform pressure should have zero net force")
    if _norm(checks["reference_uniform_total_moment_Nm"]) > 1.0e-14:
        raise AssertionError("analytic uniform pressure should have zero net moment")
    if zmax["force_N"] != (0.0, 0.0, 12.0):
        raise AssertionError("analytic zmax force drifted")
    if zmax["moment_about_pivot_Nm"] != (18.0, -12.0, 0.0):
        raise AssertionError("analytic zmax moment drifted")

    vol = {"provided": False, "path": None}
    if vol_path is not None:
        mesh = read_netgen_tri_tet_vol(vol_path)
        measured_rows = list(mesh.boundary_pressure_force_moment_rows(
            PRESSURE,
            default_pressure=0.0,
            pivot_m=PIVOT,
        ))
        comparison_rows = _compare_moment_rows(reference_rows, measured_rows)
        checks.update({
            "vol_max_force_abs_error_N": max(row.get("force_abs_error_N", math.inf) for row in comparison_rows),
            "vol_max_moment_abs_error_Nm": max(row.get("moment_abs_error_Nm", math.inf) for row in comparison_rows),
            "vol_all_passed": all(row["passed"] for row in comparison_rows),
        })
        if not checks["vol_all_passed"]:
            raise AssertionError("external .vol pressure moment comparison failed")
        vol = {
            "provided": True,
            "path": str(vol_path),
            "mesh_summary": mesh.summary(),
            "pressure_moment_rows": measured_rows,
            "comparison_rows": comparison_rows,
        }

    return {
        "kind": "build123d_box_pressure_moment_validation",
        "validation_class": True,
        "size": SIZE,
        "pivot_m": PIVOT,
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
    print("[build123d pressure force/moment]")
    print(f"  reference uniform force:  {summary['checks']['reference_uniform_total_force_N']} N")
    print(f"  reference uniform moment: {summary['checks']['reference_uniform_total_moment_Nm']} N m")
    print(f"  reference zmax force:    {summary['checks']['reference_zmax_force_N']} N")
    print(f"  reference zmax moment:   {summary['checks']['reference_zmax_moment_Nm']} N m")
    if summary["vol"]["provided"]:
        print(f"  vol max force error:  {summary['checks']['vol_max_force_abs_error_N']:.3e} N")
        print(f"  vol max moment error: {summary['checks']['vol_max_moment_abs_error_Nm']:.3e} N m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
