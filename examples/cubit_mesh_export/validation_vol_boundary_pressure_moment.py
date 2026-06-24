"""Validation-class `.vol` boundary pressure force/moment rows.

Named Cubit/Coreform sidesets exported as Netgen `.vol` already provide
oriented boundary triangle area vectors.  A constant scalar pressure on a
boundary then has

    dF = p n dS,
    dM = (r - pivot) x dF.

This example checks both the per-boundary moment rows and the generic
force/moment reducer used elsewhere in `radia_ngsolve.force`.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_pressure_moment.py
    python examples/cubit_mesh_export/validation_vol_boundary_pressure_moment.py --vol C:\\temp\\box.vol
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

from radia_mcp.radia_ngsolve.force import force_moment_resultant_summary  # noqa: E402
from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402
from validation_vol_boundary_normal_vectors import BOX_SIX_BOUNDARY_VOL  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_pressure_moment_summary.json")


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _sum_rows(rows, key):
    return [
        sum(row[key][axis] for row in rows)
        for axis in range(3)
    ]


def _resultant_from_rows(rows, pivot=(0.0, 0.0, 0.0)):
    active = [row for row in rows if row["centroid_m"] is not None]
    out = force_moment_resultant_summary(
        [row["centroid_m"] for row in active],
        [row["force_N"] for row in active],
        pivot_m=pivot,
    )
    return {
        "n_rows": out["n_rows"],
        "pivot_m": out["pivot_m"],
        "total_force": out["total_force"],
        "total_force_magnitude": out["total_force_magnitude"],
        "total_moment": out["total_moment"],
        "total_moment_magnitude": out["total_moment_magnitude"],
    }


def _record(label, mesh):
    uniform_pressure = 2.0
    uniform_rows = list(mesh.boundary_pressure_force_moment_rows({
        row["name"]: uniform_pressure
        for row in mesh.boundary_normal_summary_rows()
    }))
    zmax_rows = list(mesh.boundary_pressure_force_moment_rows({"zmax": 2.0}, default_pressure=0.0))
    shifted_rows = list(mesh.boundary_pressure_force_moment_rows(
        {"zmax": 2.0},
        default_pressure=0.0,
        pivot_m=(1.0, 1.5, 0.0),
    ))
    zmax_row = {row["name"]: row for row in zmax_rows}.get("zmax")
    zmax_shifted_row = {row["name"]: row for row in shifted_rows}.get("zmax")
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "uniform_pressure_Pa": uniform_pressure,
        "uniform_rows": uniform_rows,
        "uniform_total_force_N": _sum_rows(uniform_rows, "force_N"),
        "uniform_total_moment_Nm": _sum_rows(uniform_rows, "moment_about_pivot_Nm"),
        "uniform_resultant": _resultant_from_rows(uniform_rows),
        "zmax_row": zmax_row,
        "zmax_resultant": _resultant_from_rows(zmax_rows),
        "zmax_shifted_pivot_row": zmax_shifted_row,
        "zmax_shifted_resultant": _resultant_from_rows(shifted_rows, pivot=(1.0, 1.5, 0.0)),
    }


def build_summary(external_vol=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if external_vol is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(external_vol)))

    builtin = records[0]
    zmax = builtin["zmax_row"]
    shifted = builtin["zmax_shifted_pivot_row"]
    checks = {
        "builtin_uniform_total_force_norm_N": _norm(builtin["uniform_total_force_N"]),
        "builtin_uniform_total_moment_norm_Nm": _norm(builtin["uniform_total_moment_Nm"]),
        "builtin_zmax_centroid_m": zmax["centroid_m"],
        "builtin_zmax_force_N": zmax["force_N"],
        "builtin_zmax_moment_Nm": zmax["moment_about_pivot_Nm"],
        "builtin_zmax_shifted_moment_Nm": shifted["moment_about_pivot_Nm"],
        "builtin_zmax_resultant_force_error_N": _norm([
            builtin["zmax_resultant"]["total_force"][axis] - zmax["force_N"][axis]
            for axis in range(3)
        ]),
        "builtin_zmax_resultant_moment_error_Nm": _norm([
            builtin["zmax_resultant"]["total_moment"][axis] - zmax["moment_about_pivot_Nm"][axis]
            for axis in range(3)
        ]),
    }
    if checks["builtin_uniform_total_force_norm_N"] > 1.0e-14:
        raise AssertionError("uniform closed-box pressure should have zero net force")
    if checks["builtin_uniform_total_moment_norm_Nm"] > 1.0e-14:
        raise AssertionError("uniform closed-box pressure should have zero net moment")
    if _norm([zmax["force_N"][0], zmax["force_N"][1], zmax["force_N"][2] - 12.0]) > 1.0e-14:
        raise AssertionError("zmax force drifted")
    if _norm([
        zmax["moment_about_pivot_Nm"][0] - 18.0,
        zmax["moment_about_pivot_Nm"][1] + 12.0,
        zmax["moment_about_pivot_Nm"][2],
    ]) > 1.0e-14:
        raise AssertionError("zmax moment drifted")
    if _norm(shifted["moment_about_pivot_Nm"]) > 1.0e-14:
        raise AssertionError("zmax shifted-pivot moment should vanish")
    if checks["builtin_zmax_resultant_force_error_N"] > 1.0e-14:
        raise AssertionError("generic resultant force drifted")
    if checks["builtin_zmax_resultant_moment_error_Nm"] > 1.0e-14:
        raise AssertionError("generic resultant moment drifted")

    if external_vol is not None:
        external = records[-1]
        checks.update({
            "external_uniform_total_force_norm_N": _norm(external["uniform_total_force_N"]),
            "external_uniform_total_moment_norm_Nm": _norm(external["uniform_total_moment_Nm"]),
            "external_zmax_force_N": (
                external["zmax_row"].get("force_N") if external["zmax_row"] else None
            ),
            "external_zmax_moment_Nm": (
                external["zmax_row"].get("moment_about_pivot_Nm") if external["zmax_row"] else None
            ),
        })
        if checks["external_uniform_total_force_norm_N"] > 1.0e-8:
            raise AssertionError("external uniform pressure force did not cancel")
        if checks["external_uniform_total_moment_norm_Nm"] > 1.0e-8:
            raise AssertionError("external uniform pressure moment did not cancel")

    return {
        "kind": "netgen_vol_boundary_pressure_moment_validation",
        "validation_class": True,
        "force_learning": "boundary pressure rows reduce to force and pivot moments by triangle-centroid integration",
        "checks": checks,
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None, help="Optional external tri/tet Netgen .vol file")
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary(args.vol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    checks = summary["checks"]
    print("[boundary pressure force/moment]")
    print(f"  builtin uniform total force norm:  {checks['builtin_uniform_total_force_norm_N']:.3e} N")
    print(f"  builtin uniform total moment norm: {checks['builtin_uniform_total_moment_norm_Nm']:.3e} N m")
    print(f"  builtin zmax force:  {checks['builtin_zmax_force_N']} N")
    print(f"  builtin zmax moment: {checks['builtin_zmax_moment_Nm']} N m")
    if "external_uniform_total_force_norm_N" in checks:
        print(f"  external uniform total force norm:  {checks['external_uniform_total_force_norm_N']:.3e} N")
        print(f"  external uniform total moment norm: {checks['external_uniform_total_moment_norm_Nm']:.3e} N m")


if __name__ == "__main__":
    main()
