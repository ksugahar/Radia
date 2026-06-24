"""Validation-class `.vol` boundary vector-traction force/moment rows.

Named Cubit/Coreform sidesets exported as Netgen `.vol` provide boundary
triangles and names.  A constant vector traction on a boundary has

    dF = t dS,
    dM = (r - pivot) x dF.

This example checks the per-boundary rows and the generic force/moment reducer.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_traction_moment.py
    python examples/cubit_mesh_export/validation_vol_boundary_traction_moment.py --vol C:\\temp\\box.vol
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


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_traction_moment_summary.json")
TRACTION = (1.0, -2.0, 3.0)


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


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
    rows = list(mesh.boundary_traction_force_moment_rows(
        {"zmax": TRACTION},
        default_traction=(0.0, 0.0, 0.0),
    ))
    shifted_rows = list(mesh.boundary_traction_force_moment_rows(
        {"zmax": TRACTION},
        default_traction=(0.0, 0.0, 0.0),
        pivot_m=(1.0, 1.5, 5.0),
    ))
    zmax = {row["name"]: row for row in rows}.get("zmax")
    shifted = {row["name"]: row for row in shifted_rows}.get("zmax")
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "traction_N_per_m2": TRACTION,
        "rows": rows,
        "zmax_row": zmax,
        "zmax_resultant": _resultant_from_rows(rows),
        "zmax_shifted_pivot_row": shifted,
        "zmax_shifted_resultant": _resultant_from_rows(shifted_rows, pivot=(1.0, 1.5, 5.0)),
    }


def _row_resultant_errors(record):
    zmax = record["zmax_row"]
    resultant = record["zmax_resultant"]
    return {
        "force_error_N": _norm([
            resultant["total_force"][axis] - zmax["force_N"][axis]
            for axis in range(3)
        ]),
        "moment_error_Nm": _norm([
            resultant["total_moment"][axis] - zmax["moment_about_pivot_Nm"][axis]
            for axis in range(3)
        ]),
    }


def build_summary(external_vol=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if external_vol is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(external_vol)))

    builtin = records[0]
    zmax = builtin["zmax_row"]
    shifted = builtin["zmax_shifted_pivot_row"]
    errors = _row_resultant_errors(builtin)
    checks = {
        "builtin_zmax_centroid_m": zmax["centroid_m"],
        "builtin_zmax_force_N": zmax["force_N"],
        "builtin_zmax_moment_Nm": zmax["moment_about_pivot_Nm"],
        "builtin_zmax_shifted_moment_Nm": shifted["moment_about_pivot_Nm"],
        "builtin_zmax_resultant_force_error_N": errors["force_error_N"],
        "builtin_zmax_resultant_moment_error_Nm": errors["moment_error_Nm"],
    }
    if _norm([zmax["force_N"][0] - 6.0, zmax["force_N"][1] + 12.0, zmax["force_N"][2] - 18.0]) > 1.0e-14:
        raise AssertionError("zmax vector traction force drifted")
    if _norm([
        zmax["moment_about_pivot_Nm"][0] - 87.0,
        zmax["moment_about_pivot_Nm"][1] - 12.0,
        zmax["moment_about_pivot_Nm"][2] + 21.0,
    ]) > 1.0e-14:
        raise AssertionError("zmax vector traction moment drifted")
    if _norm(shifted["moment_about_pivot_Nm"]) > 1.0e-14:
        raise AssertionError("zmax shifted-pivot moment should vanish")
    if errors["force_error_N"] > 1.0e-14 or errors["moment_error_Nm"] > 1.0e-14:
        raise AssertionError("generic resultant drifted")

    if external_vol is not None:
        external = records[-1]
        external_errors = _row_resultant_errors(external)
        checks.update({
            "external_zmax_force_N": (
                external["zmax_row"].get("force_N") if external["zmax_row"] else None
            ),
            "external_zmax_moment_Nm": (
                external["zmax_row"].get("moment_about_pivot_Nm") if external["zmax_row"] else None
            ),
            "external_zmax_resultant_force_error_N": external_errors["force_error_N"],
            "external_zmax_resultant_moment_error_Nm": external_errors["moment_error_Nm"],
        })
        if external_errors["force_error_N"] > 1.0e-8:
            raise AssertionError("external resultant force drifted")
        if external_errors["moment_error_Nm"] > 1.0e-8:
            raise AssertionError("external resultant moment drifted")

    return {
        "kind": "netgen_vol_boundary_traction_moment_validation",
        "validation_class": True,
        "force_learning": "boundary vector traction rows reduce to force and pivot moments by triangle-centroid integration",
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
    print("[boundary vector traction force/moment]")
    print(f"  builtin zmax force:  {checks['builtin_zmax_force_N']} N")
    print(f"  builtin zmax moment: {checks['builtin_zmax_moment_Nm']} N m")
    print(f"  builtin shifted moment: {checks['builtin_zmax_shifted_moment_Nm']} N m")
    if "external_zmax_resultant_force_error_N" in checks:
        print(f"  external zmax force: {checks['external_zmax_force_N']} N")
        print(f"  external resultant force error: {checks['external_zmax_resultant_force_error_N']:.3e} N")
        print(f"  external resultant moment error: {checks['external_zmax_resultant_moment_error_Nm']:.3e} N m")


if __name__ == "__main__":
    main()
