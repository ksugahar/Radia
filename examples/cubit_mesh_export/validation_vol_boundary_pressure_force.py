"""Validation-class `.vol` boundary pressure-to-force rows.

Once a Cubit/Coreform sideset has a reliable oriented vector area, a scalar
pressure load is just

    F_boundary = pressure * vector_area.

This example checks that named Netgen `.vol` boundaries can be turned into
force rows for pressure loads, Maxwell pressure, or acoustic pressure examples.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_pressure_force.py
    python examples/cubit_mesh_export/validation_vol_boundary_pressure_force.py --vol C:\\temp\\box.vol
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402
from validation_vol_boundary_normal_vectors import BOX_SIX_BOUNDARY_VOL  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_pressure_force_summary.json")


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _sum_forces(rows):
    return [
        sum(row["force_N"][axis] for row in rows)
        for axis in range(3)
    ]


def _record(label, mesh):
    uniform_pressure = 2.0
    uniform_rows = list(mesh.boundary_pressure_force_rows({
        row["name"]: uniform_pressure for row in mesh.boundary_normal_summary_rows()
    }))
    zmax_rows = list(mesh.boundary_pressure_force_rows({"zmax": 2.0}, default_pressure=0.0))
    by_name = {row["name"]: row for row in zmax_rows}
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "uniform_pressure_Pa": uniform_pressure,
        "uniform_pressure_rows": uniform_rows,
        "uniform_total_force_N": _sum_forces(uniform_rows),
        "zmax_only_rows": zmax_rows,
        "zmax_only_total_force_N": _sum_forces(zmax_rows),
        "zmax_force_N": by_name.get("zmax", {}).get("force_N"),
    }


def build_summary(external_vol=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if external_vol is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(external_vol)))

    builtin = records[0]
    checks = {
        "builtin_uniform_total_force_norm_N": _norm(builtin["uniform_total_force_N"]),
        "builtin_zmax_force_N": builtin["zmax_force_N"],
        "builtin_zmax_force_abs_error_N": _norm([
            builtin["zmax_force_N"][0],
            builtin["zmax_force_N"][1],
            builtin["zmax_force_N"][2] - 12.0,
        ]),
    }
    if checks["builtin_uniform_total_force_norm_N"] > 1.0e-14:
        raise AssertionError("uniform pressure on closed box should sum to zero")
    if checks["builtin_zmax_force_abs_error_N"] > 1.0e-14:
        raise AssertionError("zmax pressure force drifted")

    if external_vol is not None:
        external = records[-1]
        checks.update({
            "external_uniform_total_force_norm_N": _norm(external["uniform_total_force_N"]),
            "external_zmax_force_N": external["zmax_force_N"],
            "external_zmax_force_abs_error_N": _norm([
                external["zmax_force_N"][0],
                external["zmax_force_N"][1],
                external["zmax_force_N"][2] - 12.0,
            ]) if external["zmax_force_N"] is not None else math.inf,
        })
        if checks["external_uniform_total_force_norm_N"] > 1.0e-8:
            raise AssertionError("external uniform pressure force did not cancel")
        if checks["external_zmax_force_abs_error_N"] > 1.0e-8:
            raise AssertionError("external zmax pressure force drifted")

    return {
        "kind": "netgen_vol_boundary_pressure_force_validation",
        "validation_class": True,
        "force_learning": "pressure loads use force = pressure * oriented boundary vector area",
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
    print("[boundary pressure force]")
    print(f"  builtin uniform total force norm: {checks['builtin_uniform_total_force_norm_N']:.3e} N")
    print(f"  builtin zmax force: {checks['builtin_zmax_force_N']} N")
    if "external_uniform_total_force_norm_N" in checks:
        print(f"  external uniform total force norm: {checks['external_uniform_total_force_norm_N']:.3e} N")
        print(f"  external zmax force: {checks['external_zmax_force_N']} N")


if __name__ == "__main__":
    main()
