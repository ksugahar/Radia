"""Validation-class `.vol` boundary pressure resultant summary.

Coreform/Cubit sidesets exported through Netgen `.vol` give oriented surface
triangles.  This example reduces constant scalar pressure rows to the net
force/moment summary that downstream electromagnetic and acoustic force gates
can consume directly.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_pressure_resultant.py
    python examples/cubit_mesh_export/validation_vol_boundary_pressure_resultant.py --vol C:\\temp\\box.vol
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


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_pressure_resultant_summary.json")


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _record(label, mesh):
    uniform_pressure = 2.0
    uniform = mesh.boundary_pressure_resultant_summary({}, default_pressure=uniform_pressure)
    record = {
        "label": label,
        "mesh_summary": mesh.summary(),
        "uniform_pressure_Pa": uniform_pressure,
        "uniform_resultant": uniform,
    }
    if "zmax" in set(mesh.boundary_names.values()):
        record["zmax_resultant"] = mesh.boundary_pressure_resultant_summary(
            {"zmax": 2.0},
            default_pressure=0.0,
        )
        record["zmax_shifted_pivot_resultant"] = mesh.boundary_pressure_resultant_summary(
            {"zmax": 2.0},
            default_pressure=0.0,
            pivot_m=(1.0, 1.5, 0.0),
        )
    return record


def build_summary(external_vol=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if external_vol is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(external_vol)))

    builtin = records[0]
    uniform = builtin["uniform_resultant"]
    zmax = builtin["zmax_resultant"]
    shifted = builtin["zmax_shifted_pivot_resultant"]
    checks = {
        "builtin_uniform_total_force_norm_N": uniform["total_force_magnitude_N"],
        "builtin_uniform_total_moment_norm_Nm": uniform["total_moment_magnitude_Nm"],
        "builtin_uniform_force_balance_ratio": uniform["force_balance_ratio"],
        "builtin_uniform_moment_balance_ratio": uniform["moment_balance_ratio"],
        "builtin_uniform_absolute_force_sum_N": uniform["absolute_force_sum_N"],
        "builtin_surface_vector_area_norm": uniform["surface_vector_area_norm"],
        "builtin_zmax_force_N": zmax["total_force_N"],
        "builtin_zmax_moment_Nm": zmax["total_moment_about_pivot_Nm"],
        "builtin_zmax_shifted_moment_Nm": shifted["total_moment_about_pivot_Nm"],
    }
    if checks["builtin_uniform_total_force_norm_N"] > 1.0e-14:
        raise AssertionError("uniform closed-box pressure should have zero net force")
    if checks["builtin_uniform_total_moment_norm_Nm"] > 1.0e-14:
        raise AssertionError("uniform closed-box pressure should have zero net moment")
    if checks["builtin_uniform_absolute_force_sum_N"] != 124.0:
        raise AssertionError("box absolute force sum drifted")
    if checks["builtin_surface_vector_area_norm"] > 1.0e-14:
        raise AssertionError("box surface vector area should close")
    if _norm([zmax["total_force_N"][0], zmax["total_force_N"][1], zmax["total_force_N"][2] - 12.0]) > 1.0e-14:
        raise AssertionError("zmax resultant force drifted")
    if _norm([
        zmax["total_moment_about_pivot_Nm"][0] - 18.0,
        zmax["total_moment_about_pivot_Nm"][1] + 12.0,
        zmax["total_moment_about_pivot_Nm"][2],
    ]) > 1.0e-14:
        raise AssertionError("zmax resultant moment drifted")
    if _norm(shifted["total_moment_about_pivot_Nm"]) > 1.0e-14:
        raise AssertionError("shifted pivot should remove the zmax resultant moment")

    if external_vol is not None:
        external = records[-1]["uniform_resultant"]
        checks.update({
            "external_uniform_total_force_norm_N": external["total_force_magnitude_N"],
            "external_uniform_total_moment_norm_Nm": external["total_moment_magnitude_Nm"],
            "external_uniform_force_balance_ratio": external["force_balance_ratio"],
            "external_surface_vector_area_norm": external["surface_vector_area_norm"],
        })
        if checks["external_uniform_total_force_norm_N"] > 1.0e-8:
            raise AssertionError("external uniform pressure force did not cancel")
        if checks["external_uniform_total_moment_norm_Nm"] > 1.0e-8:
            raise AssertionError("external uniform pressure moment did not cancel")

    return {
        "kind": "netgen_vol_boundary_pressure_resultant_validation",
        "validation_class": True,
        "force_learning": "closed .vol surfaces cancel uniform pressure while one-sided pressure reduces to p times oriented vector area",
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
    print("[boundary pressure resultant]")
    print(f"  builtin uniform total force norm:  {checks['builtin_uniform_total_force_norm_N']:.3e} N")
    print(f"  builtin uniform total moment norm: {checks['builtin_uniform_total_moment_norm_Nm']:.3e} N m")
    print(f"  builtin zmax force:  {checks['builtin_zmax_force_N']} N")
    print(f"  builtin zmax moment: {checks['builtin_zmax_moment_Nm']} N m")
    if "external_uniform_total_force_norm_N" in checks:
        print(f"  external uniform total force norm:  {checks['external_uniform_total_force_norm_N']:.3e} N")
        print(f"  external uniform total moment norm: {checks['external_uniform_total_moment_norm_Nm']:.3e} N m")


if __name__ == "__main__":
    main()
