"""Validation-class build123d analytic box pressure resultant.

The CAD-side reference starts from a real build123d ``Box`` measurement, then
reduces constant pressure on its planar faces to net force/moment metrics.  The
same pressure table is compared against a named Netgen `.vol` boundary mesh so
the build123d design side and Coreform/Cubit mesh side share one force summary
vocabulary.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_resultant.py
    python examples/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_resultant.py --vol C:\\temp\\box.vol
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from build123d import Box, Pos


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
CUBIT_EXAMPLES = REPO / "examples" / "cubit_mesh_export"
for path in (SRC, CUBIT_EXAMPLES):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from radia_mcp.build123d.modeling import box_face_pressure_resultant_summary, shape_measurement_row  # noqa: E402
from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402
from validation_vol_boundary_normal_vectors import BOX_SIX_BOUNDARY_VOL  # noqa: E402


SUMMARY_JSON = HERE / "validation_build123d_cubit_pressure_resultant_summary.json"
SIZE = (2.0, 3.0, 5.0)
CENTER = (1.0, 1.5, 2.5)
PRESSURE = {"zmax": 2.0}
UNIFORM_PRESSURE = 2.0


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _vector_error(a, b):
    return _norm([a[axis] - b[axis] for axis in range(3)])


def _build123d_box_reference():
    box = (Pos(*CENTER) * Box(*SIZE)).solid()
    measurement = shape_measurement_row(box, name="build123d_box")
    size = measurement["bounding_box"]["size"]
    center = measurement["bounding_box"]["center"]
    return measurement, size, center


def _summary_pair(mesh, size, center):
    reference_uniform = box_face_pressure_resultant_summary(
        size,
        {},
        center=center,
        default_pressure=UNIFORM_PRESSURE,
    )
    reference_zmax = box_face_pressure_resultant_summary(
        size,
        PRESSURE,
        center=center,
        default_pressure=0.0,
    )
    mesh_uniform = mesh.boundary_pressure_resultant_summary(
        {},
        default_pressure=UNIFORM_PRESSURE,
    )
    mesh_zmax = mesh.boundary_pressure_resultant_summary(
        PRESSURE,
        default_pressure=0.0,
    )
    return {
        "reference_uniform": reference_uniform,
        "reference_zmax": reference_zmax,
        "mesh_uniform": mesh_uniform,
        "mesh_zmax": mesh_zmax,
        "force_error_uniform_N": _vector_error(
            reference_uniform["total_force_N"],
            mesh_uniform["total_force_N"],
        ),
        "moment_error_uniform_Nm": _vector_error(
            reference_uniform["total_moment_about_pivot_Nm"],
            mesh_uniform["total_moment_about_pivot_Nm"],
        ),
        "force_error_zmax_N": _vector_error(
            reference_zmax["total_force_N"],
            mesh_zmax["total_force_N"],
        ),
        "moment_error_zmax_Nm": _vector_error(
            reference_zmax["total_moment_about_pivot_Nm"],
            mesh_zmax["total_moment_about_pivot_Nm"],
        ),
    }


def build_summary(vol_path: Path | None = None) -> dict:
    measurement, size, center = _build123d_box_reference()
    builtin_mesh = parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL)
    builtin = _summary_pair(builtin_mesh, size, center)
    checks = {
        "build123d_volume": measurement["volume"],
        "build123d_area": measurement["area"],
        "builtin_uniform_force_norm_N": builtin["reference_uniform"]["total_force_magnitude_N"],
        "builtin_uniform_moment_norm_Nm": builtin["reference_uniform"]["total_moment_magnitude_Nm"],
        "builtin_uniform_abs_force_sum_N": builtin["reference_uniform"]["absolute_force_sum_N"],
        "builtin_mesh_force_error_uniform_N": builtin["force_error_uniform_N"],
        "builtin_mesh_moment_error_uniform_Nm": builtin["moment_error_uniform_Nm"],
        "builtin_zmax_force_N": builtin["reference_zmax"]["total_force_N"],
        "builtin_zmax_moment_Nm": builtin["reference_zmax"]["total_moment_about_pivot_Nm"],
        "builtin_mesh_force_error_zmax_N": builtin["force_error_zmax_N"],
        "builtin_mesh_moment_error_zmax_Nm": builtin["moment_error_zmax_Nm"],
    }
    if abs(checks["build123d_volume"] - 30.0) > 1.0e-12:
        raise AssertionError("build123d box volume drifted")
    if abs(checks["build123d_area"] - 62.0) > 1.0e-12:
        raise AssertionError("build123d box area drifted")
    if checks["builtin_uniform_force_norm_N"] > 1.0e-14:
        raise AssertionError("uniform pressure should have zero net force")
    if checks["builtin_uniform_moment_norm_Nm"] > 1.0e-14:
        raise AssertionError("uniform pressure should have zero net moment")
    if abs(checks["builtin_uniform_abs_force_sum_N"] - 124.0) > 1.0e-14:
        raise AssertionError("uniform absolute force sum drifted")
    if _vector_error(checks["builtin_zmax_force_N"], (0.0, 0.0, 12.0)) > 1.0e-14:
        raise AssertionError("zmax pressure force drifted")
    if _vector_error(checks["builtin_zmax_moment_Nm"], (18.0, -12.0, 0.0)) > 1.0e-14:
        raise AssertionError("zmax pressure moment drifted")
    for key in (
        "builtin_mesh_force_error_uniform_N",
        "builtin_mesh_moment_error_uniform_Nm",
        "builtin_mesh_force_error_zmax_N",
        "builtin_mesh_moment_error_zmax_Nm",
    ):
        if checks[key] > 1.0e-14:
            raise AssertionError(f"build123d/.vol pressure summary mismatch: {key}")

    external = None
    if vol_path is not None:
        external_mesh = read_netgen_tri_tet_vol(vol_path)
        external = _summary_pair(external_mesh, size, center)
        checks.update({
            "external_mesh_force_error_uniform_N": external["force_error_uniform_N"],
            "external_mesh_moment_error_uniform_Nm": external["moment_error_uniform_Nm"],
            "external_mesh_force_error_zmax_N": external["force_error_zmax_N"],
            "external_mesh_moment_error_zmax_Nm": external["moment_error_zmax_Nm"],
        })
        for key in (
            "external_mesh_force_error_uniform_N",
            "external_mesh_moment_error_uniform_Nm",
            "external_mesh_force_error_zmax_N",
            "external_mesh_moment_error_zmax_Nm",
        ):
            if checks[key] > 1.0e-8:
                raise AssertionError(f"external .vol pressure summary mismatch: {key}")

    return {
        "kind": "build123d_cubit_pressure_resultant_validation",
        "validation_class": True,
        "size": size,
        "center": center,
        "pressure_by_face_Pa": PRESSURE,
        "uniform_pressure_Pa": UNIFORM_PRESSURE,
        "build123d_measurement": measurement,
        "checks": checks,
        "builtin_vol_pair": builtin,
        "external_vol_pair": external,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=SUMMARY_JSON)
    args = parser.parse_args()

    summary = build_summary(args.vol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    checks = summary["checks"]
    print("[build123d pressure resultant]")
    print(f"  build123d volume/area: {checks['build123d_volume']} / {checks['build123d_area']}")
    print(f"  uniform force norm:    {checks['builtin_uniform_force_norm_N']:.3e} N")
    print(f"  uniform moment norm:   {checks['builtin_uniform_moment_norm_Nm']:.3e} N m")
    print(f"  zmax force:            {checks['builtin_zmax_force_N']} N")
    print(f"  zmax moment:           {checks['builtin_zmax_moment_Nm']} N m")
    print(f"  builtin .vol force err: {checks['builtin_mesh_force_error_zmax_N']:.3e} N")
    print(f"  builtin .vol moment err: {checks['builtin_mesh_moment_error_zmax_Nm']:.3e} N m")
    if "external_mesh_force_error_zmax_N" in checks:
        print(f"  external .vol force err: {checks['external_mesh_force_error_zmax_N']:.3e} N")
        print(f"  external .vol moment err: {checks['external_mesh_moment_error_zmax_Nm']:.3e} N m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
