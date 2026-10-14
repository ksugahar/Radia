"""Validation-class enclosure/void region measurement through Cubit.

This example exercises the build123d helpers used before multi-region meshing:

    inner material regions -> enclosing box with bbox margin -> void/air region

The same STEP geometry is measured by headless Cubit so the build123d volume,
surface-area, and bounding-box measurements are cross-checked by an external
CAD kernel.

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_enclosure_cubit_measurement.py --require-cubit
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from build123d import Box, Pos  # noqa: E402

from radia_mcp.build123d.modeling import (  # noqa: E402
    enclosure_clearance_row,
    enclosure_difference_region,
    enclosing_box,
    shape_measurement_comparison_summary,
    shape_measurement_rows,
    tube,
)
from validation_build123d_cubit_measurement import (  # noqa: E402
    export_steps,
    find_cubit_bin,
    measure_steps_with_cubit,
)


OUT = HERE / "runs" / "validation_enclosure_cubit_measurement"
SUMMARY_JSON = HERE / "validation_enclosure_cubit_measurement_summary.json"


def _box_area(x: float, y: float, z: float) -> float:
    return 2.0 * (x * y + x * z + y * z)


def _tube_volume_area(r_in: float, r_out: float, h: float) -> tuple[float, float]:
    volume = math.pi * (r_out * r_out - r_in * r_in) * h
    area = 2.0 * math.pi * (r_out + r_in) * h
    area += 2.0 * math.pi * (r_out * r_out - r_in * r_in)
    return volume, area


def _rel_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def build_model() -> tuple[list[tuple[object, str]], dict, dict]:
    block_dims = (2.0, 3.0, 1.0)
    sleeve_dims = (0.35, 0.70, 1.20)
    margin = (0.90, 0.80, 0.70)

    block = (Pos(-1.60, 0.0, -0.40) * Box(*block_dims)).solid()
    block.label = "block"
    sleeve = (Pos(1.80, 0.0, 0.50) * tube(*sleeve_dims, label="sleeve")).solid()
    sleeve.label = "sleeve"

    enclosure = enclosing_box([block, sleeve], margin=margin, label="enclosure_box")
    void = enclosure_difference_region(enclosure, [block, sleeve], label="enclosure_void")
    clearance = enclosure_clearance_row(enclosure, [block, sleeve], name="enclosure_box")

    block_volume = block_dims[0] * block_dims[1] * block_dims[2]
    block_area = _box_area(*block_dims)
    sleeve_volume, sleeve_area = _tube_volume_area(*sleeve_dims)
    enclosure_size = clearance["inner_envelope"]["size"]
    enclosure_size = [enclosure_size[i] + 2.0 * margin[i] for i in range(3)]
    enclosure_volume = enclosure_size[0] * enclosure_size[1] * enclosure_size[2]
    enclosure_area = _box_area(*enclosure_size)

    analytic = {
        "block": {"volume": block_volume, "area": block_area},
        "sleeve": {"volume": sleeve_volume, "area": sleeve_area},
        "enclosure_box": {"volume": enclosure_volume, "area": enclosure_area},
        "enclosure_void": {
            "volume": enclosure_volume - block_volume - sleeve_volume,
            "area": enclosure_area + block_area + sleeve_area,
        },
    }
    cases = [
        (block, "block"),
        (sleeve, "sleeve"),
        (enclosure, "enclosure_box"),
        (void, "enclosure_void"),
    ]
    return cases, analytic, clearance


def analytic_comparison_rows(build_rows: list[dict], analytic: dict, rtol: float) -> list[dict]:
    rows = []
    for row in build_rows:
        ref = analytic[row["name"]]
        volume_rel_error = _rel_error(row["volume"], ref["volume"])
        area_rel_error = _rel_error(row["area"], ref["area"])
        rows.append({
            "name": row["name"],
            "reference_volume": ref["volume"],
            "reference_area": ref["area"],
            "build123d_volume": row["volume"],
            "build123d_area": row["area"],
            "volume_rel_error": volume_rel_error,
            "area_rel_error": area_rel_error,
            "rtol": rtol,
            "passed": volume_rel_error <= rtol and area_rel_error <= rtol,
        })
    return rows


def build_summary(
    out_dir: Path,
    cubit_bin_arg: str | None,
    require_cubit: bool,
    rtol: float,
    bbox_atol: float,
) -> dict:
    t0 = time.perf_counter()
    cases, analytic, clearance = build_model()
    build_rows = shape_measurement_rows(cases)
    step_paths = export_steps(cases, out_dir)

    cubit_bin = find_cubit_bin(cubit_bin_arg)
    if cubit_bin is None:
        if require_cubit:
            raise RuntimeError("Cubit was required but no Coreform Cubit bin directory was found")
        cubit = {"available": False, "rows": []}
    else:
        cubit = measure_steps_with_cubit(step_paths, cubit_bin)

    analytic_rows = analytic_comparison_rows(build_rows, analytic, rtol=rtol)
    cubit_comparison = shape_measurement_comparison_summary(
        build_rows,
        cubit.get("rows", []),
        rtol=rtol,
        measured_label="cubit",
        bbox_atol=bbox_atol,
    )
    cubit_by_name = {row["name"]: row for row in cubit.get("rows", [])}
    cubit_cmp_by_name = {row["name"]: row for row in cubit_comparison["rows"]}
    rows = []
    for row in build_rows:
        rows.append({
            "name": row["name"],
            "build123d": {
                "volume": row["volume"],
                "area": row["area"],
                "faces": row["faces"],
                "edges": row["edges"],
                "vertices": row["vertices"],
                "bbox_size": row["bounding_box"]["size"],
                "bounding_box": row["bounding_box"],
            },
            "analytic": analytic[row["name"]],
            "analytic_comparison": next(r for r in analytic_rows if r["name"] == row["name"]),
            "cubit": cubit_by_name.get(row["name"]),
            "cubit_comparison": cubit_cmp_by_name.get(row["name"]),
        })

    checks = {
        "rtol": rtol,
        "n_cases": len(rows),
        "n_analytic_passed": sum(1 for row in analytic_rows if row["passed"]),
        "n_cubit_passed": cubit_comparison["n_passed"],
        "max_analytic_volume_rel_error": max(row["volume_rel_error"] for row in analytic_rows),
        "max_analytic_area_rel_error": max(row["area_rel_error"] for row in analytic_rows),
        "max_cubit_volume_rel_error": cubit_comparison["max_volume_rel_error"],
        "max_cubit_area_rel_error": cubit_comparison["max_area_rel_error"],
        "max_cubit_bbox_abs_error": cubit_comparison["max_bbox_abs_error"],
        "n_cubit_bbox_compared": cubit_comparison["n_bbox_compared"],
        "cubit_available": bool(cubit.get("available")),
        "cubit_bin_name": cubit.get("bin_name"),
        "bbox_atol": bbox_atol,
        "min_bbox_clearance": clearance["min_clearance"],
        "nominal_void_volume": clearance["nominal_void_volume"],
        "contained_by_bbox": clearance["contained_by_bbox"],
    }
    checks["passed"] = (
        checks["n_analytic_passed"] == checks["n_cases"]
        and checks["n_cubit_passed"] == checks["n_cases"]
        and checks["contained_by_bbox"]
    )

    return {
        "kind": "build123d_enclosure_cubit_volume_area_bbox_cross_validation",
        "validation_class": True,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "checks": checks,
        "clearance": clearance,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--cubit-bin")
    parser.add_argument("--require-cubit", action="store_true")
    parser.add_argument("--rtol", type=float, default=1.0e-5)
    parser.add_argument("--bbox-atol", type=float, default=1.0e-6)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = build_summary(
        args.out_dir,
        args.cubit_bin,
        args.require_cubit,
        args.rtol,
        args.bbox_atol,
    )
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[build123d enclosure -> Cubit volume/area/bbox cross validation]")
    print(
        f"  cases={checks['n_cases']} analytic={checks['n_analytic_passed']} "
        f"cubit={checks['n_cubit_passed']} Cubit available={checks['cubit_available']} "
        f"({checks['cubit_bin_name']})"
    )
    print(
        f"  max analytic rel errors: volume={checks['max_analytic_volume_rel_error']:.3e}, "
        f"area={checks['max_analytic_area_rel_error']:.3e}"
    )
    print(
        f"  max Cubit rel errors: volume={checks['max_cubit_volume_rel_error']:.3e}, "
        f"area={checks['max_cubit_area_rel_error']:.3e}"
    )
    print(
        f"  max Cubit bbox abs error={checks['max_cubit_bbox_abs_error']:.3e} "
        f"({checks['n_cubit_bbox_compared']} bbox comparisons)"
    )
    print(
        f"  min bbox clearance={checks['min_bbox_clearance']:.6g}, "
        f"nominal void volume={checks['nominal_void_volume']:.12g}"
    )
    print(f"[OK] wrote {args.summary}")
    return 0 if checks["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
