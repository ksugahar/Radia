"""Validation-class build123d measurement health cross-check.

This example keeps the volume/area cross validation readable at assembly scale:

    labelled build123d shapes -> STEP -> headless Cubit measurement

The build123d side reports region volume fractions and bounding-box fill
fraction.  The Cubit side checks volume, surface area, and bounding boxes, then
keeps only the worst comparison rows in the health report.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_build123d_measurement_health.py --require-cubit
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
    annular_segment,
    assembly,
    racetrack_coil,
    shape_measurement_health_summary,
    shape_measurement_inventory_summary,
    shape_measurement_rows,
    tube,
)
from validation_build123d_cubit_measurement import (  # noqa: E402
    export_steps,
    find_cubit_bin,
    measure_steps_with_cubit,
)


OUT = HERE / "runs" / "validation_build123d_measurement_health"
SUMMARY_JSON = HERE / "validation_build123d_measurement_health_summary.json"


def _box_area(x: float, y: float, z: float) -> float:
    return 2.0 * (x * y + x * z + y * z)


def _tube_volume_area(r_in: float, r_out: float, h: float) -> tuple[float, float]:
    return (
        math.pi * (r_out * r_out - r_in * r_in) * h,
        2.0 * math.pi * (r_out + r_in) * h
        + 2.0 * math.pi * (r_out * r_out - r_in * r_in),
    )


def _rel_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def build_cases() -> list[tuple[object, str]]:
    core = (Pos(-3.4, 0.0, 0.0) * Box(2.5, 1.5, 1.0)).solid()
    core.label = "core_block"

    sleeve = (Pos(-0.7, 0.0, 0.1) * tube(0.20, 0.45, 1.6, label="cooling_sleeve")).solid()
    sleeve.label = "cooling_sleeve"

    sector = (Pos(1.6, 0.0, 0.0) * annular_segment(0.8, 1.25, 0.7, 10.0, 80.0, label="sector_pole")).solid()
    sector.label = "sector_pole"

    coil = (Pos(4.2, 0.0, 0.0) * racetrack_coil(
        length=3.4,
        width=2.2,
        band=0.28,
        h=0.55,
        corner_radius=0.55,
        label="racetrack_coil",
    )).solid()
    coil.label = "racetrack_coil"

    return [
        (core, core.label),
        (sleeve, sleeve.label),
        (sector, sector.label),
        (coil, coil.label),
    ]


def analytic_reference(name: str) -> dict[str, float] | None:
    if name == "core_block":
        x, y, z = 2.5, 1.5, 1.0
        return {"volume": x * y * z, "area": _box_area(x, y, z)}
    if name == "cooling_sleeve":
        volume, area = _tube_volume_area(0.20, 0.45, 1.6)
        return {"volume": volume, "area": area}
    if name == "sector_pole":
        volume = (70.0 / 360.0) * math.pi * (1.25 * 1.25 - 0.8 * 0.8) * 0.7
        return {"volume": volume}
    return None


def analytic_comparison_rows(build_rows: list[dict]) -> list[dict[str, object]]:
    rows = []
    for row in build_rows:
        ref = analytic_reference(row["name"])
        if ref is None:
            continue
        record: dict[str, object] = {
            "name": row["name"],
            "reference_volume": ref["volume"],
            "build123d_volume": row["volume"],
            "volume_rel_error": _rel_error(row["volume"], ref["volume"]),
        }
        if "area" in ref:
            record.update({
                "reference_area": ref["area"],
                "build123d_area": row["area"],
                "area_rel_error": _rel_error(row["area"], ref["area"]),
            })
        rows.append(record)
    return rows


def build_summary(
    out_dir: Path,
    cubit_bin_arg: str | None,
    require_cubit: bool,
    rtol: float,
    bbox_atol: float,
) -> dict[str, object]:
    t0 = time.perf_counter()
    cases = build_cases()
    build_rows = shape_measurement_rows(assembly(*[shape for shape, _name in cases], label="measurement_health_demo"))
    step_paths = export_steps(cases, out_dir)

    cubit_bin = find_cubit_bin(cubit_bin_arg)
    if cubit_bin is None:
        if require_cubit:
            raise RuntimeError("Cubit was required but no Coreform Cubit bin directory was found")
        cubit = {"available": False, "rows": []}
    else:
        cubit = measure_steps_with_cubit(step_paths, cubit_bin)

    health = shape_measurement_health_summary(
        build_rows,
        cubit.get("rows", []),
        rtol=rtol,
        measured_label="cubit",
        bbox_atol=bbox_atol,
        worst_limit=4,
    )
    inventory = shape_measurement_inventory_summary(build_rows)
    analytic_rows = analytic_comparison_rows(build_rows)

    checks = {
        "rtol": rtol,
        "bbox_atol": bbox_atol,
        "n_cases": len(build_rows),
        "cubit_available": bool(cubit.get("available")),
        "cubit_bin_name": cubit.get("bin_name"),
        "health_status": health["status"],
        "n_cubit_passed": health["comparison_summary"]["n_passed"],
        "max_cubit_volume_rel_error": health["comparison_summary"]["max_volume_rel_error"],
        "max_cubit_area_rel_error": health["comparison_summary"]["max_area_rel_error"],
        "max_cubit_bbox_abs_error": health["comparison_summary"]["max_bbox_abs_error"],
        "bbox_fill_fraction": inventory["bbox_fill_fraction"],
        "total_volume": inventory["total_volume"],
        "n_analytic_rows": len(analytic_rows),
        "max_analytic_volume_rel_error": max((row["volume_rel_error"] for row in analytic_rows), default=0.0),
        "max_analytic_area_rel_error": max((row.get("area_rel_error", 0.0) for row in analytic_rows), default=0.0),
    }
    checks["passed"] = bool(health["ok_for_geometry_roundtrip"])

    return {
        "kind": "build123d_measurement_health_cross_validation",
        "validation_class": True,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "checks": checks,
        "inventory": inventory,
        "analytic_comparison_rows": analytic_rows,
        "health": health,
        "build123d_rows": build_rows,
        "cubit_rows": cubit.get("rows", []),
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
    print("[build123d measurement health -> Cubit volume/area/bbox]")
    print(
        f"  cases={checks['n_cases']} cubit={checks['n_cubit_passed']} "
        f"available={checks['cubit_available']} ({checks['cubit_bin_name']}) "
        f"status={checks['health_status']}"
    )
    print(
        f"  max Cubit rel errors: volume={checks['max_cubit_volume_rel_error']:.3e}, "
        f"area={checks['max_cubit_area_rel_error']:.3e}; "
        f"bbox abs={checks['max_cubit_bbox_abs_error']:.3e}"
    )
    print(
        f"  assembly total volume={checks['total_volume']:.12g}, "
        f"bbox fill fraction={checks['bbox_fill_fraction']:.6g}"
    )
    print(f"[OK] wrote {args.summary}")
    return 0 if checks["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
