"""Validation-class build123d -> Cubit measurement cross-check.

The example strengthens the build123d API by making volume/area/bbox
measurements explicit and JSON-friendly, then checks the same STEP geometry
with Cubit's headless Python API:

    build123d shape -> STEP -> Cubit import -> volume/area/bbox

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_measurement.py --require-cubit

Without ``--require-cubit`` the script still writes the build123d measurement
summary and marks the Cubit stage as skipped when Cubit is not installed.  The
default tolerance is ``1e-5`` because this is a cross-kernel STEP round trip,
not two calls into the same OCCT object.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from build123d import Box, export_step  # noqa: E402

from radia_mcp.build123d.modeling import (  # noqa: E402
    racetrack_coil,
    shape_measurement_comparison_summary,
    shape_measurement_rows,
    tube,
)


OUT = HERE / "runs" / "validation_build123d_cubit_measurement"
SUMMARY_JSON = HERE / "validation_build123d_cubit_measurement_summary.json"


def build_cases() -> list[tuple[object, str]]:
    box = Box(2.0, 3.0, 4.0).solid()
    box.label = "box_2x3x4"

    annulus = tube(0.6, 1.0, 2.5, label="tube_r0p6_1p0_h2p5")

    coil = racetrack_coil(
        length=12.0,
        width=8.0,
        band=1.2,
        h=2.0,
        corner_radius=2.0,
        label="racetrack_12x8_band1p2",
    )
    return [(box, box.label), (annulus, annulus.label), (coil, coil.label)]


def analytic_reference(name: str) -> dict[str, float] | None:
    if name == "box_2x3x4":
        return {"volume": 24.0, "area": 52.0}
    if name == "tube_r0p6_1p0_h2p5":
        r_in = 0.6
        r_out = 1.0
        h = 2.5
        return {
            "volume": math.pi * (r_out * r_out - r_in * r_in) * h,
            "area": 2.0 * math.pi * (r_out + r_in) * h
            + 2.0 * math.pi * (r_out * r_out - r_in * r_in),
        }
    return None


def export_steps(cases: list[tuple[object, str]], out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for shape, name in cases:
        path = out_dir / f"{name}.step"
        export_step(shape, str(path))
        paths[name] = path
    return paths


def find_cubit_bin(explicit: str | None = None) -> Path | None:
    if explicit:
        path = Path(explicit)
        return path if (path / "cubit.py").is_file() else None
    env = os.environ.get("CUBIT_BIN") or os.environ.get("CORE_FORM_CUBIT_BIN")
    if env:
        path = Path(env)
        return path if (path / "cubit.py").is_file() else None
    root = Path("C:/Program Files")
    candidates = sorted(root.glob("Coreform Cubit */bin/cubit.py"))
    return candidates[-1].parent if candidates else None


def _cubit_volume_bounding_box(cubit, volume_ids: list[int]) -> dict[str, object] | None:
    if not volume_ids:
        return None
    boxes = [list(cubit.volume(vid).bounding_box()) for vid in volume_ids]
    bbox_min = [min(box[i] for box in boxes) for i in range(3)]
    bbox_max = [max(box[i + 3] for box in boxes) for i in range(3)]
    bbox_size = [hi - lo for lo, hi in zip(bbox_min, bbox_max)]
    bbox_center = [(lo + hi) / 2.0 for lo, hi in zip(bbox_min, bbox_max)]
    return {
        "min": [float(value) for value in bbox_min],
        "max": [float(value) for value in bbox_max],
        "center": [float(value) for value in bbox_center],
        "size": [float(value) for value in bbox_size],
        "diagonal": math.sqrt(sum(value * value for value in bbox_size)),
    }


def measure_steps_with_cubit(step_paths: dict[str, Path], cubit_bin: Path) -> dict[str, object]:
    sys.path.insert(0, str(cubit_bin))
    os.environ["PATH"] = str(cubit_bin) + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(cubit_bin))
        python3 = cubit_bin / "python3"
        if python3.is_dir():
            os.add_dll_directory(str(python3))

    import cubit  # noqa: PLC0415

    plugin_dir = cubit_bin / "plugins"
    cubit.init([
        "cubit",
        "-nojournal",
        "-batch",
        "-nographics",
        "-commandplugindir",
        str(plugin_dir),
    ])

    rows = []
    for index, (name, step_path) in enumerate(step_paths.items(), start=1):
        cubit.cmd("reset")
        cubit.cmd(f'import step "{step_path.as_posix()}"')
        volume_ids = list(cubit.get_entities("volume"))
        surface_ids = list(cubit.get_entities("surface"))
        rows.append({
            "index": index,
            "name": name,
            "volumes": len(volume_ids),
            "surfaces": len(surface_ids),
            "volume": sum(cubit.volume(vid).volume() for vid in volume_ids),
            "area": sum(cubit.surface(sid).area() for sid in surface_ids),
            "bounding_box": _cubit_volume_bounding_box(cubit, volume_ids),
        })
    cubit.cmd("reset")
    return {
        "available": True,
        "bin_name": cubit_bin.parent.name,
        "rows": rows,
    }


def _rel_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def build_summary(
    out_dir: Path,
    cubit_bin_arg: str | None,
    require_cubit: bool,
    rtol: float,
    bbox_atol: float,
) -> dict:
    t0 = time.perf_counter()
    cases = build_cases()
    build_rows = shape_measurement_rows(cases)
    step_paths = export_steps(cases, out_dir)

    cubit_bin = find_cubit_bin(cubit_bin_arg)
    if cubit_bin is None:
        if require_cubit:
            raise RuntimeError("Cubit was required but no Coreform Cubit bin directory was found")
        cubit = {"available": False, "rows": []}
    else:
        cubit = measure_steps_with_cubit(step_paths, cubit_bin)

    cubit_comparison = shape_measurement_comparison_summary(
        build_rows,
        cubit.get("rows", []),
        rtol=rtol,
        measured_label="cubit",
        bbox_atol=bbox_atol,
    )
    cubit_by_name = {row["name"]: row for row in cubit.get("rows", [])}
    comparison_by_name = {row["name"]: row for row in cubit_comparison["rows"]}
    rows = []
    for row in build_rows:
        ref = analytic_reference(row["name"])
        cubit_row = cubit_by_name.get(row["name"])
        comparison = comparison_by_name.get(row["name"])
        record = {
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
            "analytic": ref,
            "analytic_volume_rel_error": (
                _rel_error(row["volume"], ref["volume"]) if ref else None
            ),
            "analytic_area_rel_error": (
                _rel_error(row["area"], ref["area"]) if ref else None
            ),
            "cubit": None,
            "measurement_comparison": comparison,
            "cubit_volume_rel_error": None,
            "cubit_area_rel_error": None,
            "passed": False,
        }
        if cubit_row is not None:
            record["cubit"] = {
                "volume": cubit_row["volume"],
                "area": cubit_row["area"],
                "volumes": cubit_row["volumes"],
                "surfaces": cubit_row["surfaces"],
                "bounding_box": cubit_row["bounding_box"],
            }
            record["cubit_volume_rel_error"] = comparison["volume_rel_error"]
            record["cubit_area_rel_error"] = comparison["area_rel_error"]

        analytic_ok = (
            ref is None
            or (
                record["analytic_volume_rel_error"] <= rtol
                and record["analytic_area_rel_error"] <= rtol
            )
        )
        cubit_ok = comparison is not None and comparison["passed"]
        record["passed"] = bool(row["is_valid"] and analytic_ok and cubit_ok)
        rows.append(record)

    checks = {
        "rtol": rtol,
        "n_cases": len(rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "max_cubit_volume_rel_error": cubit_comparison["max_volume_rel_error"],
        "max_cubit_area_rel_error": cubit_comparison["max_area_rel_error"],
        "max_cubit_bbox_abs_error": cubit_comparison["max_bbox_abs_error"],
        "n_cubit_bbox_compared": cubit_comparison["n_bbox_compared"],
        "max_analytic_volume_rel_error": max(
            (row["analytic_volume_rel_error"] or 0.0) for row in rows
        ),
        "max_analytic_area_rel_error": max(
            (row["analytic_area_rel_error"] or 0.0) for row in rows
        ),
        "cubit_available": bool(cubit.get("available")),
        "cubit_bin_name": cubit.get("bin_name"),
        "bbox_atol": bbox_atol,
    }

    return {
        "kind": "build123d_cubit_volume_area_bbox_cross_validation",
        "validation_class": True,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "checks": checks,
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

    summary = build_summary(
        args.out_dir,
        args.cubit_bin,
        args.require_cubit,
        args.rtol,
        args.bbox_atol,
    )
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[build123d -> Cubit volume/area/bbox cross validation]")
    print(
        f"  cases={checks['n_cases']} passed={checks['n_passed']}, "
        f"Cubit available={checks['cubit_available']} ({checks['cubit_bin_name']})"
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
        f"  max analytic rel errors: volume={checks['max_analytic_volume_rel_error']:.3e}, "
        f"area={checks['max_analytic_area_rel_error']:.3e}"
    )
    print(f"[OK] wrote {args.summary}")
    return 0 if checks["n_passed"] == checks["n_cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
