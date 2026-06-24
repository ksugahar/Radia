"""Validation-class racetrack coil / plate / air multi-region model.

This is a moderate example, not a pytest test.  It checks a useful CAE
geometry contract for induction-heating and magnetostatic teaching models:

* build123d authors a racetrack coil, conductive plate, and surrounding air;
* analytic region volumes match the CAD volumes;
* the multi-region STEP -> Netgen -> Gmsh path preserves one volume group per
  region so solver-side materials can be assigned by name.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_racetrack_plate_air_region.py --quick
    python examples/build123d_netgen_gmsh_flow/validation_racetrack_plate_air_region.py
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

from _pipeline import run_pipeline_multi, save_record  # noqa: E402
from radia_mcp.build123d.modeling import racetrack_coil  # noqa: E402


OUT = HERE / "runs" / "validation_racetrack_plate_air_region"
SUMMARY_JSON = HERE / "validation_racetrack_plate_air_region_summary.json"
REGION_NAMES = ("plate", "coil", "air")


CASES = [
    {
        "label": "racetrack_plate_baseline",
        "plate_x": 80.0,
        "plate_y": 50.0,
        "plate_h": 4.0,
        "coil_length": 56.0,
        "coil_width": 32.0,
        "coil_band": 4.0,
        "coil_h": 3.0,
        "corner_radius": 8.0,
        "gap": 5.0,
        "air_x": 100.0,
        "air_y": 70.0,
        "air_h": 30.0,
        "air_center_z": 4.0,
        "maxh": 8.0,
    },
    {
        "label": "racetrack_plate_wide_coil",
        "plate_x": 96.0,
        "plate_y": 64.0,
        "plate_h": 5.0,
        "coil_length": 70.0,
        "coil_width": 42.0,
        "coil_band": 5.0,
        "coil_h": 3.5,
        "corner_radius": 10.0,
        "gap": 6.0,
        "air_x": 120.0,
        "air_y": 86.0,
        "air_h": 34.0,
        "air_center_z": 5.0,
        "maxh": 9.0,
    },
    {
        "label": "racetrack_plate_tall_gap",
        "plate_x": 72.0,
        "plate_y": 44.0,
        "plate_h": 3.0,
        "coil_length": 50.0,
        "coil_width": 28.0,
        "coil_band": 3.5,
        "coil_h": 2.8,
        "corner_radius": 7.0,
        "gap": 9.0,
        "air_x": 94.0,
        "air_y": 64.0,
        "air_h": 36.0,
        "air_center_z": 6.0,
        "maxh": 8.5,
    },
]


def _stage(record: dict, name: str) -> dict | None:
    for stage in record.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def _rounded_rectangle_area(length: float, width: float, radius: float) -> float:
    return length * width - (4.0 - math.pi) * radius * radius


def analytic_coil_volume(case: dict) -> float:
    r_outer = case["corner_radius"]
    band = case["coil_band"]
    r_inner = max(r_outer - band, 0.1)
    outer = _rounded_rectangle_area(case["coil_length"], case["coil_width"], r_outer)
    inner = _rounded_rectangle_area(
        case["coil_length"] - 2.0 * band,
        case["coil_width"] - 2.0 * band,
        r_inner,
    )
    return (outer - inner) * case["coil_h"]


def analytic_region_volumes(case: dict) -> dict[str, float]:
    plate = case["plate_x"] * case["plate_y"] * case["plate_h"]
    coil = analytic_coil_volume(case)
    air_box = case["air_x"] * case["air_y"] * case["air_h"]
    return {
        "plate": plate,
        "coil": coil,
        "air": air_box - plate - coil,
    }


def build_regions(case: dict):
    plate = Box(case["plate_x"], case["plate_y"], case["plate_h"]).solid()
    plate.label = "plate"
    coil_z = case["plate_h"] / 2.0 + case["gap"] + case["coil_h"] / 2.0
    coil = Pos(0, 0, coil_z) * racetrack_coil(
        case["coil_length"],
        case["coil_width"],
        case["coil_band"],
        case["coil_h"],
        case["corner_radius"],
        label="coil",
    )
    coil.label = "coil"
    air_box = Pos(0, 0, case["air_center_z"]) * Box(case["air_x"], case["air_y"], case["air_h"])
    air = (air_box - plate - coil).solid()
    air.label = "air"
    return [(plate, "plate"), (coil, "coil"), (air, "air")]


def _cad_volume_map(cad_stage: dict | None) -> dict[str, float]:
    if cad_stage is None:
        return {}
    return {row["name"]: row["volume"] for row in cad_stage.get("regions", [])}


def run_case(case: dict, out_dir: Path, volume_rtol: float) -> dict:
    record = run_pipeline_multi(build_regions(case), out_dir=out_dir, label=case["label"], maxh=case["maxh"])
    cad = _stage(record, "cad")
    mesh = _stage(record, "mesh")
    post = _stage(record, "post")
    expected = analytic_region_volumes(case)
    cad_volumes = _cad_volume_map(cad)
    rel_errors = {
        name: abs(cad_volumes[name] - expected[name]) / expected[name]
        for name in REGION_NAMES
    } if set(cad_volumes) == set(REGION_NAMES) else {}
    max_rel_error = max(rel_errors.values()) if rel_errors else None
    post_regions = post.get("regions", []) if post else []
    post_names = [row["name"] for row in post_regions]
    elem_counts = {row["name"]: row["n_elem"] for row in post_regions}
    air_positive = expected["air"] > 0.0
    region_order_ok = post_names == list(REGION_NAMES)

    record["case"] = case
    record["validation"] = {
        "kind": "racetrack_plate_air_region",
        "expected_volumes": expected,
        "cad_volumes": cad_volumes,
        "volume_rel_errors": rel_errors,
        "max_volume_rel_error": max_rel_error,
        "volume_rtol": volume_rtol,
        "post_region_names": post_names,
        "post_region_element_counts": elem_counts,
        "mesh_elements": mesh.get("ne") if mesh else None,
        "mesh_vertices": mesh.get("nv") if mesh else None,
        "air_positive": air_positive,
        "region_order_ok": region_order_ok,
        "passed": (
            record.get("status") == "ok"
            and air_positive
            and max_rel_error is not None
            and max_rel_error <= volume_rtol
            and region_order_ok
            and all(elem_counts.get(name, 0) > 0 for name in REGION_NAMES)
        ),
    }
    save_record(record, out_dir)
    return record


def _summary_row(record: dict) -> dict:
    validation = record.get("validation", {})
    return {
        "label": record.get("label"),
        "status": record.get("status"),
        "passed": validation.get("passed", False),
        "max_volume_rel_error": validation.get("max_volume_rel_error"),
        "post_region_names": validation.get("post_region_names"),
        "post_region_element_counts": validation.get("post_region_element_counts"),
        "mesh_elements": validation.get("mesh_elements"),
        "mesh_vertices": validation.get("mesh_vertices"),
        "expected_volumes": validation.get("expected_volumes"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="run only the baseline case")
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--volume-rtol", type=float, default=1.0e-7)
    args = parser.parse_args()

    cases = CASES[:1] if args.quick else CASES
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[validation_racetrack] {len(cases)} case(s), writing to {args.out_dir.resolve()}")

    t0 = time.perf_counter()
    records = []
    for i, case in enumerate(cases, 1):
        print(f"[{i}/{len(cases)}] {case['label']} ... ", end="", flush=True)
        record = run_case(case, args.out_dir, args.volume_rtol)
        records.append(record)
        validation = record.get("validation", {})
        err = validation.get("max_volume_rel_error")
        if validation.get("passed"):
            print(
                "OK  "
                f"max_rel_vol_err={err:.3e}  "
                f"regions={validation['post_region_names']}  "
                f"ne={validation['mesh_elements']}"
            )
        else:
            tail = record.get("error", "").splitlines()[-1:] or [""]
            print(f"FAIL  max_rel_vol_err={err} {tail[0]}".rstrip())

    summary = {
        "kind": "build123d_racetrack_plate_air_region_validation",
        "validation_class": True,
        "quick": args.quick,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "region_names": list(REGION_NAMES),
        "cases": [_summary_row(record) for record in records],
        "n_ok": sum(1 for record in records if record.get("validation", {}).get("passed")),
        "n_total": len(records),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        f"[validation_racetrack] {summary['n_ok']} ok / {summary['n_total']} total "
        f"in {summary['elapsed_seconds']} s"
    )
    print(f"[validation_racetrack] summary = {args.summary}")
    return 0 if summary["n_ok"] == summary["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
