"""Validation-class segmented Halbach region sweep.

This is a heavier example, not a pytest test. It checks the geometry-to-mesh
contract for labelled segmented permanent-magnet rings:

* CAD volume matches the analytic annular-sector fill, including air gaps;
* each segment label carries the Mallinson easy-axis angle;
* the multi-region STEP -> Netgen -> Gmsh path preserves one volume group per
  segment so solver-side magnetization can be assigned by label.

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_halbach_region_sweep.py --quick
    python validation_test/build123d_netgen_gmsh_flow/validation_halbach_region_sweep.py
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

from radia_mcp.build123d.pipeline import run_pipeline_multi, save_record  # noqa: E402
from radia_mcp.build123d.archetypes import (  # noqa: E402
    halbach_ring,
    parse_magnetization,
)


OUT = HERE / "runs" / "validation_halbach_region_sweep"
SUMMARY_JSON = HERE / "validation_halbach_region_sweep_summary.json"


CASES = [
    {
        "label": "halbach_dipole_8_gapless",
        "r_in": 12.0,
        "r_out": 18.0,
        "height": 6.0,
        "n_segments": 8,
        "pole_pairs": 1,
        "gap_deg": 0.0,
        "maxh": 5.0,
    },
    {
        "label": "halbach_quad_12_gap1",
        "r_in": 10.0,
        "r_out": 16.0,
        "height": 5.0,
        "n_segments": 12,
        "pole_pairs": 2,
        "gap_deg": 1.0,
        "maxh": 4.5,
    },
    {
        "label": "halbach_dipole_16_gap0p5",
        "r_in": 14.0,
        "r_out": 22.0,
        "height": 5.5,
        "n_segments": 16,
        "pole_pairs": 1,
        "gap_deg": 0.5,
        "maxh": 5.0,
    },
]


def _stage(record: dict, name: str) -> dict | None:
    for stage in record.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def _angle_error_deg(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def _expected_volume(case):
    span = 360.0 / case["n_segments"] - case["gap_deg"]
    fill = case["n_segments"] * span / 360.0
    full_annulus = math.pi * (case["r_out"] ** 2 - case["r_in"] ** 2) * case["height"]
    return full_annulus * fill


def _expected_angles(case):
    span = 360.0 / case["n_segments"] - case["gap_deg"]
    rows = []
    for k in range(case["n_segments"]):
        theta_c = k * 360.0 / case["n_segments"] + 0.5 * span
        rows.append(((case["pole_pairs"] + 1) * theta_c) % 360.0)
    return rows


def run_case(case: dict, out_dir: Path, volume_rtol: float) -> dict:
    label = case["label"]
    hb = halbach_ring(
        case["r_in"],
        case["r_out"],
        case["height"],
        case["n_segments"],
        pole_pairs=case["pole_pairs"],
        gap_deg=case["gap_deg"],
        name=label,
    )
    children = list(hb.children)
    regions = [(child, child.label) for child in children]

    record = run_pipeline_multi(regions, out_dir=out_dir, label=label, maxh=case["maxh"])
    cad = _stage(record, "cad")
    mesh = _stage(record, "mesh")
    post = _stage(record, "post")

    expected_angles = _expected_angles(case)
    got_angles = [parse_magnetization(child.label) for child in children]
    angle_errors = [
        _angle_error_deg(got, expected)
        for got, expected in zip(got_angles, expected_angles)
    ]

    cad_total = sum(row["volume"] for row in cad.get("regions", [])) if cad else None
    expected_volume = _expected_volume(case)
    volume_rel_error = (
        abs(cad_total - expected_volume) / expected_volume
        if cad_total is not None and expected_volume > 0.0
        else None
    )
    magnetization_labels = [
        child.label for child, angle in zip(children, got_angles) if angle is not None
    ]
    region_count = post.get("n_regions") if post else None

    record["case"] = case
    record["validation"] = {
        "kind": "segmented_halbach_region_mesh",
        "expected_volume": expected_volume,
        "cad_total_volume": cad_total,
        "volume_rel_error": volume_rel_error,
        "volume_rtol": volume_rtol,
        "fill_factor": expected_volume / (
            math.pi * (case["r_out"] ** 2 - case["r_in"] ** 2) * case["height"]
        ),
        "max_mallinson_angle_error_deg": max(angle_errors),
        "magnetization_labels": magnetization_labels,
        "region_count": region_count,
        "mesh_elements": mesh.get("ne") if mesh else None,
        "mesh_vertices": mesh.get("nv") if mesh else None,
        "passed": (
            record.get("status") == "ok"
            and volume_rel_error is not None
            and volume_rel_error <= volume_rtol
            and max(angle_errors) <= 1e-9
            and len(magnetization_labels) == case["n_segments"]
            and region_count == case["n_segments"]
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
        "volume_rel_error": validation.get("volume_rel_error"),
        "fill_factor": validation.get("fill_factor"),
        "max_mallinson_angle_error_deg": validation.get("max_mallinson_angle_error_deg"),
        "region_count": validation.get("region_count"),
        "mesh_elements": validation.get("mesh_elements"),
        "mesh_vertices": validation.get("mesh_vertices"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="run only the 8-segment baseline")
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--volume-rtol", type=float, default=1.0e-6)
    args = parser.parse_args()

    cases = CASES[:1] if args.quick else CASES
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[validation_halbach] {len(cases)} case(s), writing to {args.out_dir.resolve()}")

    t0 = time.perf_counter()
    records = []
    for i, case in enumerate(cases, 1):
        print(f"[{i}/{len(cases)}] {case['label']} ... ", end="", flush=True)
        record = run_case(case, args.out_dir, args.volume_rtol)
        records.append(record)
        validation = record.get("validation", {})
        err = validation.get("volume_rel_error")
        if validation.get("passed"):
            print(
                "OK  "
                f"rel_vol_err={err:.3e}  "
                f"regions={validation['region_count']}  "
                f"ne={validation['mesh_elements']}"
            )
        else:
            tail = ""
            if record.get("error"):
                tail = record["error"].splitlines()[-1]
            print(f"FAIL  rel_vol_err={err} {tail}".rstrip())

    summary = {
        "kind": "build123d_halbach_region_validation",
        "validation_class": True,
        "quick": args.quick,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "cases": [_summary_row(record) for record in records],
        "n_ok": sum(1 for record in records if record.get("validation", {}).get("passed")),
        "n_total": len(records),
    }
    summary_path = args.summary
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        f"[validation_halbach] {summary['n_ok']} ok / {summary['n_total']} total "
        f"in {summary['elapsed_seconds']} s"
    )
    print(f"[validation_halbach] summary = {summary_path}")
    return 0 if summary["n_ok"] == summary["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
