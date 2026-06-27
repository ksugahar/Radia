"""Validation-class helical conductor sweep for the build123d -> Netgen -> Gmsh flow.

This is intentionally an example/validation script, not a pytest test:
the cases are heavier than CI smoke tests and are meant to produce reusable
CAD/mesh/post records for learning and later solver cross-validation.

Run from this directory or the repository root:

    python validation_test/build123d_netgen_gmsh_flow/validation_helix_mesh_sweep.py --quick
    python validation_test/build123d_netgen_gmsh_flow/validation_helix_mesh_sweep.py
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

from build123d import Circle  # noqa: E402

from radia_mcp.build123d.pipeline import run_pipeline, save_record  # noqa: E402
from radia_mcp.build123d.modeling import coil, round_wire_helix_metrics  # noqa: E402


OUT = HERE / "runs" / "validation_helix_mesh_sweep"


CASES = [
    {
        "label": "helix_4turn_baseline",
        "radius": 2.0,
        "wire_radius": 0.16,
        "pitch": 1.2,
        "height": 4.8,
        "maxh": 0.38,
    },
    {
        "label": "helix_8turn_tight_pitch",
        "radius": 2.0,
        "wire_radius": 0.11,
        "pitch": 0.65,
        "height": 5.2,
        "maxh": 0.30,
    },
    {
        "label": "helix_thick_wire_large_radius",
        "radius": 2.7,
        "wire_radius": 0.24,
        "pitch": 1.1,
        "height": 4.4,
        "maxh": 0.45,
    },
    {
        "label": "helix_slender_long",
        "radius": 1.6,
        "wire_radius": 0.10,
        "pitch": 0.7,
        "height": 5.6,
        "maxh": 0.28,
    },
]


def _stage(record: dict, name: str) -> dict | None:
    for stage in record.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def run_case(case: dict, out_dir: Path, volume_rtol: float) -> dict:
    label = case["label"]
    part = coil(
        Circle(case["wire_radius"]),
        pitch=case["pitch"],
        height=case["height"],
        radius=case["radius"],
        label=label,
    )
    metrics = round_wire_helix_metrics(
        radius=case["radius"],
        wire_radius=case["wire_radius"],
        pitch=case["pitch"],
        height=case["height"],
    )

    record = run_pipeline(part, out_dir=out_dir, label=label, maxh=case["maxh"])
    cad = _stage(record, "cad")
    mesh = _stage(record, "mesh")
    post = _stage(record, "post")

    cad_volume = float(cad["volume"]) if cad else None
    analytic_volume = float(metrics["conductor_volume"])
    volume_rel_error = (
        abs(cad_volume - analytic_volume) / analytic_volume
        if cad_volume is not None and analytic_volume > 0.0
        else None
    )

    record["case"] = case
    record["validation"] = {
        "kind": "round_wire_helix_volume_and_mesh",
        "analytic": metrics,
        "cad_volume": cad_volume,
        "volume_rel_error": volume_rel_error,
        "volume_rtol": volume_rtol,
        "passed": (
            record.get("status") == "ok"
            and volume_rel_error is not None
            and volume_rel_error <= volume_rtol
        ),
        "mesh_elements": mesh.get("ne") if mesh else None,
        "mesh_vertices": mesh.get("nv") if mesh else None,
        "post_nodes": post.get("n_nodes") if post else None,
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
        "mesh_elements": validation.get("mesh_elements"),
        "mesh_vertices": validation.get("mesh_vertices"),
        "post_nodes": validation.get("post_nodes"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run only the baseline case for a fast validation pass",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUT,
        help="directory for .brep/.msh/_post.msh/.json outputs",
    )
    parser.add_argument(
        "--volume-rtol",
        type=float,
        default=5.0e-3,
        help="relative CAD-volume tolerance against the analytic helix volume",
    )
    args = parser.parse_args()

    cases = CASES[:1] if args.quick else CASES
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[validation_helix] {len(cases)} case(s), writing to {args.out_dir.resolve()}")
    t0 = time.perf_counter()
    records = []
    for i, case in enumerate(cases, 1):
        label = case["label"]
        print(f"[{i}/{len(cases)}] {label} ... ", end="", flush=True)
        record = run_case(case, args.out_dir, args.volume_rtol)
        records.append(record)
        validation = record.get("validation", {})
        err = validation.get("volume_rel_error")
        if validation.get("passed"):
            print(
                "OK  "
                f"rel_vol_err={err:.3e}  "
                f"nv={validation['mesh_vertices']}  "
                f"ne={validation['mesh_elements']}  "
                f"post_nodes={validation['post_nodes']}"
            )
        else:
            tail = ""
            if record.get("error"):
                tail = record["error"].splitlines()[-1]
            print(f"FAIL  rel_vol_err={err} {tail}".rstrip())

    summary = {
        "kind": "build123d_helix_mesh_validation",
        "validation_class": True,
        "quick": args.quick,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "cases": [_summary_row(record) for record in records],
        "n_ok": sum(1 for record in records if record.get("validation", {}).get("passed")),
        "n_total": len(records),
    }
    summary_path = args.out_dir / "validation_helix_mesh_sweep_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(
        f"[validation_helix] {summary['n_ok']} ok / {summary['n_total']} total "
        f"in {summary['elapsed_seconds']} s"
    )
    print(f"[validation_helix] summary = {summary_path}")
    return 0 if summary["n_ok"] == summary["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
