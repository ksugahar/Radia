"""Validation-class laminated stack for build123d -> Netgen -> Gmsh.

This moderate example checks a CAE geometry contract that appears in many
magnetic and electric teaching models:

* build123d authors touching box layers as one ordered multi-region stack;
* CAD volumes match simple analytic layer volumes;
* the STEP -> Netgen -> Gmsh bridge preserves one named volume group per
  layer, so solver-side material assignment can stay readable.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_laminated_stack_region_sweep.py --quick
    python examples/build123d_netgen_gmsh_flow/validation_laminated_stack_region_sweep.py
"""

from __future__ import annotations

import argparse
import json
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


OUT = HERE / "runs" / "validation_laminated_stack_region_sweep"
SUMMARY_JSON = HERE / "validation_laminated_stack_region_sweep_summary.json"


CASES = [
    {
        "label": "lam_stack_4layer_balanced",
        "length": 30.0,
        "width": 12.0,
        "n_steel": 4,
        "steel_thickness": 1.0,
        "insulation_thickness": 0.20,
        "maxh": 1.30,
    },
    {
        "label": "lam_stack_5layer_high_fill",
        "length": 28.0,
        "width": 14.0,
        "n_steel": 5,
        "steel_thickness": 0.90,
        "insulation_thickness": 0.08,
        "maxh": 1.10,
    },
    {
        "label": "lam_stack_3layer_thick_insulation",
        "length": 24.0,
        "width": 10.0,
        "n_steel": 3,
        "steel_thickness": 0.80,
        "insulation_thickness": 0.35,
        "maxh": 1.00,
    },
]


def _stage(record: dict, name: str) -> dict | None:
    for stage in record.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def _region_names(case: dict) -> list[str]:
    names: list[str] = []
    for i in range(case["n_steel"]):
        names.append(f"steel_{i + 1:02d}")
        if i < case["n_steel"] - 1:
            names.append(f"insulation_{i + 1:02d}")
    return names


def _total_height(case: dict) -> float:
    return (
        case["n_steel"] * case["steel_thickness"]
        + (case["n_steel"] - 1) * case["insulation_thickness"]
    )


def _solid_box(length: float, width: float, height: float, z_center: float, label: str):
    part = (Pos(0.0, 0.0, z_center) * Box(length, width, height)).solid()
    part.label = label
    return part


def build_regions(case: dict):
    regions = []
    z0 = -0.5 * _total_height(case)
    z = z0
    for i in range(case["n_steel"]):
        steel_name = f"steel_{i + 1:02d}"
        steel_t = case["steel_thickness"]
        regions.append((
            _solid_box(case["length"], case["width"], steel_t, z + 0.5 * steel_t, steel_name),
            steel_name,
        ))
        z += steel_t

        if i < case["n_steel"] - 1:
            ins_name = f"insulation_{i + 1:02d}"
            ins_t = case["insulation_thickness"]
            regions.append((
                _solid_box(case["length"], case["width"], ins_t, z + 0.5 * ins_t, ins_name),
                ins_name,
            ))
            z += ins_t
    return regions


def analytic_region_volumes(case: dict) -> dict[str, float]:
    area = case["length"] * case["width"]
    expected: dict[str, float] = {}
    for i in range(case["n_steel"]):
        expected[f"steel_{i + 1:02d}"] = area * case["steel_thickness"]
        if i < case["n_steel"] - 1:
            expected[f"insulation_{i + 1:02d}"] = area * case["insulation_thickness"]
    return expected


def lamination_fill_factor(case: dict) -> float:
    steel_h = case["n_steel"] * case["steel_thickness"]
    return steel_h / _total_height(case)


def run_case(case: dict, out_dir: Path, volume_rtol: float) -> dict:
    label = case["label"]
    region_names = _region_names(case)
    record = run_pipeline_multi(build_regions(case), out_dir=out_dir, label=label, maxh=case["maxh"])
    cad = _stage(record, "cad")
    mesh = _stage(record, "mesh")
    post = _stage(record, "post")

    expected = analytic_region_volumes(case)
    cad_regions = cad.get("regions", []) if cad else []
    got = {row["name"]: row["volume"] for row in cad_regions}
    rel_errors = (
        {name: abs(got[name] - expected[name]) / expected[name] for name in region_names}
        if set(got) == set(region_names)
        else {}
    )
    max_rel_error = max(rel_errors.values()) if rel_errors else None
    post_regions = post.get("regions", []) if post else []
    post_names = [row["name"] for row in post_regions]
    elem_counts = {row["name"]: row["n_elem"] for row in post_regions}

    stack_height = _total_height(case)
    expected_total_volume = case["length"] * case["width"] * stack_height
    cad_total_volume = sum(got.values()) if got else None
    total_volume_rel_error = (
        abs(cad_total_volume - expected_total_volume) / expected_total_volume
        if cad_total_volume is not None
        else None
    )

    record["case"] = case
    record["validation"] = {
        "kind": "laminated_touching_region_stack",
        "region_names": region_names,
        "n_regions": len(region_names),
        "stack_height": stack_height,
        "lamination_fill_factor": lamination_fill_factor(case),
        "expected_total_volume": expected_total_volume,
        "cad_total_volume": cad_total_volume,
        "total_volume_rel_error": total_volume_rel_error,
        "expected_volumes": expected,
        "cad_volumes": got,
        "volume_rel_errors": rel_errors,
        "max_volume_rel_error": max_rel_error,
        "volume_rtol": volume_rtol,
        "post_region_names": post_names,
        "post_region_element_counts": elem_counts,
        "mesh_elements": mesh.get("ne") if mesh else None,
        "mesh_vertices": mesh.get("nv") if mesh else None,
        "passed": (
            record.get("status") == "ok"
            and max_rel_error is not None
            and max_rel_error <= volume_rtol
            and total_volume_rel_error is not None
            and total_volume_rel_error <= volume_rtol
            and post_names == region_names
            and all(elem_counts.get(name, 0) > 0 for name in region_names)
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
        "n_regions": validation.get("n_regions"),
        "lamination_fill_factor": validation.get("lamination_fill_factor"),
        "stack_height": validation.get("stack_height"),
        "max_volume_rel_error": validation.get("max_volume_rel_error"),
        "total_volume_rel_error": validation.get("total_volume_rel_error"),
        "post_region_names": validation.get("post_region_names"),
        "post_region_element_counts": validation.get("post_region_element_counts"),
        "mesh_elements": validation.get("mesh_elements"),
        "mesh_vertices": validation.get("mesh_vertices"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="run only the baseline case")
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--volume-rtol", type=float, default=1.0e-8)
    args = parser.parse_args()

    cases = CASES[:1] if args.quick else CASES
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[validation_laminated] {len(cases)} case(s), writing to {args.out_dir.resolve()}")

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
                f"regions={validation['n_regions']}  "
                f"fill={validation['lamination_fill_factor']:.6f}  "
                f"ne={validation['mesh_elements']}"
            )
        else:
            tail = ""
            if record.get("error"):
                tail = record["error"].splitlines()[-1]
            print(f"FAIL  max_rel_vol_err={err} {tail}".rstrip())

    summary = {
        "kind": "build123d_laminated_stack_region_validation",
        "validation_class": True,
        "quick": args.quick,
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "cases": [_summary_row(record) for record in records],
        "n_ok": sum(1 for record in records if record.get("validation", {}).get("passed")),
        "n_total": len(records),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        f"[validation_laminated] {summary['n_ok']} ok / {summary['n_total']} total "
        f"in {summary['elapsed_seconds']} s"
    )
    print(f"[validation_laminated] summary = {args.summary}")
    return 0 if summary["n_ok"] == summary["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
