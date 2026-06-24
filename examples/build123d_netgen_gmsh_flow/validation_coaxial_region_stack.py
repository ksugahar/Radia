"""Validation-class coaxial multi-region stack for build123d -> Netgen -> Gmsh.

This is a moderate example, not a pytest test. It checks a simple but useful
CAE geometry contract:

* build123d authors disjoint touching cylindrical regions;
* CAD volumes match the analytic coaxial shell formulas;
* the multi-region STEP -> Netgen -> Gmsh path preserves one volume group per
  region so solver-side materials can be assigned by name.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_coaxial_region_stack.py --quick
    python examples/build123d_netgen_gmsh_flow/validation_coaxial_region_stack.py
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

from build123d import Cylinder  # noqa: E402

from _pipeline import run_pipeline_multi, save_record  # noqa: E402
from radia_mcp.build123d.modeling import tube  # noqa: E402


OUT = HERE / "runs" / "validation_coaxial_region_stack"
SUMMARY_JSON = HERE / "validation_coaxial_region_stack_summary.json"


CASES = [
    {
        "label": "coax_short_baseline",
        "r_inner": 2.0,
        "r_dielectric": 5.0,
        "r_outer": 6.5,
        "height": 12.0,
        "maxh": 2.2,
    },
    {
        "label": "coax_tall_slender",
        "r_inner": 1.2,
        "r_dielectric": 4.0,
        "r_outer": 5.5,
        "height": 24.0,
        "maxh": 2.0,
    },
    {
        "label": "coax_thick_shield",
        "r_inner": 2.5,
        "r_dielectric": 4.5,
        "r_outer": 8.0,
        "height": 14.0,
        "maxh": 2.4,
    },
]
REGION_NAMES = ("inner_conductor", "dielectric", "outer_conductor")


def _stage(record: dict, name: str) -> dict | None:
    for stage in record.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def _solid_cylinder(radius: float, height: float, label: str):
    part = Cylinder(radius=radius, height=height).solid()
    part.label = label
    return part


def build_regions(case: dict):
    inner = _solid_cylinder(case["r_inner"], case["height"], REGION_NAMES[0])
    dielectric = tube(case["r_inner"], case["r_dielectric"], case["height"], label=REGION_NAMES[1])
    outer = tube(case["r_dielectric"], case["r_outer"], case["height"], label=REGION_NAMES[2])
    return [(inner, REGION_NAMES[0]), (dielectric, REGION_NAMES[1]), (outer, REGION_NAMES[2])]


def analytic_region_volumes(case: dict) -> dict[str, float]:
    h = case["height"]
    ri = case["r_inner"]
    rd = case["r_dielectric"]
    ro = case["r_outer"]
    return {
        REGION_NAMES[0]: math.pi * ri * ri * h,
        REGION_NAMES[1]: math.pi * (rd * rd - ri * ri) * h,
        REGION_NAMES[2]: math.pi * (ro * ro - rd * rd) * h,
    }


def run_case(case: dict, out_dir: Path, volume_rtol: float) -> dict:
    label = case["label"]
    record = run_pipeline_multi(build_regions(case), out_dir=out_dir, label=label, maxh=case["maxh"])
    cad = _stage(record, "cad")
    mesh = _stage(record, "mesh")
    post = _stage(record, "post")

    expected = analytic_region_volumes(case)
    cad_regions = cad.get("regions", []) if cad else []
    got = {row["name"]: row["volume"] for row in cad_regions}
    rel_errors = {
        name: abs(got[name] - expected[name]) / expected[name]
        for name in REGION_NAMES
    } if set(got) == set(REGION_NAMES) else {}
    max_rel_error = max(rel_errors.values()) if rel_errors else None
    post_regions = post.get("regions", []) if post else []
    post_names = [row["name"] for row in post_regions]
    elem_counts = {row["name"]: row["n_elem"] for row in post_regions}

    record["case"] = case
    record["validation"] = {
        "kind": "coaxial_touching_region_stack",
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
            and post_names == list(REGION_NAMES)
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
    print(f"[validation_coax] {len(cases)} case(s), writing to {args.out_dir.resolve()}")

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
            tail = ""
            if record.get("error"):
                tail = record["error"].splitlines()[-1]
            print(f"FAIL  max_rel_vol_err={err} {tail}".rstrip())

    summary = {
        "kind": "build123d_coaxial_region_stack_validation",
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
        f"[validation_coax] {summary['n_ok']} ok / {summary['n_total']} total "
        f"in {summary['elapsed_seconds']} s"
    )
    print(f"[validation_coax] summary = {args.summary}")
    return 0 if summary["n_ok"] == summary["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
