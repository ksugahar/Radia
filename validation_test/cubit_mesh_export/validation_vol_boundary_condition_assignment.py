"""Validation-class boundary-condition assignment audit for tri/tet .vol meshes.

Run:

    python validation_test/cubit_mesh_export/validation_vol_boundary_condition_assignment.py

The audit maps plain boundary-condition labels onto Netgen boundary numbers or
names and reports missing/unknown assignments before a solver script turns those
labels into operators.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_boundary_condition_assignment_summary.json"

BOX_SIX_BOUNDARY_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
6
1 1 0 1 1
2 1 0 1 1
3 1 0 1 1
4 1 0 1 1
5 1 0 1 1
6 1 0 1 1
surfaceelements
12
5 5 1 0 3 1 3 2
5 5 1 0 3 1 4 3
6 6 1 0 3 5 6 7
6 6 1 0 3 5 7 8
3 3 1 0 3 1 2 6
3 3 1 0 3 1 6 5
4 4 1 0 3 4 7 3
4 4 1 0 3 4 8 7
1 1 1 0 3 1 5 8
1 1 1 0 3 1 8 4
2 2 1 0 3 2 3 7
2 2 1 0 3 2 7 6
volumeelements
12
1 4 1 3 2 9
1 4 1 4 3 9
1 4 5 6 7 9
1 4 5 7 8 9
1 4 1 2 6 9
1 4 1 6 5 9
1 4 4 7 3 9
1 4 4 8 7 9
1 4 1 5 8 9
1 4 1 8 4 9
1 4 2 3 7 9
1 4 2 7 6 9
points
9
0 0 0
2 0 0
2 3 0
0 3 0
0 0 5
2 0 5
2 3 5
0 3 5
1 1.5 2.5
pointelements
0
materials
1
1 air
bcnames
6
1 xmin
2 xmax
3 ymin
4 ymax
5 zmin
6 zmax
endmesh
"""


def build_summary() -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(
        BOX_SIX_BOUNDARY_VOL,
        source="embedded_box_six_boundary.vol",
    )
    clean = mesh.boundary_condition_assignment_summary(
        {
            "zmax": "impedance",
            "zmin": "dirichlet",
            1: "symmetry",
        },
        default_condition="open",
    )
    typo = mesh.boundary_condition_assignment_summary(
        {
            "zmax": "impedance",
            "not_a_boundary": "neumann",
        },
        default_condition="open",
    )
    missing = mesh.boundary_condition_assignment_summary({"zmax": "impedance"})

    clean_by_name = {row["name"]: row for row in clean["rows"]}
    checks = {
        "clean_ok": clean["ok"],
        "clean_boundary_count": clean["boundary_count"],
        "clean_condition_counts": clean["condition_counts"],
        "zmax_condition": clean_by_name["zmax"]["condition"],
        "zmax_trace_node_count": clean_by_name["zmax"]["trace_node_count"],
        "xmin_condition_source": clean_by_name["xmin"]["condition_source"],
        "typo_unknown_condition_keys": typo["unknown_condition_keys"],
        "typo_ok": typo["ok"],
        "missing_boundary_count": missing["missing_boundary_count"],
    }

    assert checks["clean_ok"] is True
    assert checks["clean_boundary_count"] == 6
    assert checks["clean_condition_counts"] == {
        "dirichlet": 1,
        "impedance": 1,
        "open": 3,
        "symmetry": 1,
    }
    assert checks["zmax_condition"] == "impedance"
    assert checks["zmax_trace_node_count"] == 4
    assert checks["xmin_condition_source"] == "boundary_number"
    assert checks["typo_unknown_condition_keys"] == ["not_a_boundary"]
    assert checks["typo_ok"] is False
    assert checks["missing_boundary_count"] == 5

    return {
        "kind": "vol_boundary_condition_assignment_validation",
        "validation_class": True,
        "learning_theme": (
            "named .vol boundaries should be audited against condition labels "
            "before solver operators are built"
        ),
        "mesh_summary": mesh.summary(),
        "checks": checks,
        "clean_assignment": clean,
        "typo_assignment": typo,
        "missing_assignment": missing,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[vol boundary condition assignment]")
    print(
        f"  clean_ok={checks['clean_ok']} boundary_count={checks['clean_boundary_count']} "
        f"conditions={checks['clean_condition_counts']}"
    )
    print(f"  zmax={checks['zmax_condition']} trace_nodes={checks['zmax_trace_node_count']}")
    print(f"  typo_unknown_keys={checks['typo_unknown_condition_keys']}")
    print(f"  missing_boundary_count={checks['missing_boundary_count']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
