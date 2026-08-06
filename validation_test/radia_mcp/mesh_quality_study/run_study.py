"""Mesh-quality study: netgen vs cubit, equal-h AND equal-budget.

Promoted from C:/temp/mesh_quality_study (2026-08-06) with its
committed results JSON (Data Persistence Policy). Re-run with
`python run_study.py` (requires Cubit + netgen + gmsh + build123d;
scratch meshes land in artifacts/, gitignored).
Per-route meshes are produced through the SAME machinery as
cubit_netgen_quality_compare and judged by the same gmsh minSICN
referee; this script adds:
  * a geometry matrix of lab-relevant shapes,
  * histogram-based distribution metrics (%<0.3, %>=0.9),
  * an EQUAL-BUDGET comparison: netgen maxh calibrated so its element
    count matches cubit_tet's within ~10%, answering "same element
    budget, which mesher is better?" (equal-h flattered netgen's
    economy: ~2x fewer elements at the same h).

Quality-class run (correctness, not timing) -- LAB execution allowed.
"""
import json
import os
import platform
import time
from datetime import datetime

sys.path.insert(0, r"S:\Radia\01_GitHub\packages\radia-mcp\src")

from build123d import (Box, Compound, Cylinder, Pos, Sphere,
                       export_step)
from radia_mcp.build123d.archetypes import c_core, halbach_ring, slotted_stator
from radia_mcp.cubit.server import _netgen_mesh_to_msh, _run_batch
from radia_mcp.gmsh.msh_inspect import mesh_quality

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------------
# Geometry matrix (STEP authored via the lab helpers)
# ----------------------------------------------------------------------

def build_geometries():
    geoms = {}

    def emit(name, part):
        p = os.path.join(OUT, f"{name}.step")
        export_step(part, p)
        geoms[name] = p

    emit("sphere", Sphere(1.0))
    emit("thin_plate", Box(50, 50, 1))
    emit("c_core", c_core(width=80, height=60, depth=25, leg=15, gap=8))
    emit("halbach", halbach_ring(r_in=20, r_out=30, h=10, n_segments=12))
    emit("stator", slotted_stator(r_bore=15, r_yoke=30, n_slots=6,
                                  slot_depth=8, slot_span_deg=20, h=12))
    gap_a = Box(10, 10, 5)
    gap_b = Pos(0, 0, 5.5) * Box(10, 10, 5)     # 0.5 gap between bodies
    emit("gap_pair", Compound(children=[gap_a, gap_b]))
    emit("cyl_hole", Cylinder(5, 10) - Cylinder(0.8, 22))  # through-hole
    return geoms


# Per-geometry mesh sizes (chosen ~1/8 bbox, then refined once)
SIZES = {
    "sphere":     [0.5, 0.3],
    "thin_plate": [3.0, 1.5],
    "c_core":     [8.0, 4.0],
    "halbach":    [4.0, 2.0],
    "stator":     [4.0, 2.0],
    "gap_pair":   [2.0, 1.0],
    "cyl_hole":   [1.5, 0.8],
}


# ----------------------------------------------------------------------
# Route runners (reuse the compare tool's machinery)
# ----------------------------------------------------------------------

def _hist_metrics(q):
    h = q.get("histogram") or {}
    counts = h.get("counts") or []
    total = sum(counts) or 1
    # edges are [0, .1, .3, .6, .9, 1.0]
    below03 = sum(counts[:2]) / total * 100.0
    top = counts[-1] / total * 100.0 if counts else 0.0
    return round(below03, 2), round(top, 2)


def _summarize(q):
    bts = q.get("by_type") or []
    if not bts:
        return None
    n = sum(bt["n_elements"] for bt in bts)
    mn = min(bt["min_quality"] for bt in bts)
    mean = sum(bt["mean_quality"] * bt["n_elements"] for bt in bts) / n
    neg = sum(bt["negative"] for bt in bts)
    b03, top = _hist_metrics(q)
    return {"n": n, "min": round(mn, 4), "mean": round(mean, 4),
            "neg": neg, "pct_below_0.3": b03, "pct_ge_0.9": top,
            "elements": "+".join(sorted({bt["name"] for bt in bts}))}


def run_netgen(step, maxh, tag):
    msh = os.path.join(OUT, f"{tag}_netgen.msh")
    from pathlib import Path
    _netgen_mesh_to_msh(Path(step), maxh, Path(msh))
    return _summarize(mesh_quality(msh))


def run_cubit(step, size, tag, scheme):
    msh = os.path.join(OUT, f"{tag}_cubit_{scheme}.msh")
    cmds = []
    if scheme == "tet":
        cmds.append("volume all scheme tetmesh")
    cmds += [f"volume all size {size}", "mesh volume all",
             "block 1 add volume all", 'block 1 name "mesh"',
             f'export gmsh "{msh.replace(os.sep, "/")}" overwrite']
    r = _run_batch(step, cmds, timeout_s=600)
    if r.get("status") != "ok":
        return {"error": (r.get("error") or "")[:120]}
    return _summarize(mesh_quality(msh))


def calibrate_netgen_to_count(step, start_maxh, target_n, tag,
                              tol=0.10, iters=4):
    """Scale maxh so netgen's element count matches target within tol."""
    maxh = start_maxh
    best = None
    for _ in range(iters):
        s = run_netgen(step, maxh, f"{tag}_cal")
        best = (maxh, s)
        ratio = s["n"] / target_n
        if abs(ratio - 1.0) <= tol:
            break
        maxh *= ratio ** (1.0 / 3.0)     # n ~ 1/h^3
    return best


# ----------------------------------------------------------------------
# Study
# ----------------------------------------------------------------------

def main():
    t0 = time.time()
    geoms = build_geometries()
    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "referee": "gmsh minSICN (radia_mcp.gmsh.msh_inspect)",
               "note": ("quality-class study (no timing); "
                        "equal_budget = netgen maxh calibrated to "
                        "cubit_tet element count +/-10%"),
               "cases": []}

    for name, step in geoms.items():
        for size in SIZES[name]:
            tag = f"{name}_{size}"
            row = {"geometry": name, "size": size}
            row["cubit_tet"] = run_cubit(step, size, tag, "tet")
            row["netgen_equal_h"] = run_netgen(step, size, tag)
            if name in ("thin_plate", "gap_pair"):
                row["cubit_hex"] = run_cubit(step, size, tag, "hex")
            # equal-budget netgen (calibrated to cubit_tet count)
            ct = row["cubit_tet"]
            if ct and "n" in ct:
                maxh_eq, s_eq = calibrate_netgen_to_count(
                    step, size, ct["n"], tag)
                s_eq["calibrated_maxh"] = round(maxh_eq, 4)
                row["netgen_equal_budget"] = s_eq
            results["cases"].append(row)
            print(f"[{name} size={size}]")
            for k in ("cubit_tet", "cubit_hex", "netgen_equal_h",
                      "netgen_equal_budget"):
                if k in row and row[k]:
                    s = row[k]
                    if "error" in s:
                        print(f"  {k:20s} ERROR {s['error']}")
                    else:
                        extra = (f" (maxh={s['calibrated_maxh']})"
                                 if "calibrated_maxh" in s else "")
                        print(f"  {k:20s} n={s['n']:6d} min={s['min']:.3f} "
                              f"mean={s['mean']:.3f} <0.3:{s['pct_below_0.3']:5.2f}% "
                              f">=0.9:{s['pct_ge_0.9']:5.2f}%{extra}")

    out = os.path.join(_HERE, "results_mesh_quality_study.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"\nsaved {out}  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
