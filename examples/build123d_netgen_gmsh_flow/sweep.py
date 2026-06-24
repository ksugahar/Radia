"""
sweep.py — parametric sweep driving the build123d -> Netgen -> Gmsh
pipeline many times, to exercise the lab mcp-servers on varied inputs.

Covers three geometry families (IH workpiece+coil / parametric E-core /
dipole yoke quarter) with a small grid of parameters each. Every run
writes .brep + .msh + _post.msh + .json into ./runs/ under a unique
label, so later inspection + mcp-server driven analysis can be replayed
without rebuilding the CAD.

Run:
    python sweep.py                # full sweep (~20 cases)
    python sweep.py --quick        # one case per family (smoke)
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from build123d import (
    Align, Axis, Box, BuildSketch, Cylinder, Plane, Pos, Rectangle,
    RegularPolygon, extrude, loft, revolve, CenterArc, Line,
    BuildLine, make_face, mirror,
)

from _pipeline import run_pipeline, run_pipeline_multi, save_record


OUT = Path("./runs")
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Geometry family 1: IH workpiece + coil (co-axial cylinders)
# ---------------------------------------------------------------------------

def ih_workpiece_plus_coil(r_wp=25.0, h_wp=40.0,
                           r_coil_in=32.0, r_coil_out=38.0,
                           h_coil=50.0):
    """Cylindrical steel workpiece + surrounding single-turn ring coil.
    All mm. Returns a Compound-like Part (volumes boolean-unioned)."""
    workpiece = Cylinder(radius=r_wp, height=h_wp)
    ring = Cylinder(radius=r_coil_out, height=h_coil) \
        - Cylinder(radius=r_coil_in, height=h_coil)
    # Stack axially so they don't overlap in z.
    part = workpiece + Pos(0, 0, h_wp + 2.0) * ring
    part.label = "ih_wp_coil"
    return part


# ---------------------------------------------------------------------------
# Geometry family 2: parametric E-core transformer leg
# ---------------------------------------------------------------------------

def e_core(W=60.0, H=40.0, D=20.0, t=8.0):
    """E-shape core extruded by depth D. All mm."""
    with BuildSketch() as sk:
        Rectangle(W, t, align=(Align.CENTER, Align.MIN))
        Pos(-W/2 + t/2, t) * Rectangle(t, H - t, align=(Align.CENTER, Align.MIN))
        Pos(0, t) * Rectangle(t, H - t, align=(Align.CENTER, Align.MIN))
        Pos(W/2 - t/2, t) * Rectangle(t, H - t, align=(Align.CENTER, Align.MIN))
    core = extrude(sk.sketch, amount=D)
    core.label = "e_core"
    return core


# ---------------------------------------------------------------------------
# Geometry family 3: dipole yoke quarter (accelerator magnet)
# ---------------------------------------------------------------------------

def dipole_quarter(R_bore=30.0, R_outer=150.0, L_half=200.0):
    """Quarter-annulus yoke extruded to L_half. All mm."""
    with BuildSketch(Plane.XZ) as yoke_sk:
        with BuildLine():
            Line((R_bore, 0), (R_outer, 0))
            CenterArc((0, 0), R_outer, 0, 90)
            Line((0, R_outer), (0, R_bore))
            CenterArc((0, 0), R_bore, 90, -90)
        make_face()
    yoke = extrude(yoke_sk.sketch, amount=L_half)
    yoke.label = "dipole_quarter_yoke"
    return yoke


# ---------------------------------------------------------------------------
# Sweep specs
# ---------------------------------------------------------------------------

def ih_grid():
    """Vary workpiece size + coil inner radius."""
    for r_wp in (20.0, 25.0, 30.0):
        for r_in in (r_wp + 5.0, r_wp + 10.0):
            yield {
                "builder": ih_workpiece_plus_coil,
                "kwargs": {"r_wp": r_wp, "r_coil_in": r_in,
                           "r_coil_out": r_in + 5.0},
                "label": f"ih_rwp{int(r_wp)}_rin{int(r_in)}",
                "maxh": max(3.0, r_wp / 5.0),
            }


def ecore_grid():
    """Vary leg thickness and overall dimensions."""
    for W, H in ((50.0, 35.0), (70.0, 50.0)):
        for t in (6.0, 10.0):
            yield {
                "builder": e_core,
                "kwargs": {"W": W, "H": H, "t": t},
                "label": f"ecore_W{int(W)}_H{int(H)}_t{int(t)}",
                "maxh": t * 0.8,
            }


def dipole_grid():
    """Vary bore radius."""
    for R_bore in (20.0, 30.0, 40.0):
        yield {
            "builder": dipole_quarter,
            "kwargs": {"R_bore": R_bore},
            "label": f"dipole_bore{int(R_bore)}",
            "maxh": 10.0,
        }


# ---------------------------------------------------------------------------
# Multi-region family: IH workpiece + coil + air (3 labeled regions)
# ---------------------------------------------------------------------------

def ih_multi_regions(r_wp=15.0, h_wp=20.0,
                     r_coil_in=22.0, r_coil_out=28.0, h_coil=10.0,
                     r_air=60.0, h_air=80.0):
    """Return a list of (Part, name) tuples in strict order: workpiece,
    coil, air. Order is the physical-group-id contract for the pipeline."""
    workpiece = Cylinder(radius=r_wp, height=h_wp)
    workpiece.label = "workpiece"
    coil = Pos(0, 0, h_wp + 5.0) * (
        Cylinder(radius=r_coil_out, height=h_coil)
        - Cylinder(radius=r_coil_in, height=h_coil)
    )
    coil.label = "coil"
    air = Cylinder(radius=r_air, height=h_air) - workpiece - coil
    air.label = "air"
    return [(workpiece, "workpiece"),
            (coil, "coil"),
            (air, "air")]


def ih_multi_grid():
    """Vary workpiece radius and coil inner radius."""
    for r_wp in (12.0, 15.0, 20.0):
        for r_in in (r_wp + 5.0, r_wp + 10.0):
            yield {
                "mode": "multi",
                "builder": ih_multi_regions,
                "kwargs": {"r_wp": r_wp, "r_coil_in": r_in,
                           "r_coil_out": r_in + 6.0},
                "label": f"ih_multi_rwp{int(r_wp)}_rin{int(r_in)}",
                "maxh": max(4.0, r_wp / 4.0),
            }


def all_cases(quick: bool = False):
    """Yield all sweep cases. In --quick mode, take one per family."""
    grids = [ih_grid(), ecore_grid(), dipole_grid(), ih_multi_grid()]
    for g in grids:
        for case in g:
            yield case
            if quick:
                break


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="one case per family (smoke mode)")
    args = ap.parse_args()

    cases = list(all_cases(quick=args.quick))
    print(f"[sweep] {len(cases)} case(s), writing to {OUT.resolve()}")

    n_ok = 0
    n_err = 0
    summary = []
    for i, c in enumerate(cases, 1):
        label = c["label"]
        print(f"[{i}/{len(cases)}] {label} ... ", end="", flush=True)
        try:
            built = c["builder"](**c["kwargs"])
        except Exception as e:
            print(f"CAD build failed: {e}")
            n_err += 1
            summary.append({"label": label, "status": "cad_build_error",
                            "error": str(e)})
            continue

        is_multi = c.get("mode") == "multi"
        if is_multi:
            rec = run_pipeline_multi(built, OUT, label, maxh=c["maxh"])
        else:
            rec = run_pipeline(built, OUT, label, maxh=c["maxh"])
        save_record(rec, OUT)

        if rec["status"] == "ok":
            mesh = next(s for s in rec["stages"] if s["stage"] == "mesh")
            post = next(s for s in rec["stages"] if s["stage"] == "post")
            if is_multi:
                per_region = "  ".join(
                    f"{r['name']}={r['n_elem']}" for r in post["regions"]
                )
                print(f"OK  nv={mesh['nv']:>5}  ne={mesh['ne']:>5}  "
                      f"regions[{per_region}]")
            else:
                print(f"OK  nv={mesh['nv']:>5}  ne={mesh['ne']:>5}  "
                      f"post_nodes={post['n_nodes']}")
            n_ok += 1
        else:
            print(f"FAIL  ({rec.get('error','').splitlines()[-1] if rec.get('error') else '?'})")
            n_err += 1
        summary.append({"label": label, "status": rec["status"],
                        "mode": "multi" if is_multi else "single",
                        "kwargs": c["kwargs"], "maxh": c["maxh"]})

    print()
    print(f"[sweep] {n_ok} ok / {n_err} err / {len(cases)} total")
    (OUT / "sweep_summary.json").write_text(
        __import__("json").dumps(summary, indent=2))
    return 0 if n_err == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
