"""build123d-generated lofted coil: PEEC vs BEM-A cross-validation.

This script demonstrates the workflow to hand off to a downstream user:

  1. Build a lofted single-turn coil with build123d
     (multi-profile ``Solid.make_loft`` + small angular gap).
  2. Save as STEP.
  3. Compute self-inductance via Radia PEEC (``calc_peec_inductance``).
  4. Compute self-inductance via the new intree BEM-A
     (``radia.bem.efie_rwg.solve_inductance_source_sink_intree``).
  5. Compare.

The geometry mirrors the v4.24.0 PEEC golden fixture
``tests/panels/golden/rect_torus_lofted_united.step`` (R = 50 mm,
8 mm radial x 6 mm vertical rect cross-section, 355 deg arc, 13
loft profiles -> 12 ruled lofts -> united into one solid).
PEEC golden L = 138.16 nH (band [120, 160]).

Run::

    python loft_coil_peec_vs_bem.py
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from typing import List, Tuple

import numpy as np

from build123d import (
    Solid, Polyline, Plane, Vector, Wire, Face, Compound,
    Pos, Rot,
)


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


# ============================================================
# Step 1: build123d lofted coil
# ============================================================

def build_loft_coil_step(
        out_path: str,
        R_major: float = 0.050,
        rect_radial: float = 0.008,
        rect_vertical: float = 0.006,
        gap_deg: float = 5.0,
        n_lofts: int = 12) -> str:
    """Build a 355 deg lofted rect-cross-section single-turn coil.

    Uses build123d ``Solid.make_loft`` between (n_lofts + 1) cross-section
    wires placed along an arc spine, then ``fuse`` to a single solid.
    """
    sweep_deg = 360.0 - gap_deg
    n_profiles = n_lofts + 1
    angles_deg = np.linspace(0.0, -sweep_deg, n_profiles)

    half_w = 0.5 * rect_radial
    half_h = 0.5 * rect_vertical

    profiles_w: List[Wire] = []
    for theta_deg in angles_deg:
        theta = math.radians(theta_deg)
        center = Vector(R_major * math.cos(theta),
                         R_major * math.sin(theta),
                         0.0)
        radial = Vector(math.cos(theta), math.sin(theta), 0.0)
        vertical = Vector(0.0, 0.0, 1.0)
        # Build the 4 corners of the rectangle in 3D directly.
        corners_3d = []
        for sign_w, sign_h in [(-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1)]:
            corners_3d.append(
                center + (sign_w * half_w) * radial
                       + (sign_h * half_h) * vertical)
        poly = Polyline(*corners_3d)
        wire = Wire(poly.edges())
        profiles_w.append(wire)

    # Loft consecutive pairs (ruled=True for clean planar segments)
    segments: List[Solid] = []
    for i in range(n_lofts):
        s = Solid.make_loft([profiles_w[i], profiles_w[i + 1]], ruled=True)
        segments.append(s)

    # Fuse into one solid
    fused = segments[0]
    for s in segments[1:]:
        fused = fused.fuse(s)

    # Export STEP via the module-level helper (build123d 0.10).
    from build123d import export_step
    export_step(fused, out_path)
    return out_path


# ============================================================
# Step 2: PEEC inductance via existing CLI
# ============================================================

def compute_peec_L(step_path: str, n_peri: int = 16,
                    frequency_Hz: float = 50e3,
                    current_A: float = 1.0,
                    sigma: float = 5.8e7) -> dict:
    """Run ``src/radia/panels/calc_peec_inductance.py`` as a subprocess
    on the STEP and parse the JSON result."""
    calc = os.path.join(ROOT, "src", "radia", "panels",
                         "calc_peec_inductance.py")
    cmd = [sys.executable, calc,
           "--peec-step", step_path,
           "--peec-n-peri", str(n_peri),
           "--frequency", str(frequency_Hz),
           "--current", str(current_A),
           "--coil-sigma", str(sigma)]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(
            f"PEEC calc returned {proc.returncode}\n"
            f"STDERR:\n{proc.stderr[-2000:]}")
    json_lines = [ln for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("{")]
    return json.loads(json_lines[-1])


# ============================================================
# Step 3: BEM-A via NGSolve OCC + intree saddle solver
# ============================================================

def compute_bem_a_L(step_path: str, maxh: float = 0.012,
                    regular_quad_degree: int = 5,
                    singular_n_q: int = 6) -> dict:
    """Mesh the STEP via NGSolve OCC ``Glue(faces)``, label the two
    smallest-mass faces as source/sink, run the intree BEM-A
    saddle-point EFIE solver."""
    from netgen.occ import OCCGeometry, Glue
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh, BND
    from radia.bem.efie_rwg import solve_inductance_source_sink_intree

    geo = OCCGeometry(step_path)
    shape = geo.shape

    # Cap classifier: 2 smallest mass = caps; source = y >= 0, sink = y < 0
    fs_sorted = sorted(enumerate(shape.faces), key=lambda x: x[1].mass)
    ca, fa = fs_sorted[0]
    cb, fb = fs_sorted[1]
    if abs(fa.center[1]) < abs(fb.center[1]):
        fa.name, fb.name = "source", "sink"
    else:
        fa.name, fb.name = "sink", "source"
    for i, f in enumerate(shape.faces):
        if i not in (ca, cb):
            f.name = "body"

    ngmesh = OCCGeometry(Glue(list(shape.faces))).GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=1.0,
                              segmentsperedge=1))
    mesh = Mesh(ngmesh)

    n_v = mesh.nv
    verts = np.empty((n_v, 3), dtype=np.float64)
    for i in range(n_v):
        verts[i] = [mesh.vertices[i].point[k] for k in range(3)]
    bnd_names = list(mesh.GetBoundaries())
    src_idx = {i for i, n in enumerate(bnd_names) if n == "source"}
    snk_idx = {i for i, n in enumerate(bnd_names) if n == "sink"}
    tris, src_mask, snk_mask = [], [], []
    for el in mesh.Elements(BND):
        tris.append([int(v.nr) for v in el.vertices])
        src_mask.append(el.index in src_idx)
        snk_mask.append(el.index in snk_idx)
    tris = np.array(tris, dtype=np.int64)
    src_mask = np.array(src_mask, dtype=bool)
    snk_mask = np.array(snk_mask, dtype=bool)

    t0 = time.perf_counter()
    res = solve_inductance_source_sink_intree(
        verts, tris, src_mask, snk_mask,
        regular_quad_degree=regular_quad_degree,
        singular_n_q=singular_n_q)
    res["t_total"] = time.perf_counter() - t0
    res["mesh_nv"] = mesh.nv
    res["mesh_n_tris"] = len(tris)
    res["mesh_n_src_tris"] = int(src_mask.sum())
    res["mesh_n_snk_tris"] = int(snk_mask.sum())
    return res


# ============================================================
# Step 4: main
# ============================================================
def main():
    step_path = os.path.join(HERE, "loft_coil_b123d.step")
    print("=" * 70)
    print("build123d loft coil -- PEEC vs BEM-A cross-validation")
    print("=" * 70)
    print()
    print("[1] Building lofted coil with build123d "
          "(R=50mm, 8x6 rect, 355deg, 12 lofts)...")
    t0 = time.perf_counter()
    build_loft_coil_step(step_path)
    print(f"    saved: {step_path}  ({time.perf_counter()-t0:.1f}s)")
    print()

    print("[2] PEEC inductance (calc_peec_inductance.py, n_peri=16, "
          "50 kHz, Cu)...")
    t0 = time.perf_counter()
    peec = compute_peec_L(step_path)
    L_peec_nH = peec.get("L_coil_nH") or peec.get("L_nH") or (
        peec.get("L_coil_H", 0.0) * 1e9)
    print(f"    PEEC L = {L_peec_nH:.3f} nH  "
          f"({time.perf_counter()-t0:.1f}s)")
    print()

    print("[3] BEM-A intree saddle-point (maxh=12mm, RWG order=0)...")
    t0 = time.perf_counter()
    bem = compute_bem_a_L(step_path, maxh=0.012)
    L_bem_nH = bem["L"] * 1e9
    print(f"    BEM-A L = {L_bem_nH:.3f} nH  ({bem['t_total']:.1f}s)")
    print(f"    n_J = {bem['n_J']}, residual = {bem['residual']:.2e}")
    print()

    print("[4] Comparison")
    rel = (L_bem_nH - L_peec_nH) / L_peec_nH * 100.0
    print(f"    L_PEEC  = {L_peec_nH:.3f} nH")
    print(f"    L_BEM_A = {L_bem_nH:.3f} nH "
          f"(diff: {rel:+.2f}%)")
    print()
    print("    Note: ~10 % spread is expected at this mesh density.")
    print("    PEEC uses thin-skin perimeter filaments (n_peri=16);")
    print("    BEM-A solves the EFIE on the conductor surface in the")
    print("    PEC limit. The two converge with mesh + n_peri refinement.")
    print()


if __name__ == "__main__":
    main()
