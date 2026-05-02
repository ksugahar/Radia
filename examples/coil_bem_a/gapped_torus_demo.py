"""BEM-A coil self-inductance demo on a gapped torus (Weggler stabilized
EFIE = Lucy decomposition).

Geometry: 355deg gapped torus, R_maj = 30 mm, r_min = 3 mm.
Reference: PEEC golden L = 85.10 nH (band [80, 92]).

This demo runs the saddle-point EFIE formulation:

    [ SL  D^T ] [ J ]   [ 0 ]
    [ D   0   ] [ p ] = [ g ]

at HDivSurface(0) (= RWG, lowest order), gives L = mu_0 * J^T SL J.

See README.md in this folder for the math and the rationale for choosing
BEM-A over scalar BEM-SIBC.

Run::

    python gapped_torus_demo.py
"""
from __future__ import annotations

import math
import os
import sys
import time

# Make the reference solver importable.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(
    ROOT, "examples", "induction_heating", "bem_reference"))

from bem_inductance import compute_inductance_source_sink, MU_0


def make_gapped_torus_surface_mesh(R_major=0.030, r_minor=0.003,
                                    gap_deg=5.0, maxh=None):
    """OCC gapped-torus surface-only mesh with source / sink labelled."""
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    if maxh is None:
        maxh = r_minor

    wp = WorkPlane(Axes(p=Pnt(R_major, 0, 0), n=Dir(0, 1, 0),
                         h=Dir(0, 0, 1)))
    profile = wp.Circle(r_minor).Face()
    sweep_deg = 360.0 - gap_deg
    torus = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), sweep_deg)

    for f in torus.faces:
        f.name = "conductor"
    # The two flat caps are the smallest-area faces (= pi*r_minor^2 each).
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"

    surf = Glue(torus.faces)
    geo = OCCGeometry(surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=2.0,
                              segmentsperedge=2))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    return mesh


def main():
    R_maj = 0.030
    r_min = 0.003
    gap_deg = 5.0
    L_neumann_nH = (MU_0 * R_maj
                    * (math.log(8 * R_maj / r_min) - 2.0)) * 1e9
    L_PEEC_nH = 85.10

    print("=" * 70)
    print("BEM-A coil self-inductance via Weggler stabilized EFIE")
    print("=" * 70)
    print(f"Geometry: gapped torus  R = {R_maj*1e3:.0f} mm, "
          f"r = {r_min*1e3:.0f} mm, gap = {gap_deg:.0f} deg")
    print()
    print(f"Reference values:")
    print(f"  Maxwell analytical (no gap correction): "
          f"L = {L_neumann_nH:.2f} nH")
    print(f"  PEEC golden (band [80, 92]):            "
          f"L = {L_PEEC_nH:.2f} nH")
    print()

    print("[1] Building surface mesh (OCC + Glue)...")
    t0 = time.perf_counter()
    mesh = make_gapped_torus_surface_mesh(R_maj, r_min, gap_deg, maxh=r_min)
    t_mesh = time.perf_counter() - t0
    from ngsolve import BND
    nse = mesh.GetNE(BND)
    print(f"    {nse} BND triangles, {mesh.nv} vertices, t = {t_mesh:.2f}s")
    print()

    print("[2] Solving saddle-point EFIE (HDivSurface order=0, dense LU)...")
    print("    NOTE: ngsolve.bem.LaplaceSL assembly is slow "
          "(~100s at N_J ~ 5000)")
    t0 = time.perf_counter()
    res = compute_inductance_source_sink(
        mesh, source_label="source", sink_label="sink",
        fes_order=0, solver="lu")
    t_solve = time.perf_counter() - t0

    if "error" in res:
        print(f"    ERROR: {res['error']}")
        return 1

    print(f"    n_J = {res['n_J']} (HDivSurface DOFs)")
    print(f"    n_f = {res['n_f']} (SurfaceL2 DOFs)")
    print(f"    A_source = {res['A_source']*1e6:.3f} mm^2  "
          f"(expected pi*r^2 = {math.pi*r_min**2*1e6:.3f} mm^2)")
    print(f"    A_sink   = {res['A_sink']*1e6:.3f} mm^2")
    print(f"    div(J) residual = {res['residual']:.2e}")
    print(f"    t_assembly = {res['t_assembly']:.1f}s, "
          f"t_solve = {res['t_solve']:.2f}s, total = {t_solve:.1f}s")
    print()

    L_nH = res["L"] * 1e9
    rel_err_pct = (L_nH - L_PEEC_nH) / L_PEEC_nH * 100.0

    print("[3] Result")
    print(f"    L_BEM_A  = {L_nH:.3f} nH")
    print(f"    L_PEEC   = {L_PEEC_nH:.2f} nH "
          f"(rel diff: {rel_err_pct:+.2f}%)")
    print(f"    Maxwell  = {L_neumann_nH:.2f} nH "
          f"(BEM should land between PEEC and Maxwell)")

    in_PEEC_band = 80.0 <= L_nH <= 92.0
    print()
    print(f"    PEEC golden hard band [80, 92] : "
          f"{'PASS' if in_PEEC_band else 'FAIL'}")
    return 0 if in_PEEC_band else 2


if __name__ == "__main__":
    sys.exit(main())
