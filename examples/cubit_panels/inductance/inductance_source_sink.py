"""
BEM inductance via source/sink saddle point EFIE.

Standalone test script using OCC gapped torus.
Core solver: src/radia/bem_inductance.py

Usage:
    python inductance_source_sink.py
"""

import math
import os
import sys
import numpy as np

# Import core solver from src/radia/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from bem_inductance import compute_inductance_source_sink, MU_0
from ngsolve import *
from ngsolve import TaskManager


def make_gapped_torus_mesh(R, a, gap_deg=5, maxh=None):
    """Create OCC gapped torus surface mesh with source/sink labels.

    Uses Glue(faces) to create a surface-only mesh (no volume elements),
    which gives correct HDivSurface DOF count for BEM.

    Args:
        R: Major radius [m]
        a: Minor radius [m]
        gap_deg: Gap angle [degrees]
        maxh: Max mesh size (default: a/2)

    Returns:
        NGSolve Mesh with boundary labels "conductor", "source", "sink"
    """
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters

    if maxh is None:
        maxh = a / 2

    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)

    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"

    torus_surf = Glue(torus.faces)
    geo = OCCGeometry(torus_surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=2.0,
                             segmentsperedge=2))
    mesh = Mesh(ngmesh)
    with TaskManager():
        mesh.Curve(2)

        return mesh


if __name__ == "__main__":
    R, a = 0.05, 0.005
    L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)

    print(f"Torus: R = {R*1e3:.0f} mm, a = {a*1e3:.0f} mm")
    print(f"Neumann formula: L = {L_neumann*1e9:.2f} nH")
    print()

    mesh = make_gapped_torus_mesh(R, a, gap_deg=5, maxh=a)
    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne_edges = sum(1 for e in mesh.edges)

    bnd_names = list(set(mesh.GetBoundaries()))
    print(f"Boundary labels: {bnd_names}")
    print(f"Mesh: {nse} faces, {nv} vertices, {ne_edges} edges")

    result = compute_inductance_source_sink(mesh)

    if 'error' in result:
        print(f"ERROR: {result['error']}")
    else:
        L = result['L']
        err = (L - L_neumann) / L_neumann * 100
        print(f"\nResult: L = {L*1e9:.4f} nH  (err = {err:+.2f}%)")
        print(f"Neumann reference: {L_neumann*1e9:.2f} nH")
        print(f"BEM assembly: {result['t_assembly']}s, "
              f"LU solve: {result['t_solve']}s, "
              f"total: {result['t_total']}s")
