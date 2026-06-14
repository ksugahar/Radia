"""cube_alpha_sweep.py -- Generate a 5mm Cu cube .vol and sweep alpha(s).

Demonstrates the radia.levitation.mixed_galerkin API on the simplest case.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from netgen.occ import Box, OCCGeometry, Pnt
from ngsolve import Mesh, TaskManager

from radia.levitation.mixed_galerkin import (
    bulk_foster_via_eigen,
    K_SIBC_total,
    measure_total_area,
    Y_mixed,
    alpha_from_Y,
    cad_topology_c1,
    cad_topology_total_area,
)


def main():
    SIGMA_CU = 5.8e7
    MU_0 = 4 * math.pi * 1e-7
    L = 5e-3

    # Build cube via OCC
    box = Box(Pnt(0, 0, 0), Pnt(L, L, L))
    box.mat("Cu")
    for f in box.faces:
        f.name = "outer"
    geo = OCCGeometry(box)
    ng = geo.GenerateMesh(maxh=L/6)
    mesh = Mesh(ng)
    ng.Save("cube_5mm.vol")

    # CAD-direct c_1 + area (mesh-independent, IGA-style -- the recommended,
    # dependency-free route; exact for polyhedra).
    c1_cad, L_total_cad, n_edges_cad = cad_topology_c1(box, MU_0)
    S_cad = cad_topology_total_area(box)
    print(f"CAD-direct: {n_edges_cad} edges, total L = {L_total_cad*1e3:.3f} mm, "
          f"S = {S_cad*1e6:.3f} mm^2")
    print(f"            c_1 = {c1_cad:.4e}")
    print(f"  (cube closed form c_1 = -16(3L)/(pi mu) = {-16*3*L/(math.pi*MU_0):.4e})")

    # Cross-check: CAD area vs mesh-integrated area (both must agree).
    with TaskManager():
        lam, tau, g_n, V = bulk_foster_via_eigen(mesh, SIGMA_CU, MU_0, n_eigen=80)
        S_mesh = measure_total_area(mesh)
    print(f"Mesh-integrated area S = {S_mesh*1e6:.3f} mm^2 "
          f"(CAD-direct {S_cad*1e6:.3f} mm^2)")
    K_SIBC = K_SIBC_total(S_cad, SIGMA_CU, MU_0)
    print()

    print(f"V = {V*1e9:.3f} mm^3, K_SIBC = {K_SIBC:.4e}")
    print(f"Foster modes: {len(lam)}")
    print()
    print(f"  {'f (Hz)':>10}  {'|alpha|/V':>12}  {'Re(alpha)/V':>14}  {'-Im(alpha)/V':>14}")
    for f in [1e3, 1e4, 1e5, 1e6, 1e9]:
        s = 1j * 2 * math.pi * f
        Y = Y_mixed(s, lam, tau, g_n, K_SIBC, c1_cad)   # use CAD-direct c_1
        a = alpha_from_Y(Y, V, SIGMA_CU)
        print(f"  {f:10.2e}  {abs(a)/V:12.4f}  {a.real/V:+14.4e}  {-a.imag/V:+14.4e}")


if __name__ == "__main__":
    main()
