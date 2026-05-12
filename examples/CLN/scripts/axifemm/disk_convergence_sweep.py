"""
disk_convergence_sweep.py — Probe axifemm convergence to v22b (208.32 us).

Two error sources to characterize:
  1. Mesh resolution (maxh_disk, maxh_air)
  2. Truncation domain size (R_AIR, Z_AIR) — currently uses Dirichlet box,
     not Kelvin transform.

Strategy: vary one at a time, report tau_1 and check ratio to v22b 208.32 us.
"""

from __future__ import annotations

import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from axifemm_core import assemble, MU0
from sigma_mass import assemble_sigma_mass

R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
V22B_TAU1 = 208.32  # microseconds

TARGET_RATIO_OK = 1.01   # within 1%


def build_mesh(R_air, Z_air, maxh_disk, maxh_air):
    from ngsolve import Mesh
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y

    air_box = MoveTo(0, -Z_air).Rectangle(R_air, 2 * Z_air).Face()
    disk = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    disk.faces.name = "conductor"
    disk.maxh = maxh_disk
    air_box.edges.Max(X).name = "right"
    air_box.edges.Min(X).name = "axis"
    air_box.edges.Max(Y).name = "top"
    air_box.edges.Min(Y).name = "bot"
    shape = Glue([air_box, disk])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air))


def extract(mesh):
    ng = mesh.ngmesh
    nodes = np.array([(p.p[0], p.p[1]) for p in ng.Points()])
    tri = np.array([[v.nr - 1 for v in el.vertices] for el in ng.Elements2D()])
    mats = np.array([ng.GetMaterial(el.index) for el in ng.Elements2D()])
    return nodes, tri, mats


def boundary_idx(nodes, R_air, Z_air, tol=1e-6):
    r, z = nodes[:, 0], nodes[:, 1]
    bnd = (r < tol) | (r > R_air - tol) | (z > Z_air - tol) | (z < -Z_air + tol)
    return np.where(bnd)[0]


def solve_tau1(R_air, Z_air, maxh_disk, maxh_air, n_eigs=4):
    mesh = build_mesh(R_air, Z_air, maxh_disk, maxh_air)
    nodes, tri, mats = extract(mesh)
    Nelem = tri.shape[0]
    mu_r = np.full(Nelem, MU0)
    mu_z = np.full(Nelem, MU0)
    sigma = np.zeros(Nelem)
    sigma[mats == "conductor"] = SIGMA_CU

    K, _ = assemble(nodes, tri, mu_r, mu_z)
    M_sigma = assemble_sigma_mass(nodes, tri, sigma)

    bnd = boundary_idx(nodes, R_air, Z_air)
    free = np.setdiff1d(np.arange(nodes.shape[0]), bnd)

    K_red = -np.pi * K[free[:, None], free[None, :]]
    M_red = M_sigma[free[:, None], free[None, :]]
    K_red = (K_red + K_red.T) / 2.0
    M_red = (M_red + M_red.T) / 2.0

    eigs, _ = spla.eigsh(K_red, k=n_eigs, M=M_red, sigma=0.0, which="LM",
                          tol=1e-10, maxiter=3000)
    taus = np.sort(1.0 / eigs)[::-1] * 1e6  # us, descending
    return taus, nodes.shape[0], Nelem, len(free)


def main():
    print(f"=== axifemm convergence sweep vs v22b tau_1 = {V22B_TAU1} us ===\n")

    # Sweep 1: truncation domain (mesh fixed at moderate)
    print("[1] Truncation domain sweep (maxh_disk=0.2mm, maxh_air=10mm)")
    print(f"  {'R_air mm':>10} {'Z_air mm':>10} {'Nnode':>8} {'Nelem':>8} {'tau_1 us':>12} {'ratio':>10}")
    for R_air, Z_air in [(50e-3, 50e-3), (100e-3, 100e-3), (200e-3, 200e-3),
                         (500e-3, 500e-3), (1.0, 1.0)]:
        try:
            taus, Nn, Ne, _ = solve_tau1(R_air, Z_air, 0.2e-3, 10e-3)
            t1 = taus[0]
            ratio = t1 / V22B_TAU1
            print(f"  {R_air*1e3:>10.0f} {Z_air*1e3:>10.0f} {Nn:>8} {Ne:>8} "
                  f"{t1:>12.4f} {ratio:>10.4f}")
        except Exception as e:
            print(f"  {R_air*1e3:>10.0f} {Z_air*1e3:>10.0f} ERROR: {e}")

    # Sweep 2: mesh refinement (truncation fixed large)
    print("\n[2] Mesh refinement (R_air=500mm, Z_air=500mm)")
    print(f"  {'maxh_disk':>10} {'maxh_air':>10} {'Nnode':>8} {'Nelem':>8} {'tau_1 us':>12} {'ratio':>10}")
    for h_disk, h_air in [(0.4e-3, 20e-3), (0.2e-3, 10e-3), (0.1e-3, 5e-3),
                          (0.05e-3, 5e-3)]:
        try:
            taus, Nn, Ne, _ = solve_tau1(0.5, 0.5, h_disk, h_air)
            t1 = taus[0]
            ratio = t1 / V22B_TAU1
            print(f"  {h_disk*1e3:>10.3f} {h_air*1e3:>10.1f} {Nn:>8} {Ne:>8} "
                  f"{t1:>12.4f} {ratio:>10.4f}")
        except Exception as e:
            print(f"  {h_disk*1e3:>10.3f} ERROR: {e}")


if __name__ == "__main__":
    main()
