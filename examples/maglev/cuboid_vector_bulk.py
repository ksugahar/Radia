"""cuboid_vector_bulk.py -- vector (HCurl) eddy-current Foster spectrum.

Frontier B of the radia-maglev mixed-Galerkin stack: the de-Rham VECTOR
partner of the scalar bulk (HDiv = demagnetisation, HCurl = eddy current).
Instead of the scalar diffusion Laplacian, this solves the curl-curl
generalized eigenproblem S w = lam M w (S = curl-curl/mu, M = sigma mass) on
an HCurl(nograds) space with A x n = 0 on the conductor surface, and projects
the three uniform-field drives A_ext_k = (1/2) e_k x (r-r_c) onto the modes.

A non-cubic Cu box (5 x 2 x 1 mm) has three DISTINCT leading eddy time
constants -- the shape split: a field along z drives currents in the 5x2 mm
cross-section (largest tau), along x in the 2x1 mm cross-section (smallest).
This script reproduces the analytic interior-PEC TE-mode tau

    tau_k = mu sigma / (pi^2 (1/La^2 + 1/Lb^2))

across two mesh densities, demonstrating convergence to <2.4% at h = a/28.

SCOPE (honest): these are the INTERIOR-PEC eddy modes (A x n = 0), the vector
analog of the scalar H1-Dirichlet bulk -- a MODEL, not the physical
exterior-matched (free-decay / Stoll) spectrum.  The verified PHYSICAL
polarizability tensor (with the exterior reaction field) is the per-frequency
3D HCurl solve in ellipsoid/ellipsoid_alpha_tensor_3d.py.  The method is
mesh-hungry and the tree-cotree gauge / leading-mode selection is
resolution-sensitive; see mixed_galerkin/vector_bulk.py.

Run:  python cuboid_vector_bulk.py    (h=0.25 then h=0.18; ~1-2 min total)
"""
from __future__ import annotations

import json
import math
import os

from netgen.occ import Box, OCCGeometry, Pnt
from ngsolve import Mesh, TaskManager

from radia.maglev.mixed_galerkin import bulk_foster_vector_via_eigen

HERE = os.path.dirname(os.path.abspath(__file__))
MU_0 = 4 * math.pi * 1e-7
SIGMA_CU = 5.8e7
AX, AY, AZ = 5e-3, 2e-3, 1e-3


def tau_te(L1, L2):
    """Analytic interior-PEC TE-mode time constant."""
    return MU_0 * SIGMA_CU / (math.pi**2 * (1.0 / L1**2 + 1.0 / L2**2))


def main():
    tau_ref = {0: tau_te(AY, AZ), 1: tau_te(AX, AZ), 2: tau_te(AX, AY)}
    nm = {0: "x", 1: "y", 2: "z"}
    print("Cu box 5x2x1 mm -- vector (HCurl) eddy Foster vs analytic TE tau")
    print(f"  analytic: tau_x={tau_ref[0]*1e6:.3f}  tau_y={tau_ref[1]*1e6:.3f}  "
          f"tau_z={tau_ref[2]*1e6:.3f} us\n")

    runs = {}
    for h in (0.25e-3, 0.18e-3):
        box = Box(Pnt(-AX/2, -AY/2, -AZ/2), Pnt(AX/2, AY/2, AZ/2))
        box.mat("conductor").bc("conductor_surface")
        box.maxh = h
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=h))
        with TaskManager():
            lam, tau_n, G_n, V, lead = bulk_foster_vector_via_eigen(
                mesh, SIGMA_CU, MU_0, n_per_dir=12, order=2,
                conductor_bnd="conductor_surface")
        print(f"h = {h*1e3:.2f} mm   ne={mesh.ne}  Nmodes={len(lam)}")
        rec = {}
        for k in range(3):
            err = (lead[k] - tau_ref[k]) / tau_ref[k] * 100
            print(f"   tau_{nm[k]} = {lead[k]*1e6:8.4f} us   "
                  f"(analytic {tau_ref[k]*1e6:7.4f},  err {err:+6.2f}%)")
            rec[nm[k]] = {"tau_us": lead[k]*1e6, "err_pct": err}
        print(f"   shape-split ordering tau_z>tau_y>tau_x: "
              f"{lead[2] > lead[1] > lead[0]}\n")
        runs[f"h_{h*1e3:.2f}mm"] = {"ne": mesh.ne, "n_modes": int(len(lam)),
                                    "leading_tau": rec}

    out = {
        "geometry": "Cu cuboid 5x2x1 mm", "sigma": SIGMA_CU, "mu": MU_0,
        "dims_m": [AX, AY, AZ],
        "analytic_tau_te_us": {nm[k]: tau_ref[k]*1e6 for k in range(3)},
        "runs": runs,
        "scope_note": ("interior-PEC eddy modes (A x n = 0), the de-Rham HCurl "
                       "partner of the scalar H1 bulk; physical exterior-matched "
                       "tensor is ellipsoid/ellipsoid_alpha_tensor_3d.py"),
    }
    with open(os.path.join(HERE, "cuboid_vector_bulk_results.json"), "w") as fp:
        json.dump(out, fp, indent=2)
    print("wrote cuboid_vector_bulk_results.json")


if __name__ == "__main__":
    main()
