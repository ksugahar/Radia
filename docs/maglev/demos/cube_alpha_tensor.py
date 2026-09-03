"""cube_alpha_tensor.py -- multi-port (matrix) Mixed-Galerkin admittance Y(s).

Frontier C + D of the radia-maglev mixed-Galerkin stack, demonstrated together
on a Cu cube with FOUR ports: the monopole f=1 plus the three centered dipole
(coordinate) drives {x-c, y-c, z-c}.

  D (multi-port matrix-CLN, multipole admittance):
    Each port is projected onto the SAME shared Dirichlet-Laplacian eigenbasis;
    the per-pole residue becomes the rank-1 matrix G_n = sigma V b_n b_n^T (the
    matrix-form Kameari / multi-port CLN structure, Matsuo 2017/2018c).  The
    matrix Y(s)_{pq} is a MULTIPOLE expansion of the conductor's scalar eddy
    response: the monopole port 1 <-> uniform external field, the dipoles
    x,y,z <-> field gradients.  Port 0 (monopole) reproduces the verified
    SCALAR eddy admittance exactly (its alpha_00(s) is the physical 0->1
    response in the mid-high band); the dipole ports add the multipole
    structure (isotropic for the symmetric cube).

  C (per-face / multi-port SIBC envelope):
    The scalar wetted-area tail K_SIBC = S sqrt(sigma/mu) generalizes to the
    surface MOMENT matrix K_mat[p,q] = sqrt(sigma/mu) integral_S f_p f_q dS,
    and the scalar edge coefficient c_1 to the per-edge moment matrix C1_mat
    (cad_edges.edge_moment_matrix, CAD-direct).  Both reduce to the scalar
    values for the monopole port.

The 4-port LTI is exported MIMO via simulink.export.build_state_space_mimo.

HONEST SCOPE: this is the multipole eddy-ADMITTANCE matrix of the SCALAR
diffusion model.  The monopole port and the scalar alpha(s) it produces are
validated (matches the scalar cube_alpha_sweep + the analytic high-f PEC
limit).  The dipole-port multipole admittances are structurally sane
(symmetric PSD residues, cube-isotropic) but are NOT validated against an
external multipole reference.  The physical 3x3 VECTOR polarizability tensor
(transverse m=1 + triaxial shape anisotropy) is the full 3D HCurl solve in
ellipsoid/ellipsoid_alpha_tensor_3d.py, not this.

Run:  python cube_alpha_tensor.py
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
from netgen.occ import Box, OCCGeometry, Pnt
from ngsolve import Mesh, TaskManager, x, y, z

from radia.maglev.mixed_galerkin import (
    bulk_foster_matrix_via_eigen,
    K_SIBC_matrix,
    K_SIBC_total,
    Y_matrix_mixed,
    alpha_matrix_from_Y,
    cad_topology_total_area,
    edge_moment_matrix,
)
from radia.maglev.simulink import build_state_space_mimo

from _validation_output import validation_output

HERE = os.path.dirname(os.path.abspath(__file__))
SIGMA_CU = 5.8e7
MU_0 = 4 * math.pi * 1e-7
L = 5e-3


def main():
    box = Box(Pnt(0, 0, 0), Pnt(L, L, L))
    box.mat("Cu")
    for f in box.faces:
        f.name = "outer"
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=L / 5))

    cen = L / 2.0
    drives_cf = [1.0, x - cen, y - cen, z - cen]               # monopole + 3 dipoles
    drives_fn = [lambda p: 1.0] + [(lambda i: (lambda p: p[i] - cen))(i)
                                   for i in range(3)]

    with TaskManager():
        lam, tau, G_n, V = bulk_foster_matrix_via_eigen(
            mesh, SIGMA_CU, MU_0, drives_cf, n_eigen=60)
        K_mat = K_SIBC_matrix(mesh, drives_cf, SIGMA_CU, MU_0)

    S_cad = cad_topology_total_area(box)
    C1_mat = edge_moment_matrix(box, drives_fn, MU_0)

    print(f"cube L = {L*1e3:.1f} mm, V = {V*1e9:.3f} mm^3, "
          f"S = {S_cad*1e6:.3f} mm^2,  4 ports (1, x, y, z)")
    print(f"K port-0 (monopole) = {K_mat[0,0]:.4e}  "
          f"(scalar K_SIBC_total = {K_SIBC_total(S_cad, SIGMA_CU, MU_0):.4e})")
    print(f"C1 port-0 (monopole) = {C1_mat[0,0]:.4e}")
    # bulk-mode completeness captured by the truncated Foster sum (rest -> tail)
    Yb0 = G_n.sum(axis=0)
    print(f"bulk Foster captures {Yb0[0,0]/(SIGMA_CU*V)*100:.0f}% of the "
          f"monopole DC drive; the SIBC tail completes the residual")
    print()

    # ---- monopole port: the PHYSICAL eddy admittance alpha_00(s)/V (0 -> 1) --
    print("monopole port alpha_00(s)/V  (physical, mid-high band):")
    print(f"  {'f (Hz)':>9}  {'Re':>9} {'Im':>9}    dipole|Y_11|  iso(spread)")
    sweep = {}
    for f in (1e3, 1e4, 1e5, 1e6, 1e9):
        s = 1j * 2 * math.pi * f
        Y = Y_matrix_mixed(s, lam, tau, G_n, K_mat, C1_mat)
        a = alpha_matrix_from_Y(Y, V, SIGMA_CU)
        a00 = a[0, 0] / V
        dip = np.abs(np.diag(Y)[1:])
        iso = (dip.max() - dip.min()) / dip.mean()
        print(f"  {f:9.1e}  {a00.real:+9.4f} {a00.imag:+9.4f}    "
              f"{dip.mean():.3e}   {iso*100:5.1f}%")
        sweep[f"{f:.0e}"] = {
            "alpha00_re_over_V": float(a00.real),
            "alpha00_im_over_V": float(a00.imag),
            "dipole_Yii_mean": float(dip.mean()),
            "dipole_isotropy_spread": float(iso),
        }

    # ---- MIMO LTI export + transfer-function cross-check ---------------------
    A_, B_, C_, D_, n_f, n_w, n_i = build_state_space_mimo(
        G_n, tau, V, SIGMA_CU, K_mat, C1_mat, n_warburg_rungs=30)
    n_states = A_.shape[0]
    eye = np.eye(n_states)
    worst = 0.0
    for f in (1e3, 1e4, 1e5):
        s = 1j * 2 * math.pi * f
        H = C_ @ np.linalg.solve(s * eye - A_, B_) + D_
        Ycf = Y_matrix_mixed(s, lam, tau, G_n, K_mat, C1_mat)
        worst = max(worst, np.linalg.norm(H - Ycf) / np.linalg.norm(Ycf))
    print(f"\nMIMO LTI: {n_states} states ({n_f} Foster + {n_w} Warburg "
          f"+ {n_i} integrator), 4 in / 4 out")
    print(f"transfer-matrix vs Y_matrix_mixed: worst {worst*100:.3f}% over "
          f"1e3-1e5 Hz")

    out = {
        "geometry": "Cu cube", "L_m": L, "sigma": SIGMA_CU, "mu": MU_0,
        "V_m3": V, "S_m2": S_cad, "n_foster": len(lam),
        "ports": ["monopole_1", "dipole_x", "dipole_y", "dipole_z"],
        "K_port0": float(K_mat[0, 0]),
        "C1_port0": float(C1_mat[0, 0]),
        "bulk_monopole_dc_fraction": float(Yb0[0, 0] / (SIGMA_CU * V)),
        "mimo_states": {"total": int(n_states), "foster": int(n_f),
                        "warburg": int(n_w), "integrator": int(n_i)},
        "mimo_transfer_worst_relerr": float(worst),
        "sweep": sweep,
        "scope_note": ("multipole eddy-admittance matrix of the SCALAR "
                       "diffusion model; monopole port = validated scalar "
                       "alpha(s); physical vector tensor is in "
                       "ellipsoid/ellipsoid_alpha_tensor_3d.py"),
    }
    output = validation_output("cube_alpha_tensor_results.json", HERE)
    with open(output, "w") as fp:
        json.dump(out, fp, indent=2, allow_nan=False)
    print(f"\nwrote {output}")


if __name__ == "__main__":
    main()
