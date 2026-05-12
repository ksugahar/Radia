"""
nmr_validate.py — Reproduce the FEMM NMR axisymmetric problem with axifemm_core
                  and compare against FEMM.mat reference data.

Geometry (matches axisymmetric_mixed.py and FEMM.fem):
  - Domain: rectangle r in [0, 70mm], z in [-50mm, 50mm]
  - Two ferrite magnets:
      mag1: r in [0, 40mm], z in [-20mm, -10mm], M = +z, H_c = 795774.7 A/m
      mag2: r in [0, 40mm], z in [+10mm, +20mm], M = +z, H_c = 795774.7 A/m
  - Air everywhere else.
  - Dirichlet phi = 0 on all outer boundaries (top, bottom, right) and on axis.

Comparison: sample (B_z, A_phi) along z = 0, r in [0, 70mm], 1000 points.
"""

from __future__ import annotations

import os
import sys
import numpy as np
import scipy.io
import scipy.sparse as sp

# Adjust path
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from axifemm_core import (
    assemble,
    apply_dirichlet,
    solve,
    field_in_element,
    find_element,
    MU0,
)


def build_nmr_mesh(maxh: float = 0.001):
    """Generate the NMR axisymmetric mesh with NGSolve OCC, return arrays."""
    from ngsolve import Mesh
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y

    dom = MoveTo(0, -0.05).Rectangle(0.07, 0.10).Face()
    dom.edges.Max(X).name = "right"
    dom.edges.Min(X).name = "left"
    dom.edges.Max(Y).name = "top"
    dom.edges.Min(Y).name = "bot"

    mag1 = MoveTo(0, -0.02).Rectangle(0.04, 0.01).Face()
    mag2 = MoveTo(0, 0.01).Rectangle(0.04, 0.01).Face()
    mag1.faces.name = "mag1"
    mag2.faces.name = "mag2"

    shape = Glue([dom, mag1, mag2])
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh))

    # Extract nodes and triangles from netgen mesh
    ng_mesh = mesh.ngmesh
    points = np.array([(p.p[0], p.p[1]) for p in ng_mesh.Points()])  # (Nnode, 2): (r, z)
    nodes = points

    triangles = []
    materials = []  # per-element material name
    for el in ng_mesh.Elements2D():
        verts = [v.nr - 1 for v in el.vertices]  # 0-based
        triangles.append(verts)
        # material name lookup
        mat_idx = el.index
        mat_name = ng_mesh.GetMaterial(mat_idx)
        materials.append(mat_name)

    triangles = np.array(triangles)
    materials = np.array(materials)

    print(f"  Mesh: {nodes.shape[0]} nodes, {triangles.shape[0]} triangles")
    print(f"  Materials present: {set(materials)}")
    return nodes, triangles, materials, mesh


def boundary_node_indices(nodes, tol=1e-6):
    """Identify nodes on boundary (axis, right, top, bottom) for Dirichlet."""
    r = nodes[:, 0]
    z = nodes[:, 1]
    on_axis = r < tol
    on_right = r > 0.07 - tol
    on_top = z > 0.05 - tol
    on_bot = z < -0.05 + tol
    bnd = on_axis | on_right | on_top | on_bot
    return np.where(bnd)[0]


def main():
    print("=== Phase 1a: Triangle P1 axifemm vs FEMM NMR ===")

    # Mesh
    print("\n[1] Building mesh...")
    nodes, triangles, materials, ng_mesh = build_nmr_mesh(maxh=0.002)
    Nelem = triangles.shape[0]

    # Material setup
    print("\n[2] Setting up material properties...")
    mu_r_arr = np.full(Nelem, MU0)
    mu_z_arr = np.full(Nelem, MU0)
    H_c_arr = np.zeros(Nelem)
    M_dir_arr = np.zeros(Nelem)

    H_C_MAG = 795774.71545947669  # A/m, from FEMM.fem
    for ie in range(Nelem):
        mat = materials[ie]
        if mat in ("mag1", "mag2"):
            H_c_arr[ie] = H_C_MAG
            M_dir_arr[ie] = 90.0  # +z direction (matches FEMM block label MagDir=90)

    print(f"  Magnetized elements: {(H_c_arr != 0).sum()}/{Nelem}")

    # Assemble
    print("\n[3] Assembling stiffness matrix and RHS...")
    K, b = assemble(nodes, triangles, mu_r_arr, mu_z_arr,
                    mat_H_c=H_c_arr, mat_M_dir_deg=M_dir_arr)
    print(f"  K: shape={K.shape}, nnz={K.nnz}")
    print(f"  b: range=[{b.min():.4e}, {b.max():.4e}]")

    # Apply Dirichlet
    print("\n[4] Applying Dirichlet BC (phi = 0 on boundary)...")
    bnd = boundary_node_indices(nodes)
    print(f"  Dirichlet nodes: {len(bnd)}")
    K_bc, b_bc = apply_dirichlet(K, b, bnd, np.zeros(len(bnd)))

    # Solve
    print("\n[5] Solving linear system...")
    V = solve(K_bc, b_bc)
    print(f"  V range=[{V.min():.4e}, {V.max():.4e}]")

    # Convert V -> phi at each node.
    # Empirically the V we solve for has units of A_phi directly (Wb/m), not
    # A_phi/mu_0 — the K matrix and RHS scale accordingly in pure SI.
    # So phi_j = 2*pi*r_j*V_j with V_j = A_phi at node j.
    phi_nodes = 2.0 * np.pi * nodes[:, 0] * V

    # Sample along z=0 axis
    print("\n[6] Sampling field along z=0...")
    rr = np.linspace(0.0, 0.07, 1000)
    A_axifemm = np.zeros_like(rr)
    Bz_axifemm = np.zeros_like(rr)
    Br_axifemm = np.zeros_like(rr)
    failed = 0
    for i, r in enumerate(rr):
        # Find element containing (r, 0)
        if r < 1e-9:
            r_eval = 1e-7
        else:
            r_eval = r
        ie = find_element(nodes, triangles, r_eval, 0.0)
        if ie < 0:
            failed += 1
            continue
        idx = triangles[ie]
        rn = nodes[idx, 0]
        zn = nodes[idx, 1]
        phi_e = phi_nodes[idx]
        A, Br, Bz = field_in_element(rn, zn, phi_e, r_eval, 0.0)
        A_axifemm[i] = A
        Br_axifemm[i] = Br
        Bz_axifemm[i] = Bz
    print(f"  Sampled {len(rr) - failed}/{len(rr)} points (failures: {failed})")

    # Load reference
    print("\n[7] Loading reference data...")
    femm = scipy.io.loadmat("W:/00_CAE/NGSolve/01_菅原/2024_08_01_軸対称/NMR/FEMM.mat")
    r_femm = femm["r"].flatten()       # mm
    Bz_femm = femm["Bz"].flatten()     # T
    At_femm = femm["At"].flatten()     # 2*pi*r*A_phi (Wb)

    mixed = scipy.io.loadmat("W:/00_CAE/NGSolve/01_菅原/2024_08_01_軸対称/NMR/axisymmetric_mixed.mat")
    r_mixed = mixed["r"].flatten() * 1000  # m -> mm
    Bz_mixed = mixed["Bz0"].flatten()
    At_mixed = mixed["At"].flatten()    # r*A_phi (NOT 2*pi*r*A_phi)

    # Save our results in the same format
    out = {
        "r": rr,
        "Bz": Bz_axifemm,
        "Br": Br_axifemm,
        "A_phi": A_axifemm,
        "rA_phi": rr * A_axifemm,  # to compare with mixed's At = r*A_phi
        "method": "axifemm Phase 1a (Triangle {1, r^2, z}, FEMM-style)",
    }
    out_path = os.path.join(HERE, "nmr_axifemm_p1a.mat")
    scipy.io.savemat(out_path, out)
    print(f"\n  Saved: {out_path}")

    # Numerical comparison at a few representative r values
    print("\n[8] Numerical comparison at sample points:")
    print(f"  {'r (mm)':>8} {'B_z FEMM':>14} {'B_z mixed':>14} {'B_z axifemm':>16} {'ratio':>10}")
    for r_mm in [0.5, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0]:
        i_femm = np.argmin(np.abs(r_femm - r_mm))
        i_mixed = np.argmin(np.abs(r_mixed - r_mm))
        i_ax = np.argmin(np.abs(rr * 1000 - r_mm))
        ratio = (Bz_femm[i_femm] / Bz_axifemm[i_ax]) if abs(Bz_axifemm[i_ax]) > 1e-30 else float("nan")
        print(f"  {r_mm:>8.2f} {Bz_femm[i_femm]:>14.6e} {Bz_mixed[i_mixed]:>14.6e} {Bz_axifemm[i_ax]:>16.6e} {ratio:>10.3e}")

    # Try a quick matplotlib plot if available
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))
        ax1.plot(r_femm, Bz_femm, "k-", label="FEMM, p=2")
        ax1.plot(r_mixed, Bz_mixed, "r--", label="NGSolve mixed, p=2")
        ax1.plot(rr * 1000, Bz_axifemm, "b-", label="axifemm Phase 1a")
        ax1.set_xlabel("r (mm)"); ax1.set_ylabel("B_z (T)")
        ax1.set_xlim(0, 10); ax1.set_ylim(0.16, 0.164)
        ax1.legend(); ax1.grid(True)

        ax2.plot(r_femm, At_femm / (2 * np.pi), "k-", label="FEMM, p=2")
        ax2.plot(r_mixed, At_mixed, "r--", label="NGSolve mixed, p=2")
        ax2.plot(rr * 1000, rr * A_axifemm, "b-", label="axifemm Phase 1a")
        ax2.set_xlabel("r (mm)"); ax2.set_ylabel("r*A_phi (Wb)")
        ax2.set_xlim(0, 10)
        ax2.legend(); ax2.grid(True)

        out_png = os.path.join(HERE, "nmr_axifemm_p1a.png")
        plt.tight_layout(); plt.savefig(out_png, dpi=120)
        print(f"\n  Plot saved: {out_png}")
    except ImportError:
        print("  (matplotlib not available, skipping plot)")


if __name__ == "__main__":
    main()
