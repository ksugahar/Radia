"""Coil field synthesis with the (ACA+)+TSVD least-norm solver.

Stream-function-style coil design: place N thin current loops (the basis) on a
plane and solve for their currents `phi` so the axial field over a target
region matches a desired profile `B`.  The field-coupling matrix

    A(i, j) = Bz at observation point i  from unit-current loop j

is built from Radia's OWN Biot-Savart (radia.ObjFlmCur + radia.Fld) via
`radia_field_kernel` -- no coil-specific field code here.  The underdetermined
system `A phi = B` (M < N) is regularised with the truncated-SVD pseudo-inverse,
accelerated by ACA+ recompression (HACApK) so only the small factors are SVD'd.

Run:  python demo_coil_field_synthesis.py
"""
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

import radia as rad
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel


def build_loops(n_side, half_extent, loop_half_side, z_plane, current=1.0):
    """N = n_side^2 square filament loops on the z = z_plane plane."""
    xs = np.linspace(-half_extent, half_extent, n_side)
    ys = np.linspace(-half_extent, half_extent, n_side)
    handles = []
    a = loop_half_side
    for cx in xs:
        for cy in ys:
            pts = [[cx - a, cy - a, z_plane],
                   [cx + a, cy - a, z_plane],
                   [cx + a, cy + a, z_plane],
                   [cx - a, cy + a, z_plane],
                   [cx - a, cy - a, z_plane]]
            handles.append(rad.ObjFlmCur(pts, current))
    return handles


def main():
    rad.UtiDelAll()
    try:
        # --- geometry --------------------------------------------------------
        n_side = 8                                   # N = 64 basis loops
        loops = build_loops(n_side, half_extent=0.10,
                            loop_half_side=0.012, z_plane=0.0)
        N = len(loops)

        m_side = 5                                   # M = 25 field points
        oxs = np.linspace(-0.06, 0.06, m_side)
        oys = np.linspace(-0.06, 0.06, m_side)
        ox, oy = np.meshgrid(oxs, oys, indexing="ij")
        obs = np.column_stack([ox.ravel(), oy.ravel(),
                               np.full(ox.size, 0.05)])
        M = obs.shape[0]

        # --- matrix entry from Radia's Biot-Savart ---------------------------
        entry = radia_field_kernel(obs, loops, component=2, field="b")

        # dense A only for residual evaluation / reporting (not needed by ACA+)
        A = np.array([[entry(i, j) for j in range(N)] for i in range(M)])

        # --- (ACA+)+TSVD factorization ---------------------------------------
        res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                       aca_eps=1.0e-8)

        print("=== Coil field synthesis: (ACA+)+TSVD ===")
        print(f"M (field points) = {M},  N (loops) = {N}")
        print(f"ACA+ rank k_aca  = {res.k_aca}  (of min(M,N) = {min(M, N)})")
        print(f"singular values S[:6] = "
              f"{np.array2string(res.S[:6], precision=4)}")

        # --- target: uniform Bz over the region ------------------------------
        B_target = np.full(M, float(np.mean(np.abs(A))))  # ~ scale of the field

        # --- TSVD L-curve: residual vs number of modes -----------------------
        print("\n  modes |  ||A phi - B|| / ||B||  |   ||phi||")
        print("  ------+-------------------------+-------------")
        ks = sorted(set([1, 2, 3, 5, 8, res.k_aca]))
        ks = [k for k in ks if 1 <= k <= res.k_aca]
        resid, solnorm = {}, {}
        for k in ks:
            phi = pseudo_inverse_solve(res, B_target, k_mode=k)
            r = np.linalg.norm(A @ phi - B_target) / np.linalg.norm(B_target)
            resid[k] = r
            solnorm[k] = float(np.linalg.norm(phi))
            print(f"  {k:5d} |        {r:.3e}        |  {solnorm[k]:.3e}")

        # --- optional plot ---------------------------------------------------
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams["pdf.fonttype"] = 42
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
            ax1.semilogy(np.arange(1, res.k_aca + 1), res.S[:res.k_aca], "o-")
            ax1.set_xlabel("mode index"); ax1.set_ylabel("singular value")
            ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=10)
            ax1.grid(True, which="both", alpha=0.3)
            kk = sorted(resid)
            ax2.semilogy(kk, [resid[k] for k in kk], "s-", label="residual")
            ax2.semilogy(kk, [solnorm[k] / max(solnorm.values()) for k in kk],
                         "^--", label="||phi|| (norm.)")
            ax2.set_xlabel("modes kept")
            ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
            ax2.legend(); ax2.grid(True, which="both", alpha=0.3)
            out = HERE / "demo_coil_field_synthesis.png"
            fig.tight_layout(); fig.savefig(out, dpi=110)
            print(f"\nSaved plot: {out.name}")
        except ImportError:
            print("\n(matplotlib not installed; skipped plot)")
    finally:
        rad.UtiDelAll()


if __name__ == "__main__":
    main()
