"""Same (ACA+)+TSVD solver on a permanent-magnet array (magnetic material).

This is the SAME solver as `demo_coil_field_synthesis.py`, but the basis
sources are permanent-magnet cubes (`radia.ObjRecMag`) instead of current
loops.  The matrix entry now comes from Radia's MMM / MSC field (surface-charge
analytical formulas) -- again via `radia_field_kernel` + `radia.Fld`, with zero
changes to the solver.  This demonstrates that the (ACA+)+TSVD machinery is
kernel-agnostic: it serves coils and magnetic materials alike.

Run:  python demo_magnet_array.py
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import radia as rad
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel


def main():
    rad.UtiDelAll()
    try:
        # --- N permanent-magnet cubes on a plane -----------------------------
        n_side = 7                                   # N = 49 magnets
        xs = np.linspace(-0.06, 0.06, n_side)
        ys = np.linspace(-0.06, 0.06, n_side)
        sources = []
        for cx in xs:
            for cy in ys:
                # 10 mm cube, Br ~ 1.2 T  (M = Br/mu0 ~ 954930 A/m), +z.
                sources.append(rad.magnet_box([float(cx), float(cy), 0.0],
                                             [0.01, 0.01, 0.01],
                                             [0.0, 0.0, 954930.0]))
        N = len(sources)

        # --- M observation points above the array ---------------------------
        m_side = 4                                   # M = 16 field points
        oxs = np.linspace(-0.04, 0.04, m_side)
        oys = np.linspace(-0.04, 0.04, m_side)
        ox, oy = np.meshgrid(oxs, oys, indexing="ij")
        obs = np.column_stack([ox.ravel(), oy.ravel(),
                               np.full(ox.size, 0.06)])
        M = obs.shape[0]

        # --- matrix entry from Radia's MMM/MSC field -------------------------
        entry = radia_field_kernel(obs, sources, component=2, field="b")
        A = np.array([[entry(i, j) for j in range(N)] for i in range(M)])

        # --- (ACA+)+TSVD factorization ---------------------------------------
        res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                       aca_eps=1.0e-8, method=3)
        Alr = res.U @ np.diag(res.S) @ res.V.T
        recon = np.linalg.norm(A - Alr) / np.linalg.norm(A)

        print("=== Permanent-magnet array: (ACA+)+TSVD (magnetic material) ===")
        print(f"M (field points) = {M},  N (magnets) = {N}")
        print(f"ACA+ rank k_aca  = {res.k_aca}  (of min(M,N) = {min(M, N)})")
        print(f"singular values S[:6] = "
              f"{np.array2string(res.S[:6], precision=4)}")
        print(f"reconstruction ||A - U S V^T|| / ||A|| = {recon:.3e}")

        # --- least-norm solve: recover a target in range(A) ------------------
        rng = np.random.default_rng(0)
        phi0 = rng.standard_normal(N)
        B = A @ phi0
        phi = pseudo_inverse_solve(res, B, k_mode=res.k_aca)
        solve_err = np.linalg.norm(A @ phi - B) / np.linalg.norm(B)
        print(f"least-norm solve ||A phi - B|| / ||B|| = {solve_err:.3e} "
              f"(||phi||={np.linalg.norm(phi):.3e}, ||phi0||={np.linalg.norm(phi0):.3e})")
        print("\nSame solver, magnetic-material sources -- no coil-specific code.")

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams["pdf.fonttype"] = 42
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.semilogy(np.arange(1, res.k_aca + 1), res.S[:res.k_aca], "o-")
            ax.set_xlabel("mode index"); ax.set_ylabel("singular value")
            ax.grid(True, which="both", alpha=0.3)
            out = HERE / "demo_magnet_array.png"
            fig.tight_layout(); fig.savefig(out, dpi=110)
            print(f"Saved plot: {out.name}")
        except ImportError:
            print("(matplotlib not installed; skipped plot)")
    finally:
        rad.UtiDelAll()


if __name__ == "__main__":
    main()
