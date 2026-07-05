"""SF coil design: H1 FE psi + Radia (ACA+)+TSVD via HACApK.

Same target as ``demo_planar_uniform_fem_psi.py`` (uniform Bz=B0 at z=h
over a square target above a square planar source), and same H1 FE
formulation, but the matrix-vector machinery now goes through
``radia.stream_function.aca_tsvd`` -- HACApK's ACA+ + Method-3 TSVD.

Why this matters for the SFM pipeline:

  - **ACA+ is kernel-agnostic** (CLAUDE.md "FMM Removed (2026-03-06)" /
    "SCOPE CLARIFICATION 2026-05-30").  Validating that it works on a
    FE-direct matrix opens the door to MATERIAL KERNELS (Radia
    HDiv-VIM iron yoke, shielded coil, SIBC workpiece, ...) where each
    entry costs a Radia container solve and ACA+ is essential.
  - The Path-A compensated iteration RE-USES the cached TSVD
    factorisation across iters (one ACA+ build + N back-substitutions).
    Confirms the speed argument carries over from basis-loop to FE-
    direct without changing the outer loop.
  - When ngsolve.bem ships its native H-matrix (Joachim's intent
    confirmed 2026-05-30 「JOACHIMは、H-matrixを実装したがっていた」),
    the only change in this file will be the matrix-assembly section:
    replace ``build_fem_matrix`` (full LinearForm per target) with an
    H-matrix-backed ``ngsolve.bem`` operator.  The ``entry`` callback
    becomes an operator-application; everything below stays.

Pipeline:
  [0] Same H1 mesh / FES as the LinearForm demo
  [1] Build M x N_free matrix (interim -- replace with ngsbem hmatrix)
  [2] Wrap as ``entry(i, j) -> float`` for radia.stream_function.aca_tsvd
  [3] ACA+ + Method-3 TSVD; cache the factorisation
  [4] Pseudo-inverse solve: psi_free = V diag(1/S) U^T B_target
  [5] Sample psi on grid, extract contours, single-stroke spiral chain
  [6] Optional Path-A iteration -- delta_psi = pseudo_inverse_solve(res, residual)
  [7] Field verification

Run: python demo_planar_uniform_fem_psi_aca.py [--order 2] [--maxh 0.025]
                                                 [--aca-eps 1e-8]
                                                 [--compensated-iter N]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.stream_function import aca_tsvd, pseudo_inverse_solve
from radia.biot_savart import h_segments_batch, MU0

from demo_planar_uniform_coil import (
    contour_polylines_xy, single_stroke_spiral_xy, chain_xy_to_3d, bz_at,
)
from demo_planar_uniform_fem_psi import (
    build_fem_matrix, sample_psi_grid,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plane-half", type=float, default=0.25)
    ap.add_argument("--maxh", type=float, default=0.025)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--target-half", type=float, default=0.05)
    ap.add_argument("--target-z", type=float, default=0.10)
    ap.add_argument("--n-target", type=int, default=5)
    ap.add_argument("--B0", type=float, default=1.0e-3)
    ap.add_argument("--n-sample", type=int, default=81)
    ap.add_argument("--nlevels", type=int, default=12)
    ap.add_argument("--aca-eps", type=float, default=1.0e-10,
                    help="ACA+ stopping tolerance.  Below the singular-value "
                         "floor of A_free for a clean reconstruction.")
    ap.add_argument("--aca-modes", type=int, default=0,
                    help="TSVD modes (0 = auto = min(M, k_aca)).")
    ap.add_argument("--compensated-iter", type=int, default=0)
    ap.add_argument("--compensated-step", type=float, default=1.0)
    args = ap.parse_args()

    print("=== SF -> planar uniform Bz (H1 FE psi via ACA+TSVD on FE matrix) ===")

    # --- targets ---
    gt = np.linspace(-args.target_half, args.target_half, args.n_target)
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    targets = np.column_stack([Xt.ravel(), Yt.ravel(),
                                args.target_z * np.ones(Xt.size)])
    M = targets.shape[0]
    B_target = args.B0 * np.ones(M)
    print(f"[0] target: {M} obs at z={args.target_z*1e3:.0f} mm, "
          f"square half={args.target_half*1e3:.0f} mm, B0={args.B0*1e3:.2f} mT")

    # --- FE matrix (interim path; replace with ngsolve.bem operator later) ---
    print("[1] FE matrix assembly via LinearForm (interim path)")
    t0 = time.time()
    A, fes, mesh, n_free = build_fem_matrix(
        args.plane_half, args.maxh, args.order, "boundary", targets)
    free_idx = np.where(np.array(fes.FreeDofs()))[0]
    A_free = A[:, free_idx]
    t_assembly = time.time() - t0
    print(f"    A_free: M={M} x N_free={len(free_idx)}, assembled in {t_assembly:.2f} s")

    # --- (ACA+)+TSVD on FE matrix via callback ---
    print("[2] (ACA+)+TSVD via radia.stream_function.aca_tsvd")
    def entry(i, j):
        return float(A_free[i, j])
    modes = args.aca_modes if args.aca_modes > 0 else M
    t0 = time.time()
    res = aca_tsvd(M, len(free_idx), entry, modes=modes,
                   kmax=min(M, len(free_idx)),
                   aca_eps=args.aca_eps)
    t_aca = time.time() - t0
    print(f"    k_aca = {res.k_aca} (of min(M,N)={min(M, len(free_idx))}),"
          f" TSVD modes used = {res.modes}, factorised in {t_aca:.3f} s")

    # --- pseudo-inverse solve ---
    psi_free = pseudo_inverse_solve(res, B_target, k_mode=res.modes)
    psi_vec = np.zeros(A.shape[1])
    psi_vec[free_idx] = psi_free
    A_psi = A_free @ psi_free
    cont_res = float(np.linalg.norm(A_psi - B_target) /
                     (np.linalg.norm(B_target) + 1e-30))
    print(f"    continuous SF: ||A psi - B0||/||B0|| = {cont_res:.3e}")

    # --- sanity check vs direct lstsq baseline ---
    psi_lstsq, _, rank, _ = np.linalg.lstsq(A_free, B_target, rcond=None)
    psi_norm_diff = float(np.linalg.norm(psi_free - psi_lstsq) /
                          (np.linalg.norm(psi_lstsq) + 1e-30))
    print(f"    sanity: ||psi_ACA - psi_lstsq||/||psi_lstsq|| = {psi_norm_diff:.3e}"
          f" (lstsq rank = {rank})")

    # --- sample on grid ---
    psi_grid, g_sample = sample_psi_grid(
        psi_vec, fes, mesh, args.plane_half, args.n_sample)
    polylines, dI = contour_polylines_xy(
        psi_grid, g_sample, g_sample, args.nlevels)
    n_wire = len(polylines)
    if n_wire == 0:
        print("[3] no usable contours; aborting")
        return
    print(f"[3] contours: {n_wire} closed curves, dI = {dI:.4g}")

    # --- Path-A iteration with cached ACA+TSVD factorisation ---
    if args.compensated_iter > 0:
        print(f"[3a] Path-A compensated iteration "
              f"({args.compensated_iter} iters, step={args.compensated_step})")
        best_psi = psi_vec.copy()
        best_polylines = polylines
        best_res_norm = float("inf")
        for it in range(args.compensated_iter):
            chain_it = single_stroke_spiral_xy(polylines)
            path_it = chain_xy_to_3d(chain_it)
            Bz_unit_it = bz_at(path_it, 1.0, targets)
            denom_it = float(np.dot(Bz_unit_it, Bz_unit_it))
            I_w_it = float(np.dot(Bz_unit_it, B_target) / denom_it) \
                if denom_it > 0 else 0.0
            residual = B_target - I_w_it * Bz_unit_it
            res_norm = float(np.linalg.norm(residual) /
                             (np.linalg.norm(B_target) + 1e-30))
            tag = ""
            if res_norm < best_res_norm:
                best_res_norm = res_norm
                best_psi = psi_vec.copy()
                best_polylines = polylines
                tag = " <-- best"
            print(f"     iter {it+1}: I_w = {I_w_it:.3e},"
                  f" residual = {res_norm:.3e}{tag}")
            # cached factorisation re-used here -- no new ACA+ build
            delta_psi_free = pseudo_inverse_solve(res, residual,
                                                   k_mode=res.modes)
            psi_vec[free_idx] += args.compensated_step * delta_psi_free
            psi_grid, g_sample = sample_psi_grid(
                psi_vec, fes, mesh, args.plane_half, args.n_sample)
            polylines, dI = contour_polylines_xy(
                psi_grid, g_sample, g_sample, args.nlevels)
            if not polylines:
                print("     iter aborted: psi became flat")
                break
        psi_vec = best_psi
        polylines = best_polylines
        psi_grid, g_sample = sample_psi_grid(
            psi_vec, fes, mesh, args.plane_half, args.n_sample)
        print(f"     final: best residual = {best_res_norm:.3e}")

    # --- single-stroke spiral chain ---
    chain_xy = single_stroke_spiral_xy(polylines)
    path = chain_xy_to_3d(chain_xy)
    seg_lens = np.linalg.norm(np.diff(path, axis=0), axis=1)
    length = float(seg_lens.sum())
    print(f"[4] spiral chain: {len(path)} pts, {len(path)-1} segs,"
          f" length = {length:.3f} m")

    # --- field eval ---
    Bz_unit = bz_at(path, 1.0, targets)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    Bz_mean = float(np.mean(Bz_dsv))
    p2p = (float(np.max(Bz_dsv)) - float(np.min(Bz_dsv))) / (abs(Bz_mean) + 1e-30)
    print(f"[6] best-fit I_w = {I_w:.4g} A,"
          f" target-plane RMS = {rms:.3e}")
    print(f"    Bz over target: mean = {Bz_mean:.3e} T,"
          f" peak-to-peak / mean = {p2p:.3e}")
    print(f"    (target B0 = {args.B0:.3e} T)")

    zv = np.linspace(0.5 * args.target_z, 1.5 * args.target_z, 41)
    obs_axis = np.column_stack([np.zeros_like(zv), np.zeros_like(zv), zv])
    Bz_axis = bz_at(path, I_w, obs_axis)
    print(f"[7] on-axis Bz at z={args.target_z*1e3:.0f}mm: {Bz_axis[20]:.3e} T")

    # --- plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig = plt.figure(figsize=(13, 4))

        ax1 = fig.add_subplot(1, 3, 1)
        cf = ax1.contourf(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                          levels=20, cmap="RdBu_r")
        lo, hi = psi_grid.min(), psi_grid.max()
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            ax1.contour(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                        levels=lv, colors="k", linewidths=0.5)
        ax1.set_xlabel("$x$ (mm)"); ax1.set_ylabel("$y$ (mm)")
        ax1.set_aspect("equal")
        ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=10)
        fig.colorbar(cf, ax=ax1, fraction=0.045)

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, "b-", lw=0.5)
        ax2.set_xlabel("$x$ (mm)"); ax2.set_ylabel("$y$ (mm)")
        ax2.set_aspect("equal")
        ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
        ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.plot(zv * 1e3, Bz_axis * 1e3, "o-", ms=3, label="single-stroke coil")
        ax3.axhline(args.B0 * 1e3, color="k", ls="--",
                    label=f"target {args.B0*1e3:.2f} mT")
        ax3.axvline(args.target_z * 1e3, color="g", ls=":",
                    alpha=0.5, label=f"target z={args.target_z*1e3:.0f}mm")
        ax3.set_xlabel("$z$ (mm)"); ax3.set_ylabel("$B_z$ (mT)")
        ax3.text(0.02, 0.96, "(c)", transform=ax3.transAxes, fontsize=10)
        ax3.legend(loc="best", fontsize=8); ax3.grid(alpha=0.3)

        out = HERE / "demo_planar_uniform_fem_psi_aca.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
