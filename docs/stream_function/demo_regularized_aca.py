"""All 5 SF regularisations routed through the cached ``RegularizedTSVD``.

Demonstrates the regularisation-folded ACA+TSVD factorisation

    psi = (S^-1 V) . W^-1 . diag(1/Sigma) . U^T B,    W = V^T S^-1 V (k x k),

implemented in ``radia.stream_function.RegularizedTSVD``.  The factor
``S^-1 V`` is precomputed once; the (k x k) ``W`` is precomputed once;
each subsequent ``solve(B)`` is then ``O(k * (M + N))`` -- no further
sparse / dense solves.

This is the SAME factorisation re-used across:

    L2 (S = I)               -- equivalent to ``pseudo_inverse_solve``
    H1 (S = stiffness)       -- min surface current density L2 norm
    H1-sigma (1/sigma weighted) -- true ohmic dissipation
    Inductance-diag (h_cf weighted) -- self-inductance proxy
    L-inf via IRLS           -- iteratively rebuilds S^(k), reuses base

For Path-A iteration on the same regularisation, the ACA+TSVD base
factorisation runs ONCE and the same ``RegularizedTSVD`` cache is
reused across all iterations -- the per-iter cost drops from an O(N^2)
linear solve to a (k x N) + (k x k) + (k,) matvec chain (~ tens of us
at our N = 1773 / k = 25 scale).

Run:
    python demo_regularized_aca.py                   # default 5-mode sweep
    python demo_regularized_aca.py --order 3         # higher H1 order
    python demo_regularized_aca.py --skip-linf       # quicker, skip IRLS
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))

from radia.stream_function import (
    aca_tsvd, pseudo_inverse_solve, RegularizedTSVD,
)
from demo_planar_uniform_fem_psi import (
    build_fem_matrix, build_h1_stiffness, sample_psi_grid,
)
from demo_planar_uniform_coil import (
    contour_polylines_xy, single_stroke_spiral_xy, chain_xy_to_3d, bz_at,
)
from demo_planar_uniform_fem_psi_advanced import (
    solve_h1, solve_inductance_diagonal, solve_linf_irls,
)


def build_l2_stiffness(n_free):
    """Identity stiffness (S = I) -- gives L2 / Euclidean min-norm psi."""
    return np.eye(n_free)


def build_h1_sigma_stiffness(fes, free_idx, sigma_expr="1.0"):
    """1/sigma weighted H1 stiffness via NGSolve CF expression."""
    from ngsolve import (BilinearForm, grad, dx, TaskManager,
                          x as x_cf, y as y_cf, IfPos, exp, sqrt as ngs_sqrt,
                          sin, cos, log)
    u, v = fes.TnT()
    sigma_cf = eval(sigma_expr,
                    {"x": x_cf, "y": y_cf, "IfPos": IfPos,
                     "exp": exp, "sqrt": ngs_sqrt, "sin": sin, "cos": cos,
                     "log": log, "__builtins__": {}})
    a_stiff = BilinearForm(fes, symmetric=True)
    a_stiff += (1.0 / sigma_cf) * grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    return S[np.ix_(free_idx, free_idx)]


def build_inductance_diag_stiffness(fes, free_idx):
    """Element-size weighted H1 stiffness -- self-inductance proxy."""
    from ngsolve import (BilinearForm, grad, dx, TaskManager, specialcf)
    u, v = fes.TnT()
    h_cf = specialcf.mesh_size
    a_stiff = BilinearForm(fes, symmetric=True)
    a_stiff += h_cf * grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    return S[np.ix_(free_idx, free_idx)]


def field_metrics_from_psi(psi, A, fes, mesh, B_target, plane_half,
                            n_sample, nlevels, targets):
    """psi -> contours -> single-stroke chain -> field RMS metrics."""
    psi_grid, g_sample = sample_psi_grid(
        psi, fes, mesh, plane_half, n_sample)
    polylines, dI = contour_polylines_xy(psi_grid, g_sample, g_sample, nlevels)
    if not polylines:
        return {"cont_res": float("inf"), "rms": float("inf"),
                "p2p_mean": float("inf"), "n_wire": 0,
                "wire_length": 0.0, "I_w": 0.0}
    chain_xy = single_stroke_spiral_xy(polylines)
    path = chain_xy_to_3d(chain_xy)
    Bz_unit = bz_at(path, 1.0, targets)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    Bz_mean = float(np.mean(Bz_dsv))
    p2p = (float(np.max(Bz_dsv)) - float(np.min(Bz_dsv))) / (abs(Bz_mean) + 1e-30)
    cont_res = float(np.linalg.norm(A @ psi - B_target) /
                      (np.linalg.norm(B_target) + 1e-30))
    seg_lens = np.linalg.norm(np.diff(path, axis=0), axis=1)
    return {"cont_res": cont_res, "rms": rms, "p2p_mean": p2p,
            "n_wire": len(polylines), "wire_length": float(seg_lens.sum()),
            "I_w": I_w}


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
    ap.add_argument("--sigma-cf", type=str,
                    default="1.0 + 5.0 * exp(-(x*x + y*y) / 0.02)",
                    help="sigma(x, y) CF expression for h1_sigma")
    ap.add_argument("--jmax", type=float, default=4000.0,
                    help="L-inf cap (linf_irls only)")
    ap.add_argument("--skip-linf", action="store_true",
                    help="skip the IRLS mode (12 iters; slowest)")
    ap.add_argument("--aca-eps", type=float, default=1.0e-12,
                    help="ACA+ tolerance (tighter = more accurate "
                         "regularised solve)")
    args = ap.parse_args()

    print("=== Regularization-folded ACA+TSVD: 5-mode sweep ===")
    print(f"  geometry: square plane half = {args.plane_half} m,"
          f" target half = {args.target_half} m at z = {args.target_z} m")
    print(f"  FE: H1 order = {args.order}, max edge = {args.maxh} m")
    print(f"  target: B0 = {args.B0} T over {args.n_target}x{args.n_target}"
          f" grid")
    print()

    # ------- mesh + A matrix --------------------------------------------
    gt = np.linspace(-args.target_half, args.target_half, args.n_target)
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    targets = np.column_stack([Xt.ravel(), Yt.ravel(),
                                args.target_z * np.ones(Xt.size)])
    M = targets.shape[0]
    B_target = args.B0 * np.ones(M)

    t0 = time.time()
    A, fes, mesh, n_free = build_fem_matrix(
        args.plane_half, args.maxh, args.order, "boundary", targets)
    t_fem = time.time() - t0
    free_idx = np.where(np.array(fes.FreeDofs()))[0]
    A_free = A[:, free_idx]
    M_dim, N_dim = A_free.shape
    print(f"[Setup] FE matrix M = {M_dim} x N_free = {N_dim},"
          f" assembled in {t_fem:.2f} s")

    # ------- ACA+TSVD factorise A_free ONCE ------------------------------
    t0 = time.time()
    res = aca_tsvd(M_dim, N_dim, lambda i, j: float(A_free[i, j]),
                   modes=M_dim, kmax=min(M_dim, N_dim),
                   aca_eps=args.aca_eps)
    t_aca = time.time() - t0
    print(f"[ACA+TSVD] rank = {res.k_aca}/{res.modes},"
          f" {t_aca*1000:.1f} ms")
    print()

    # ------- per-mode sweep ---------------------------------------------
    modes = []
    modes.append(("l2",              "S = I (Euclidean min-norm)"))
    modes.append(("h1",              "S = stiffness (min |grad psi|^2)"))
    modes.append(("h1_sigma",        f"S = 1/sigma stiffness ({args.sigma_cf})"))
    modes.append(("inductance_diag", "S = h_cf * stiffness (inductance proxy)"))
    if not args.skip_linf:
        modes.append(("linf_irls", f"IRLS L-inf cap, jmax = {args.jmax}"))

    print(f"{'mode':<18}{'t_fold[ms]':<12}{'t_solve[ms]':<14}"
          f"{'cont_res':<12}{'RMS':<10}{'p2p/mean':<12}{'n_wire':<8}{'I_w[A]'}")
    print("-" * 96)

    results = []
    for name, desc in modes:
        # ---- build stiffness ----
        if name == "l2":
            S_free = build_l2_stiffness(N_dim)
        elif name == "h1":
            S_free = build_h1_stiffness(fes, free_idx)
        elif name == "h1_sigma":
            S_free = build_h1_sigma_stiffness(fes, free_idx, args.sigma_cf)
        elif name == "inductance_diag":
            S_free = build_inductance_diag_stiffness(fes, free_idx)
        elif name == "linf_irls":
            S_free = None     # IRLS rebuilds per-iter
        else:
            raise ValueError(f"unknown mode {name!r}")

        # ---- fold into cached ACA+TSVD form ----
        t0 = time.time()
        if name == "linf_irls":
            # IRLS does its own re-folding each iter; we report the FIRST iter cost
            psi, info = solve_linf_irls(A, B_target, fes, jmax=args.jmax,
                                         n_iter=12, verbose=False)
            t_fold = float("nan")    # not a single fold
        else:
            reg = RegularizedTSVD.from_stiffness(res, S_free)
        t1 = time.time()
        t_fold_ms = (t1 - t0) * 1000.0

        # ---- solve ----
        if name == "linf_irls":
            t_solve_ms = (t1 - t0) * 1000.0   # IRLS includes all 12 iters
            t_fold_ms = float("nan")
        else:
            t0 = time.time()
            psi_free = reg.solve(B_target)
            psi = np.zeros(A.shape[1])
            psi[free_idx] = psi_free
            t_solve_ms = (time.time() - t0) * 1000.0

        # ---- field metrics ----
        m = field_metrics_from_psi(
            psi, A, fes, mesh, B_target, args.plane_half,
            args.n_sample, args.nlevels, targets)

        fold_str = f"{t_fold_ms:>10.2f}" if not np.isnan(t_fold_ms) else "       (IRLS)"
        print(f"{name:<18}{fold_str:<12}{t_solve_ms:>11.2f}   "
              f"{m['cont_res']:<12.2e}{m['rms']*100:<8.2f}% "
              f"{m['p2p_mean']*100:<10.2f}% {m['n_wire']:<8}{m['I_w']:.2f}")

        results.append({"mode": name, "desc": desc,
                        "t_fold_ms": t_fold_ms, "t_solve_ms": t_solve_ms,
                        **m})

    print()
    print(f"NOTE: ACA+TSVD factorisation (rank {res.k_aca}) ran ONCE in"
          f" {t_aca*1000:.1f} ms; all modes above re-use it.")
    print(f"      Per-mode 'fold' cost is the one-time S^-1 V + W = V^T(S^-1 V)"
          f" + W_inv precomputation.")
    print(f"      Subsequent solve(B) calls (e.g. in Path-A iteration)"
          f" cost only O(k * N) at ~{results[1]['t_solve_ms']:.2f} ms.")

    # ------- JSON dump ---------------------------------------------------
    import json
    out = HERE / "demo_regularized_aca_results.json"
    payload = {
        "ndof": int(A.shape[1]),
        "n_free": int(N_dim),
        "M": int(M_dim),
        "order": args.order,
        "aca_rank": int(res.k_aca),
        "aca_eps": float(args.aca_eps),
        "t_aca_ms": t_aca * 1000.0,
        "results": results,
    }
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Saved sweep results: {out.name}")


if __name__ == "__main__":
    main()
