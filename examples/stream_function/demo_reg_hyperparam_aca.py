"""Optimise the REGULARISATION SHAPE (not the geometry) via Optuna, with
the ACA+TSVD factorisation reused across every trial.

This is the complement to ``demo_planar_uniform_fem_psi_advanced.py --deform``:

  - deform loop:  the SURFACE moves  ->  A changes every trial  ->  the
    whole RegularizedTSVD (ACA+TSVD base + S^-1 V + W) is rebuilt per trial.
  - THIS demo:    the surface is FIXED, only sigma(x, y) varies  ->  A is
    constant  ->  the ACA+TSVD base (U, Sigma, V) is computed ONCE and
    REUSED; each trial only rebuilds the cheap fold S^-1 V + W = V^T S^-1 V.

That is exactly the amortisation ``RegularizedTSVD`` was designed for.  At
the FE-direct scale (cheap matrix entries) the ACA+ base is only ~30 ms, so
the wall-clock win is modest; the POINT is the pattern: when each A(i, j) is
expensive (material kernel = Radia rad.Solve() + rad.Fld() per entry), the
ACA+ base dominates and reusing it across an Optuna sweep is the whole game.

What is optimised: a parameterised conductivity feature

    sigma(x, y) = 1 + amp * exp(-((x - cx)^2 + (y - cy)^2) / width)

The inner SF solve uses the 1/sigma-weighted H1 seminorm (true ohmic
dissipation form): min INT (1/sigma) |grad psi|^2 s.t. A psi = B.  The
sigma feature redistributes the surface current -- a dip (amp < 0) forbids
current near (cx, cy); a bump (amp > 0) attracts it.  Optuna searches the
(amp, cx, cy, width) that minimises the chosen cost.

Cost options (--cost):
  - ``chain_rms`` (default): single-stroke chain field RMS on the target
    plane.  The sigma redistribution changes the contour topology and
    hence the manufactured chain's field.
  - ``wire_length``: total single-stroke wire length (shorter = cheaper),
    with a hard penalty if chain_rms exceeds --rms-cap.
  - ``peak_grad``: peak |grad psi| (hot-spot suppression), penalised by
    chain_rms > --rms-cap.

Run:
    python demo_reg_hyperparam_aca.py --trials 30
    python demo_reg_hyperparam_aca.py --trials 30 --cost wire_length --rms-cap 0.03
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))

from radia.stream_function import aca_tsvd, RegularizedTSVD
from demo_planar_uniform_fem_psi import build_fem_matrix, sample_psi_grid
from demo_planar_uniform_coil import (
    contour_polylines_xy, single_stroke_spiral_xy, chain_xy_to_3d, bz_at,
)


def build_sigma_stiffness(fes, free_idx, amp, cx, cy, width):
    """1/sigma-weighted H1 stiffness with a Gaussian sigma feature.

    sigma(x, y) = 1 + amp * exp(-((x-cx)^2 + (y-cy)^2)/width)
    (amp > -1 keeps sigma > 0 -> 1/sigma finite -> S SPD).
    """
    from ngsolve import (BilinearForm, grad, dx, TaskManager,
                          x as x_cf, y as y_cf, exp as ngs_exp)
    u, v = fes.TnT()
    sigma_cf = 1.0 + amp * ngs_exp(
        -((x_cf - cx) ** 2 + (y_cf - cy) ** 2) / width)
    a_stiff = BilinearForm(fes, symmetric=True)
    a_stiff += (1.0 / sigma_cf) * grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    return S[np.ix_(free_idx, free_idx)]


def peak_gradient(psi, fes, free_idx):
    """Max |grad psi| sampled at mesh vertices."""
    from ngsolve import GridFunction, grad
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = psi
    grad_gf = grad(gf)
    mesh = fes.mesh
    peak = 0.0
    for vk in range(mesh.nv):
        mp = mesh(*mesh.vertices[vk].point)
        gv = grad_gf(mp)
        g = float(np.hypot(float(gv[0]), float(gv[1])))
        if g > peak:
            peak = g
    return peak


def evaluate_psi(psi, A, fes, mesh, B_target, plane_half, n_sample,
                 nlevels, targets):
    """psi -> contours -> single-stroke chain -> field metrics."""
    psi_grid, g_sample = sample_psi_grid(psi, fes, mesh, plane_half, n_sample)
    polylines, dI = contour_polylines_xy(psi_grid, g_sample, g_sample, nlevels)
    if not polylines:
        return {"chain_rms": float("inf"), "wire_length": float("inf"),
                "n_wire": 0, "I_w": 0.0}
    chain_xy = single_stroke_spiral_xy(polylines)
    path = chain_xy_to_3d(chain_xy)
    Bz_unit = bz_at(path, 1.0, targets)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    seg = np.linalg.norm(np.diff(path, axis=0), axis=1)
    return {"chain_rms": rms, "wire_length": float(seg.sum()),
            "n_wire": len(polylines), "I_w": I_w}


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
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--cost", choices=["chain_rms", "wire_length", "peak_grad"],
                    default="chain_rms")
    ap.add_argument("--rms-cap", type=float, default=0.03,
                    help="field RMS cap for wire_length / peak_grad cost "
                         "(penalty if exceeded)")
    ap.add_argument("--aca-eps", type=float, default=1.0e-12)
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    print("=== Reg-hyperparameter Optuna (fixed surface, ACA+ reused) ===")
    print(f"  cost = {args.cost}, {args.trials} trials,"
          f" sigma(x,y) = 1 + amp*exp(-r^2/width)")

    # ------- target + FIXED plane -> A built ONCE -----------------------
    gt = np.linspace(-args.target_half, args.target_half, args.n_target)
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    targets = np.column_stack([Xt.ravel(), Yt.ravel(),
                                args.target_z * np.ones(Xt.size)])
    M = targets.shape[0]
    B_target = args.B0 * np.ones(M)

    A, fes, mesh, n_free = build_fem_matrix(
        args.plane_half, args.maxh, args.order, "boundary", targets)
    free_idx = np.where(np.array(fes.FreeDofs()))[0]
    A_free = A[:, free_idx]
    M_dim, N_dim = A_free.shape

    # ------- ACA+TSVD of A: computed ONCE, reused every trial -----------
    t0 = time.time()
    res = aca_tsvd(M_dim, N_dim, lambda i, j: float(A_free[i, j]),
                   modes=M_dim, kmax=min(M_dim, N_dim),
                   aca_eps=args.aca_eps, method=3)
    t_aca = time.time() - t0
    print(f"  [ONCE] ACA+TSVD of A: rank {res.k_aca}/{res.modes},"
          f" {t_aca*1000:.1f} ms  <-- reused across all {args.trials} trials")

    timing = {"fold_ms": [], "eval_ms": []}

    def objective(trial):
        amp = trial.suggest_float("sigma_amp", -0.9, 5.0)
        cx = trial.suggest_float("sigma_cx", -0.18, 0.18)
        cy = trial.suggest_float("sigma_cy", -0.18, 0.18)
        width = trial.suggest_float("sigma_width", 0.003, 0.05, log=True)

        S_free = build_sigma_stiffness(fes, free_idx, amp, cx, cy, width)
        # *** ACA+ base `res` REUSED here -- only the fold is per-trial ***
        t1 = time.time()
        reg = RegularizedTSVD.from_stiffness(res, S_free)
        psi_free = reg.solve(B_target)
        timing["fold_ms"].append((time.time() - t1) * 1000.0)

        psi = np.zeros(A.shape[1])
        psi[free_idx] = psi_free

        t2 = time.time()
        m = evaluate_psi(psi, A, fes, mesh, B_target, args.plane_half,
                         args.n_sample, args.nlevels, targets)
        timing["eval_ms"].append((time.time() - t2) * 1000.0)

        trial.set_user_attr("chain_rms", m["chain_rms"])
        trial.set_user_attr("wire_length", m["wire_length"])
        trial.set_user_attr("n_wire", m["n_wire"])

        if args.cost == "chain_rms":
            return m["chain_rms"]
        # wire_length / peak_grad with rms cap penalty
        viol = max(0.0, m["chain_rms"] / args.rms_cap - 1.0)
        if args.cost == "wire_length":
            base = m["wire_length"]
        else:  # peak_grad
            base = peak_gradient(psi, fes, free_idx)
            trial.set_user_attr("peak_grad", base)
        return base * (1.0 + 100.0 * viol ** 2)

    sampler = optuna.samplers.CmaEsSampler(seed=42)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    t0 = time.time()
    study.optimize(objective, n_trials=args.trials)
    elapsed = time.time() - t0

    # ------- baseline: uniform sigma (= plain h1) for comparison --------
    S_base = build_sigma_stiffness(fes, free_idx, 0.0, 0.0, 0.0, 0.01)
    reg_base = RegularizedTSVD.from_stiffness(res, S_base)
    psi_base = np.zeros(A.shape[1])
    psi_base[free_idx] = reg_base.solve(B_target)
    m_base = evaluate_psi(psi_base, A, fes, mesh, B_target, args.plane_half,
                          args.n_sample, args.nlevels, targets)

    print(f"\n=== finished {args.trials} trials in {elapsed:.1f} s ===")
    print(f"  best cost      = {study.best_value:.4e}")
    print(f"  best sigma     = amp={study.best_params['sigma_amp']:+.3f},"
          f" c=({study.best_params['sigma_cx']:+.3f},"
          f" {study.best_params['sigma_cy']:+.3f}),"
          f" width={study.best_params['sigma_width']:.4f}")
    bt = study.best_trial
    print(f"  best chain RMS = {bt.user_attrs['chain_rms']*100:.2f}%"
          f"  (uniform-sigma baseline: {m_base['chain_rms']*100:.2f}%)")
    print(f"  best wire len  = {bt.user_attrs['wire_length']:.3f} m"
          f"  (baseline: {m_base['wire_length']:.3f} m)")

    # ------- timing breakdown: prove ACA+ amortisation ------------------
    fold = np.array(timing["fold_ms"])
    ev = np.array(timing["eval_ms"])
    print(f"\n  Per-trial cost breakdown ({len(fold)} trials):")
    print(f"    S build + fold (S^-1 V + W + solve): mean {fold.mean():.1f} ms")
    print(f"    chain extract + field eval:          mean {ev.mean():.1f} ms")
    print(f"  ACA+TSVD base: {t_aca*1000:.1f} ms ONCE"
          f" (would be {t_aca*1000*args.trials/1000:.1f} s if rebuilt per trial)")
    print(f"  At FE scale the entry is cheap so the saving is small; with a")
    print(f"  MATERIAL kernel (rad.Solve()+rad.Fld() per entry, ~100x slower)")
    print(f"  the ACA+ base would dominate and this reuse is the whole point.")

    # ------- JSON dump --------------------------------------------------
    import json
    out = HERE / "demo_reg_hyperparam_results.json"
    payload = {
        "cost": args.cost,
        "n_trials": args.trials,
        "order": args.order,
        "ndof": int(A.shape[1]),
        "n_free": int(N_dim),
        "M": int(M_dim),
        "aca_rank": int(res.k_aca),
        "t_aca_ms": t_aca * 1000.0,
        "t_aca_reused_across_trials": True,
        "mean_fold_ms": float(fold.mean()),
        "mean_eval_ms": float(ev.mean()),
        "best_cost": float(study.best_value),
        "best_params": dict(study.best_params),
        "best_chain_rms": float(bt.user_attrs["chain_rms"]),
        "baseline_uniform_sigma_chain_rms": float(m_base["chain_rms"]),
        "baseline_uniform_sigma_wire_length": float(m_base["wire_length"]),
        "trials": [
            {"chain_rms": t.user_attrs.get("chain_rms"),
             "wire_length": t.user_attrs.get("wire_length"),
             "n_wire": t.user_attrs.get("n_wire"),
             "params": dict(t.params)}
            for t in study.trials
            if t.user_attrs.get("chain_rms") is not None
        ],
    }
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\n  Saved: {out.name}")

    # ------- plot: sigma map + best psi contours + convergence ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42

        bp = study.best_params
        amp, cx, cy, width = (bp["sigma_amp"], bp["sigma_cx"],
                              bp["sigma_cy"], bp["sigma_width"])
        # best psi
        S_best = build_sigma_stiffness(fes, free_idx, amp, cx, cy, width)
        reg_best = RegularizedTSVD.from_stiffness(res, S_best)
        psi_best = np.zeros(A.shape[1])
        psi_best[free_idx] = reg_best.solve(B_target)
        psi_grid, g = sample_psi_grid(psi_best, fes, mesh, args.plane_half,
                                       args.n_sample)

        gg = np.linspace(-args.plane_half, args.plane_half, args.n_sample)
        GX, GY = np.meshgrid(gg, gg, indexing="ij")
        sigma_map = 1.0 + amp * np.exp(-((GX - cx) ** 2 + (GY - cy) ** 2) / width)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

        c0 = axes[0].contourf(GX * 1e3, GY * 1e3, sigma_map, levels=20,
                              cmap="viridis")
        axes[0].set_xlabel("$x$ (mm)"); axes[0].set_ylabel("$y$ (mm)")
        axes[0].set_aspect("equal")
        axes[0].text(0.02, 0.96, "(a)", transform=axes[0].transAxes, fontsize=10)
        fig.colorbar(c0, ax=axes[0], fraction=0.045)

        c1 = axes[1].contourf(g * 1e3, g * 1e3, psi_grid.T, levels=20,
                              cmap="RdBu_r")
        lo, hi = psi_grid.min(), psi_grid.max()
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            axes[1].contour(g * 1e3, g * 1e3, psi_grid.T, levels=lv,
                            colors="k", linewidths=0.5)
        axes[1].set_xlabel("$x$ (mm)"); axes[1].set_ylabel("$y$ (mm)")
        axes[1].set_aspect("equal")
        axes[1].text(0.02, 0.96, "(b)", transform=axes[1].transAxes, fontsize=10)
        fig.colorbar(c1, ax=axes[1], fraction=0.045)

        costs = [t.user_attrs.get("chain_rms") for t in study.trials
                 if t.user_attrs.get("chain_rms") is not None]
        running = np.minimum.accumulate(costs)
        axes[2].plot(np.array(costs) * 100, "o", ms=4, alpha=0.5,
                     label="trial chain RMS")
        axes[2].plot(np.array(running) * 100, "-", c="C3", lw=1.5,
                     label="running best")
        axes[2].axhline(m_base["chain_rms"] * 100, color="k", ls="--",
                        label=f"uniform-sigma {m_base['chain_rms']*100:.2f}%")
        axes[2].set_xlabel("trial"); axes[2].set_ylabel("chain RMS (%)")
        axes[2].text(0.02, 0.96, "(c)", transform=axes[2].transAxes, fontsize=10)
        axes[2].legend(loc="best", fontsize=8); axes[2].grid(alpha=0.3)

        out_png = HERE / "demo_reg_hyperparam_plot.png"
        fig.tight_layout()
        fig.savefig(out_png, dpi=110)
        plt.close(fig)
        print(f"  Saved: {out_png.name}")
    except ImportError:
        print("  (matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
