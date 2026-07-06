r"""Excitation-INVARIANT flux lines: keep the SAME field-line pattern as the drive
current rises (NOT a cyclotron, where the field is meant to change with radius).

THE QUESTION (the user's "same flux lines even when the current is increased")
------------------------------------------------------------------------------
A current-driven electromagnet is a LINEAR magnetostatic system as long as the iron
stays below its knee: scaling the excitation by alpha scales B everywhere by alpha, so
the flux-LINE pattern (the streamlines b_hat = B/|B|) is IDENTICAL -- only the amplitude
grows.  "Same flux lines when you turn up the current" is therefore AUTOMATIC in the
linear regime; the ONLY thing that can change the pattern is the nonlinearity, i.e. iron
SATURATION (mu(|B|) dropping non-uniformly, first at the pole-tip corner), which
redistributes the flux and rotates the field lines.

So the design target is: make the flux-line pattern stay invariant DEEP into saturation.
The measurable quantity is the flux-line DIRECTION drift, as a function of excitation,

    D_dir(I) = rms_x || b_hat(x; I) - b_hat(x; I_lin) ||   over an AIR sampling region,

where I_lin is a low (linear-regime) drive.  D_dir = 0 while linear; it grows once the
iron saturates.  (The MEDIAN-plane longitudinal profile p(s) = B_z(s)/B_body does NOT
show this -- it is gap-reluctance-robust and stays ~1e-3 for any end shape; see
bending_endpack_saturation_opt.py.  The AIR-region field-line DIRECTION near the pole
tip / fringe is where saturation actually moves the pattern, so that is what we measure.)

WHAT IS COMPUTED
----------------
1. air_grid(): a fixed set of (s,z) points in the AIR gap + fringe near the pole end.
2. solve_bhat(geom, B_drive, saturate): the nonlinear scalar-potential solve (Froehlich
   mu(|B|) Picard) -> the unit field-line direction b_hat = B/|B| at those air points.
   saturate=False forces mu = MUR0 (the linear reference: b_hat is drive-INDEPENDENT).
3. flux_line_drift(depth, exponent, B_lin, B_sat): D_dir at the saturated drive vs the
   low (linear) drive -- 2 solves; the OPTIMIZATION OBJECTIVE (D_dir is monotone in the
   drive, so the top-drive value is the worst case).
4. invariance_curve(depth, exponent, drives): the full D_dir(B_drive) SATURATED sweep,
   plus the LINEAR control (saturate=False -> D_dir ~ 0 for every drive: the proof that
   linearity => invariance and that saturation is the sole breaker).
5. optimize(): minimize the saturated flux-line direction drift over the end chamfer
   (depth, exponent) -- the shape that keeps the flux lines invariant deepest into
   saturation -- vs the flat-cut baseline.

The honest result: the flux lines are ALREADY remarkably excitation-invariant for this
magnet class (sub-degree direction drift even for a hard flat cut, because the high-mu
pole face stays equipotential -- gap-reluctance robustness).  The residual drift is
corner-dominated, and relieving the pole-tip corner with an end chamfer makes the flux
lines several times MORE invariant.  So "same flux lines as the current rises" is largely
automatic; what you optimize is the corner that would otherwise break it under saturation.

run:  python excitation_invariant_field.py           # optimize + flat-cut compare
      python excitation_invariant_field.py --fig        # + figure
      python excitation_invariant_field.py --trials 40  # Optuna budget
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Reuse the bending end-pack geometry + Froehlich constants (same (s,z) small-gap pole
# terminating at +-s_pole with a parametrized end chamfer).  This file adds the
# EXCITATION-SWEEP flux-line-direction invariance metric on top of that geometry.
import bending_endpack_saturation_opt as be  # noqa: E402


def air_grid():
    """Fixed (s,z) points GUARANTEED in air for any (depth, exponent) in the search
    range: under the gap face (z <= 0.78*G2 stays below the >=G2 pole face) plus a
    fringe patch just past the pole end (where the flux lines bend the most, so where
    saturation rotates them the most)."""
    pts = []
    for s in np.linspace(0.0, 1.32 * be.S_POLE, 34):
        for z in np.linspace(0.30 * be.G2, 0.78 * be.G2, 4):
            pts.append((s, z))
    for s in np.linspace(1.02 * be.S_POLE, 1.34 * be.S_POLE, 8):
        for z in np.linspace(0.30 * be.G2, 1.5 * be.G2, 5):
            pts.append((s, z))
    return np.array(pts)


def solve_bhat(mesh, B_drive, pts, order=3, saturate=True, relax=0.5, tol=1e-4, maxit=60):
    """Nonlinear scalar-potential solve (Froehlich mu(|B|) Picard).  Return the unit
    field-line direction b_hat = B/|B| = -grad(psi)/|grad(psi)| (mu_r=1 in air) at the
    air sampling points.  saturate=False forces mu = MUR0 (the linear reference, whose
    b_hat is DRIVE-INDEPENDENT by linearity)."""
    from ngsolve import (H1, L2, BilinearForm, GridFunction, grad, dx, CF, Norm,
                         Integrate, TaskManager)
    mmf = B_drive * be.GAP / (2.0 * be.MU0)
    with TaskManager():
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        gfu = GridFunction(fes)
        bccf = mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0)
        fes_mu = L2(mesh, order=0)
        mu_gf = GridFunction(fes_mu)
        mu_gf.Set(mesh.MaterialCF({"iron": be.MUR0}, default=1.0))
        iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
        it, resid = 0, 1.0
        for it in range(1, maxit + 1):
            a = BilinearForm(fes)
            a += mu_gf * grad(u) * grad(v) * dx
            a.Assemble()
            gfu.Set(bccf, definedon=mesh.Boundaries("median|irontop"))
            r = gfu.vec.CreateVector()
            r.data = -a.mat * gfu.vec
            gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
            if not saturate:
                break                                     # linear reference (mu = MUR0)
            B = be.MU0 * mu_gf * Norm(grad(gfu))
            froh = 1.0 + (be.MUR0 - 1.0) / (1.0 + (B / be.BK) ** 2)
            mu_t = (1.0 - iron_ind) * CF(1.0) + iron_ind * froh
            mu_n = GridFunction(fes_mu)
            mu_n.Set(mu_t)
            d = mu_n.vec.CreateVector()
            d.data = mu_n.vec - mu_gf.vec
            resid = d.Norm() / (mu_gf.vec.Norm() or 1.0)
            mu_gf.vec.data += relax * d
            if resid < tol:
                break
        mur_mean = 0.0
        iv = float(Integrate(iron_ind, mesh))
        if iv > 0:
            mur_mean = float(Integrate(iron_ind * mu_gf, mesh)) / iv
        g = grad(gfu)
        bh = np.zeros((len(pts), 2))
        for i, (s, z) in enumerate(pts):
            gg = g(mesh(float(s), float(z)))
            b = np.array([-gg[0], -gg[1]])                 # B ~ -grad(psi) in air
            n = np.hypot(b[0], b[1])
            bh[i] = b / n if n > 1e-30 else 0.0
    return bh, float(mur_mean), int(it)


def _drift(bh, ref):
    """rms flux-line direction difference (radians, small-angle ~ chord length)."""
    return float(np.sqrt(np.mean(np.sum((bh - ref) ** 2, axis=1))))


def flux_line_drift(depth, exponent, B_lin=0.15, B_sat=1.70, order=3, maxh=None,
                    pts=None, mesh=None):
    """The OBJECTIVE: the flux-line DIRECTION drift between a low (linear-regime) drive
    and a high (saturated) drive.  2 solves.  Small = the flux lines keep their shape as
    the current is turned up.  Also returns the saturated iron <mu_r> (how deep it went)
    and the pole-tip corner kappa (context: the corner is what drives the drift)."""
    if mesh is None:
        mesh = be.build_endpack(depth=depth, exponent=exponent, maxh=maxh)
    if pts is None:
        pts = air_grid()
    bh_ref, mur_ref, _ = solve_bhat(mesh, B_lin, pts, order=order, saturate=True)
    bh_sat, mur_sat, _ = solve_bhat(mesh, B_sat, pts, order=order, saturate=True)
    prof = be.solve_profile(mesh, B_sat, order=order, saturate=True)  # corner context
    return {
        "depth_m": float(depth), "exponent": float(exponent),
        "D_dir": _drift(bh_sat, bh_ref),                   # the objective
        "mur_ref": mur_ref, "mur_sat": mur_sat,
        "corner_kappa_sat": float(prof["corner_kappa"]),
        "_bh_ref": bh_ref, "_bh_sat": bh_sat, "_pts": pts,
    }


def invariance_curve(depth, exponent, drives=(0.15, 0.45, 0.75, 1.05, 1.35, 1.70),
                     order=3, maxh=None):
    """The full D_dir(B_drive) SATURATED sweep (reporting / figure), plus the LINEAR
    control (mu forced constant -> D_dir ~ 0 at every drive: linearity => invariant flux
    lines; saturation is the sole breaker).  Reference = the lowest drive."""
    mesh = be.build_endpack(depth=depth, exponent=exponent, maxh=maxh)
    pts = air_grid()
    drives = list(drives)
    bh_ref_sat, _, _ = solve_bhat(mesh, drives[0], pts, order=order, saturate=True)
    bh_ref_lin, _, _ = solve_bhat(mesh, drives[0], pts, order=order, saturate=False)
    D_sat, D_lin, murs = [], [], []
    for B in drives:
        bh_s, mur_s, _ = solve_bhat(mesh, B, pts, order=order, saturate=True)
        bh_l, _, _ = solve_bhat(mesh, B, pts, order=order, saturate=False)
        D_sat.append(_drift(bh_s, bh_ref_sat))
        D_lin.append(_drift(bh_l, bh_ref_lin))
        murs.append(mur_s)
    return {
        "depth_m": float(depth), "exponent": float(exponent),
        "drives_T": drives, "D_dir_saturated": D_sat, "D_dir_linear_control": D_lin,
        "mur_sat_sweep": murs, "D_dir_max": float(max(D_sat)),
        "linear_control_max": float(max(D_lin)),
    }


def optimize(trials=30, B_lin=0.15, B_sat=1.70, order=3, seed=0):
    """Minimize the saturated flux-line DIRECTION drift D_dir(B_lin -> B_sat) over the
    end chamfer (depth, exponent) -- the shape that keeps the flux lines invariant
    deepest into saturation -- vs the flat-cut baseline.  Optuna (TPE) if available,
    else a coarse grid + local refine.  D_dir is monotone in the drive, so the
    top-drive drift is the worst-case invariance error."""
    depth_lo, depth_hi = 0.0, 0.95 * be.G2
    exp_lo, exp_hi = 0.4, 3.0
    pts = air_grid()
    flat = flux_line_drift(0.0, 1.0, B_lin=B_lin, B_sat=B_sat, order=order, pts=pts)
    history = []

    def obj(depth, exponent):
        r = flux_line_drift(depth, exponent, B_lin=B_lin, B_sat=B_sat, order=order,
                            pts=pts)
        history.append({"depth_m": r["depth_m"], "exponent": r["exponent"],
                        "D_dir": r["D_dir"], "corner_kappa_sat": r["corner_kappa_sat"]})
        return r

    def cost(r):
        return r["D_dir"]

    best, used = None, "grid"
    try:
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        used = "optuna_tpe"

        def _obj(trial):
            d = trial.suggest_float("depth", depth_lo, depth_hi)
            e = trial.suggest_float("exponent", exp_lo, exp_hi)
            return cost(obj(d, e))

        study = optuna.create_study(direction="minimize",
                                    sampler=optuna.samplers.TPESampler(seed=seed))
        study.optimize(_obj, n_trials=trials, show_progress_bar=False)
        bp = study.best_params
        best = flux_line_drift(bp["depth"], bp["exponent"], B_lin=B_lin, B_sat=B_sat,
                               order=order, pts=pts)
    except Exception:
        used = "grid_refine"
        depths = np.linspace(depth_lo, depth_hi, 5)
        exps = np.linspace(exp_lo, exp_hi, 5)
        grid = [obj(float(d), float(e)) for d in depths for e in exps]
        best = min(grid, key=cost)
        d0, e0 = best["depth_m"], best["exponent"]
        for dd in (d0 - 0.15 * be.G2, d0, d0 + 0.15 * be.G2):
            for ee in (e0 - 0.4, e0, e0 + 0.4):
                if depth_lo <= dd <= depth_hi and exp_lo <= ee <= exp_hi:
                    r = obj(float(dd), float(ee))
                    if cost(r) < cost(best):
                        best = r

    inv_factor = (flat["D_dir"] / best["D_dir"]) if best["D_dir"] > 0 else float("inf")
    _REPORT = ("depth_m", "exponent", "D_dir", "mur_sat", "corner_kappa_sat")
    return {
        "optimizer": used, "n_evals": len(history),
        "B_lin_T": float(B_lin), "B_sat_T": float(B_sat), "B_K_iron_T": float(be.BK),
        "flat_cut": {k: flat[k] for k in _REPORT},
        "optimized": {k: best[k] for k in _REPORT},
        "invariance_factor": float(inv_factor),
        "flux_lines_more_invariant": bool(best["D_dir"] < flat["D_dir"]),
        "flat_flux_lines_already_invariant": bool(flat["D_dir"] < 1e-2),  # sub ~0.6 deg
        "history": history, "_flat": flat, "_best": best,
    }


def _figure(res, flat_curve, best_curve):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.3), dpi=140)
    fc, oc = res["flat_cut"], res["optimized"]

    # LEFT: flux-line direction drift vs EXCITATION -- flat vs optimized, + linear control
    dr = np.array(flat_curve["drives_T"])
    ax[0].plot(dr, np.array(flat_curve["D_dir_saturated"]) * 1e3, "C3-o", lw=1.6, ms=4,
               label="flat cut (saturated)")
    ax[0].plot(dr, np.array(best_curve["D_dir_saturated"]) * 1e3, "C0-o", lw=1.6, ms=4,
               label="optimized (saturated)")
    ax[0].plot(dr, np.array(flat_curve["D_dir_linear_control"]) * 1e3, "k:", lw=1.2,
               label="linear control (mu=const)")
    ax[0].axvline(be.BK, color="0.6", lw=0.8, ls="--")
    ax[0].text(be.BK, ax[0].get_ylim()[1] * 0.92, " iron knee", fontsize=7.5, color="0.4")
    ax[0].set_xlabel("drive  $B_{body}$ [T]  (excitation)")
    ax[0].set_ylabel("flux-line direction drift  $D_{dir}$  [mrad]")
    ax[0].set_title("Flux lines stay invariant while linear,\n"
                    "drift once the iron saturates (control ~ 0)")
    ax[0].legend(fontsize=7.5)

    # MIDDLE: the invariance at the top drive -- flat vs optimized (the lever)
    x = np.arange(2)
    dd = [fc["D_dir"] * 1e3, oc["D_dir"] * 1e3]
    bars = ax[1].bar(x, dd, color=["C3", "C2"])
    for b_, d_ in zip(bars, dd):
        ax[1].text(b_.get_x() + b_.get_width() / 2, d_, f"{d_:.2f}\nmrad",
                   ha="center", va="bottom", fontsize=9)
    ax[1].set_xticks(x); ax[1].set_xticklabels(["FLAT cut", "OPTIMIZED"])
    ax[1].set_ylabel(f"$D_{{dir}}$ at $B$={res['B_sat_T']:.1f} T  [mrad]")
    ax[1].set_ylim(0, max(dd) * 1.28)
    ax[1].set_title(f"End chamfer keeps flux lines invariant deeper\n"
                    f"({res['invariance_factor']:.1f}x smaller drift; "
                    f"d={oc['depth_m']*1e3:.1f} mm, p={oc['exponent']:.2f})")

    # RIGHT: the optimization -- D_dir(depth) colored by exponent
    h = res["history"]
    d = np.array([x["depth_m"] for x in h]) * 1e3
    e = np.array([x["exponent"] for x in h])
    yy = np.array([x["D_dir"] for x in h]) * 1e3
    sc = ax[2].scatter(d, yy, c=e, cmap="viridis", s=28)
    ax[2].axhline(fc["D_dir"] * 1e3, color="C3", lw=1, ls="--",
                  label=f"flat cut {fc['D_dir']*1e3:.2f} mrad")
    ax[2].scatter([oc["depth_m"] * 1e3], [oc["D_dir"] * 1e3], marker="*", s=220,
                  color="C1", edgecolor="k", zorder=5, label="optimum")
    ax[2].set_xlabel("chamfer depth [mm]")
    ax[2].set_ylabel("$D_{dir}$ at top drive [mrad] (lower = more invariant)")
    ax[2].set_title(f"Optimization ({res['optimizer']}, {res['n_evals']} evals):\n"
                    "the flux-line-invariance-maximizing end chamfer")
    ax[2].legend(fontsize=8)
    fig.colorbar(sc, ax=ax[2], label="chamfer exponent")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--trials", type=int, default=30)
    args = ap.parse_args()

    print("=" * 78)
    print("EXCITATION-INVARIANT flux lines: same field-line shape as the current rises")
    print("=" * 78)

    res = optimize(trials=args.trials)
    f, o = res["flat_cut"], res["optimized"]
    print(f"\niron knee B_K = {res['B_K_iron_T']:.2f} T; gap = {be.GAP*1e3:.0f} mm; "
          f"linear drive {res['B_lin_T']:.2f} T, saturated drive {res['B_sat_T']:.2f} T")
    print(f"optimizer = {res['optimizer']} ({res['n_evals']} flux-line-drift evals)")
    print(f"\n  {'':<12}{'D_dir [mrad]':>14}{'iron <mu_r> sat':>17}"
          f"{'corner kappa':>14}")
    print(f"  {'FLAT cut':<12}{f['D_dir']*1e3:>14.3f}{f['mur_sat']:>17.0f}"
          f"{f['corner_kappa_sat']:>14.2f}")
    print(f"  {'OPTIMIZED':<12}{o['D_dir']*1e3:>14.3f}{o['mur_sat']:>17.0f}"
          f"{o['corner_kappa_sat']:>14.2f}")
    print(f"\n  => optimized end chamfer: depth = {o['depth_m']*1e3:.2f} mm, "
          f"exponent = {o['exponent']:.2f}")
    print(f"  => flux lines stay invariant {res['invariance_factor']:.1f}x deeper into "
          f"saturation (direction drift {f['D_dir']*1e3:.3f} -> {o['D_dir']*1e3:.3f} mrad "
          f"at B={res['B_sat_T']:.1f} T)")
    print(f"  => even the FLAT cut is already nearly invariant "
          f"(D_dir {f['D_dir']*1e3:.2f} mrad < 0.6 deg): "
          f"{res['flat_flux_lines_already_invariant']}")
    print(f"\n  Physics: below the iron knee the magnet is LINEAR, so scaling the current")
    print(f"  scales B everywhere -- the flux LINES are identical (the linear control")
    print(f"  D_dir ~ 0 at every drive).  Saturation is the ONLY thing that rotates them,")
    print(f"  starting at the pole-tip corner; the end chamfer relieves that corner and")
    print(f"  keeps the field-line pattern invariant deeper into saturation.")

    flat_curve = invariance_curve(f["depth_m"], f["exponent"])
    best_curve = invariance_curve(o["depth_m"], o["exponent"])

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    save = {k: v for k, v in res.items() if not k.startswith("_")}
    save["flat_invariance_curve"] = flat_curve
    save["optimized_invariance_curve"] = best_curve
    with open(jpath, "w") as fh:
        json.dump(save, fh, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(res, flat_curve, best_curve)


if __name__ == "__main__":
    main()
