"""Push the (homogeneity, peak-J) Pareto front with the GEOMETRY lever, and
trace the joint front with NSGA-II.

Builds on demo_pareto_tikhonov_aca.py (the folded Tikhonov + ACA+TSVD routine).
The Tikhonov alpha-sweep moves ALONG the front for a FIXED geometry; here the
conductor geometry (planar former half-extent) becomes an extra design dof that
pushes the whole front DOWN -- a Pareto improvement A itself, not a
redistribution within a fixed A.

Two modes:
  --mode geometry : sweep former sizes; each former's alpha-sweep traces its
      front (A + ACA rebuilt per former).  The non-dominated envelope over all
      formers = the geometry-pushed front.  Bigger former -> more area for the
      current -> lower peak, with diminishing returns.
  --mode nsga     : NSGA-II over (former_size, log10 alpha) with THREE
      objectives (homogeneity, peak, former_size).  The optimiser autonomously
      traces the (homogeneity, peak, size) trade-off SURFACE -- lower peak at a
      given homogeneity REQUIRES a larger former, and the 3-objective front
      exposes exactly that coupling.  (optuna NSGAIISampler; A+ACA cached per
      distinct former.)
  --mode both     : run both (default).

Outputs JSON + PNG next to this script.
"""
import os
import json
import argparse

import numpy as np

from _validation_output import validation_json_for_basename

_HERE = os.path.dirname(os.path.abspath(__file__))

from demo_pareto_tikhonov_aca import FoldedPareto, pareto_nondominated


# ----------------------------------------------------------------------
# Mode geometry: former-size envelope (H1 Tikhonov front per former)
# ----------------------------------------------------------------------
def geometry_envelope(formers, lams, target="gradient", maxh=0.03, order=3):
    per_former, allpts = {}, []
    for ph in formers:
        P = FoldedPareto(plane_half=ph, maxh=maxh, order=order, target=target)
        fE, _ = P.front_energy(lams)                 # H1 Tikhonov front (fast)
        per_former[ph] = fE
        allpts += fE
        nd = pareto_nondominated(fE)
        print(f"former={ph:.2f} m  exact-homog peak={nd[0]['peak']:.3e}  "
              f"(misfit {nd[0]['misfit']*100:.1e}%)", flush=True)
    return per_former, pareto_nondominated(allpts)


# ----------------------------------------------------------------------
# Mode nsga: NSGA-II over (former, log alpha), 3 objectives
# ----------------------------------------------------------------------
def nsga_joint(former_lo, former_hi, n_trials, pop, target="gradient",
               maxh=0.03, order=3, quant=0.005, seed=0):
    import optuna
    from optuna.samplers import NSGAIISampler
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    cache = {}

    def engine(fr):
        if fr not in cache:
            P = FoldedPareto(plane_half=fr, maxh=maxh, order=order,
                             target=target)
            cache[fr] = (P, P.fold(P.S0))            # H1 fold once per former
        return cache[fr]

    def objective(trial):
        former = trial.suggest_float("former", former_lo, former_hi)
        log_lam = trial.suggest_float("log_lam", -5.0, 1.5)
        fr = round(former / quant) * quant
        P, reg = engine(fr)
        psi = P.solve(10.0 ** log_lam, reg)
        misfit, peak, _ = P.objectives(psi)
        return misfit, peak, former

    study = optuna.create_study(
        directions=["minimize", "minimize", "minimize"],
        sampler=NSGAIISampler(population_size=pop, seed=seed))
    study.optimize(objective, n_trials=n_trials, show_progress_bar=False)
    front = sorted([(t.values[0], t.values[1], t.values[2])
                    for t in study.best_trials], key=lambda r: r[0])
    print(f"NSGA-II: {len(study.trials)} trials, {len(front)} non-dominated, "
          f"{len(cache)} distinct formers built", flush=True)
    return front


def _lower_envelope(front):
    """best (lowest) peak at each misfit -- the 2D projection of the 3-obj set."""
    s = sorted(front, key=lambda r: (r[0], r[1]))
    out, best = [], float("inf")
    for mf, pk, fo in s:
        if pk < best - 1e-15:
            out.append((mf, pk, fo)); best = pk
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["geometry", "nsga", "both"],
                    default="both")
    ap.add_argument("--target", choices=["gradient", "uniform"],
                    default="gradient")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--maxh", type=float, default=0.03)
    ap.add_argument("--nsga-trials", type=int, default=360)
    ap.add_argument("--nsga-pop", type=int, default=24)
    ap.add_argument("--out", default=os.path.join(_HERE,
                                                  "pareto_geometry_nsga"))
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["pdf.fonttype"] = 42

    lams = np.concatenate([[0.0], np.logspace(-4, 1.5, 11)])
    result = {"target": args.target}

    do_geo = args.mode in ("geometry", "both")
    do_nsga = args.mode in ("nsga", "both")
    fig, axes = plt.subplots(1, 2 if (do_geo and do_nsga) else 1,
                             figsize=(11.5 if (do_geo and do_nsga) else 6.6, 5.0),
                             squeeze=False)
    axes = axes[0]
    ax_i = 0

    if do_geo:
        formers = [0.18, 0.22, 0.27, 0.34, 0.42]
        per_former, env = geometry_envelope(formers, lams, args.target,
                                            args.maxh, args.order)
        result["geometry_per_former"] = {str(k): v
                                         for k, v in per_former.items()}
        result["geometry_envelope"] = env
        ax = axes[ax_i]; ax_i += 1
        cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(formers)))
        for ph, col in zip(formers, cmap):
            nd = pareto_nondominated(per_former[ph])
            ax.plot([p["misfit"] * 100 for p in nd],
                    [p["peak"] / 1e6 for p in nd], "o-", color=col, ms=4,
                    lw=1.2, label=f"former {ph*100:.0f} cm")
        ax.plot([p["misfit"] * 100 for p in env],
                [p["peak"] / 1e6 for p in env], "k--", lw=2.0,
                label="envelope (geometry-pushed)")
        ax.set_xscale("log"); ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax.set_xlabel("inhomogeneity ||A psi-B||/||B|| (%)")
        ax.set_ylabel(r"peak $|\nabla\psi|$ ($\times10^6$ A/m)")
        if do_geo and do_nsga:
            ax.text(0.02, 0.96, "(a)", transform=ax.transAxes, fontsize=10)

    if do_nsga:
        front = nsga_joint(0.16, 0.45, args.nsga_trials, args.nsga_pop,
                           args.target, args.maxh, args.order)
        result["nsga_front"] = [{"misfit": m, "peak": p, "former": f}
                                for m, p, f in front]
        env2 = _lower_envelope(front)
        result["nsga_lower_envelope"] = [{"misfit": m, "peak": p, "former": f}
                                         for m, p, f in env2]
        ax = axes[ax_i]
        mf = np.array([r[0] * 100 for r in front])
        pk = np.array([r[1] / 1e6 for r in front])
        fo = np.array([r[2] * 100 for r in front])
        sc = ax.scatter(mf, pk, c=fo, cmap="plasma", s=26, edgecolor="k",
                        lw=0.3)
        ax.plot([r[0] * 100 for r in env2], [r[1] / 1e6 for r in env2],
                "k--", lw=2.0, label="lower envelope")
        cb = fig.colorbar(sc, ax=ax); cb.set_label("former size (cm)")
        ax.set_xscale("log"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax.set_xlabel("inhomogeneity ||A psi-B||/||B|| (%)")
        ax.set_ylabel(r"peak $|\nabla\psi|$ ($\times10^6$ A/m)")
        if do_geo and do_nsga:
            ax.text(0.02, 0.96, "(b)", transform=ax.transAxes, fontsize=10)

    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=130)
    json_path = validation_json_for_basename(args.out, "pareto_geometry_nsga.json")
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print(f"\nsaved {json_path} / {args.out}.png", flush=True)


if __name__ == "__main__":
    main()
