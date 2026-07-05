"""Pareto optimisation via folded Tikhonov + ACA+TSVD.

ONE ACA+TSVD factorisation of the Biot-Savart matrix A, then the WHOLE
Pareto front of

    f1 = field inhomogeneity   ||A psi - B|| / ||B||      (data misfit)
    f2 = peak current density   max_v |grad psi|           (= ||K||_inf)

is swept through the folded Tikhonov core

    psi(alpha) = S^-1 V . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B,
                 W = V^T S^-1 V   (k x k),     alpha = lam * mean(diag(A^T A)).

WHY THIS IS PARETO-OPTIMAL.  Tikhonov
    min ||A psi - B||^2 + alpha psi^T S psi
is the weighted-sum scalarisation of the convex bi-objective
(misfit, seminorm).  Sweeping alpha in (0, inf) traces the ENTIRE convex
Pareto front of (misfit, seminorm).  With S = H1 stiffness the seminorm
is the L2 surface-current energy; the peak (L-inf) follows along.

WHY FOLD IT INTO ACA+TSVD.  The expensive parts -- the ACA+ factorisation
A = U Sigma V^T and the fold S^-1 V, W -- are computed ONCE.  Each alpha
point is then a single k x k core solve (~50 us), so the whole L-curve is
sub-millisecond of linear algebra.  When each A(i, j) is a material-kernel
evaluation (Radia HDiv-VIM solve + rad.Fld), the ACA+ factorisation
dominates and reusing it across the entire alpha-sweep is the whole
performance argument.

TWO FRONTS:
  * --front energy : pure Tikhonov alpha-sweep, S = H1.  The textbook
    L-curve = exact Pareto front of (misfit, energy).  Fold computed once.
  * --front peak   : Tikhonov alpha-sweep with an L-inf IRLS-reweighted
    seminorm S^(w) -> minimises peak instead of energy -> DOMINATES the H1
    front in peak (median ~ -18 % peak at matched misfit).  Re-uses the
    SAME cached ACA (U, Sigma, V); only the fold S^-1 V is rebuilt.
  * --front both   : overlay (default).

See docs/stream_function/regularization.md  "Tikhonov is the SAME core
with + alpha I"  for the derivation, and the verification that the folded
form equals the dense (A^T A + alpha S)^-1 A^T B (machine precision for
S = I, ACA tolerance for general S).

Outputs JSON (all sweep points + non-dominated set + solve timing) and PNG
next to this script.

Usage:
    python demo_pareto_tikhonov_aca.py                 # both fronts, gradient
    python demo_pareto_tikhonov_aca.py --front energy  # pure Tikhonov L-curve
    python demo_pareto_tikhonov_aca.py --target uniform --n-lam 20
"""
import os
import sys
import json
import argparse

import numpy as np

from radia_mcp.figure import lab_figure, save_lab_figure

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                    # sibling demos
sys.path.insert(0, os.path.join(_HERE, "..", "..", "src"))   # radia package

from demo_planar_uniform_fem_psi import build_fem_matrix, build_h1_stiffness
from radia.stream_function import aca_tsvd, RegularizedTSVD
from ngsolve import (BilinearForm, GridFunction, grad, dx, TaskManager,
                     H1 as H1_space)


# ----------------------------------------------------------------------
# Engine: ACA+TSVD once, folded Tikhonov core re-used across alpha / S
# ----------------------------------------------------------------------
class FoldedPareto:
    """Folded Tikhonov + ACA+TSVD Pareto engine (one factorisation)."""

    def __init__(self, plane_half=0.25, maxh=0.025, order=3, target="gradient",
                 G=1.0, th=0.05, tz=0.10, n_dsv=5, aca_eps=1.0e-10):
        gg = np.linspace(-th, th, n_dsv)
        Xg, Yg = np.meshgrid(gg, gg, indexing="ij")
        self.targets = np.column_stack(
            [Xg.ravel(), Yg.ravel(), tz * np.ones(Xg.size)])
        if target == "gradient":
            self.B = G * self.targets[:, 0]            # Bz = G * x (hot spots)
        else:
            self.B = G * np.ones(len(self.targets))    # uniform Bz

        # Caller wraps TaskManager (per Radia TaskManager-Only policy).
        with TaskManager():
            A, self.fes, self.mesh, _ = build_fem_matrix(
                plane_half, maxh, order, "boundary", self.targets)
            self.fi = np.where(np.array(self.fes.FreeDofs()))[0]
            self.Af = A[:, self.fi]
            self.M, self.N = self.Af.shape
            self.S0 = build_h1_stiffness(self.fes, self.fi)    # H1 stiffness
        self.AtA = self.Af.T @ self.Af
        self.md = float(np.mean(np.diag(self.AtA)))            # alpha scale
        self.nv = self.mesh.nv
        self.vpts = np.array([list(v.point) for v in self.mesh.vertices])

        # ---- the ONE ACA+TSVD factorisation ----
        base = aca_tsvd(self.M, self.N, lambda i, j: float(self.Af[i, j]),
                        modes=self.M, kmax=min(self.M, self.N),
                        aca_eps=aca_eps)
        self.k = base.modes
        self.U = base.U[:, :self.k]
        self.Sig = base.S[:self.k]
        self.V = base.V[:, :self.k]
        self.UtB = self.U.T @ self.B
        self.base = base

        # NGSolve scratch for weighted stiffness / gradient sampling
        self.u_t, self.v_t = self.fes.TnT()
        self.fes_w = H1_space(self.mesh, order=1)
        self.gf_w = GridFunction(self.fes_w)
        self.gf = GridFunction(self.fes)

    # ---- fold a seminorm S into the library routine (one S^-1 V + W) ----
    def fold(self, S):
        # RegularizedTSVD.from_stiffness is THE library routine: it does
        # the S^-1 V solve and caches W = V^T S^-1 V on the shared base.
        return RegularizedTSVD.from_stiffness(self.base, S)

    # ---- folded Tikhonov core via the library: ~50 us per alpha ----
    def solve(self, lam, reg):
        # reg.solve(B, alpha) = Sinv_V . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B
        return reg.solve(self.B, alpha=lam * self.md)   # psi_free

    # ---- objectives ----
    def grad_vert(self, psi_free):
        self.gf.vec.FV().NumPy()[:] = 0.0
        self.gf.vec.FV().NumPy()[self.fi] = psi_free
        gv = np.array(grad(self.gf)(self.mesh(self.vpts[:, 0],
                                              self.vpts[:, 1])))
        return np.sqrt(gv[:, 0] ** 2 + gv[:, 1] ** 2)

    def objectives(self, psi_free):
        misfit = float(np.linalg.norm(self.Af @ psi_free - self.B)
                       / (np.linalg.norm(self.B) + 1e-30))
        peak = float(self.grad_vert(psi_free).max())
        energy = float(psi_free @ (self.S0 @ psi_free))
        return misfit, peak, energy

    # ---- weighted H1 stiffness for the L-inf IRLS reweight ----
    def weighted_stiffness(self, w_vert):
        self.gf_w.vec.FV().NumPy()[:self.nv] = w_vert
        a = BilinearForm(self.fes, symmetric=True)
        a += self.gf_w * grad(self.u_t) * grad(self.v_t) * dx
        a += 1e-12 * self.u_t * self.v_t * dx
        with TaskManager():
            a.Assemble()
        Sd = np.array(a.mat.ToDense())
        return Sd[np.ix_(self.fi, self.fi)]

    # ---- Mode A: pure Tikhonov L-curve (S = H1) ----
    def front_energy(self, lams):
        reg = self.fold(self.S0)                        # fold ONCE (library)
        pts = []
        for lam in lams:
            psi = self.solve(lam, reg)
            mf, pk, en = self.objectives(psi)
            pts.append(dict(lam=float(lam), misfit=mf, peak=pk, energy=en))
        return pts, reg

    # ---- micro-bench: folded SOLVE only (the actual Tikhonov+ACA cost) ----
    def bench_solve(self, reg, reps=2000):
        import time
        self.solve(1.0, reg)                            # warm
        t0 = time.perf_counter()
        for _ in range(reps):
            self.solve(1.0, reg)
        return (time.perf_counter() - t0) / reps

    # ---- Mode B: Tikhonov + L-inf IRLS reweighted S (min-peak front) ----
    def front_peak(self, lams, n_iter=14, ratio_max=20.0, damping=0.5,
                   hot_q=0.85, hot_boost=3.0):
        pts = []
        for lam in lams:
            w = np.ones(self.nv)
            psi = None
            for _ in range(n_iter):
                reg = self.fold(self.weighted_stiffness(w))
                psi = self.solve(lam, reg)
                g = self.grad_vert(psi)
                peak = g.max() + 1e-30
                rel = np.clip(g / peak, 1.0 / ratio_max, 1.0)
                w_new = rel ** 2.0
                w_new = np.where(g >= np.quantile(g, hot_q),
                                 w_new * hot_boost, w_new)
                w_new /= w_new.mean()
                w = (1.0 - damping) * w + damping * w_new
            mf, pk, en = self.objectives(psi)
            pts.append(dict(lam=float(lam), misfit=mf, peak=pk, energy=en))
        return pts


def pareto_nondominated(pts, x="misfit", y="peak"):
    """Lower-left non-dominated set (minimise both x and y)."""
    s = sorted(pts, key=lambda p: (p[x], p[y]))
    out, best = [], float("inf")
    for p in s:
        if p[y] < best - 1e-15:
            out.append(p)
            best = p[y]
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", choices=["gradient", "uniform"],
                    default="gradient",
                    help="DSV target: gradient Bz=G*x (hot spots, default) "
                         "or uniform Bz=G")
    ap.add_argument("--front", choices=["energy", "peak", "both"],
                    default="both",
                    help="energy = pure Tikhonov L-curve (S=H1); peak = "
                         "L-inf IRLS-reweighted (min-peak); both = overlay")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--maxh", type=float, default=0.025)
    ap.add_argument("--plane-half", type=float, default=0.25)
    ap.add_argument("--lam-min", type=float, default=1e-4,
                    help="min ridge lam (alpha = lam*mean(diag(A^T A)))")
    ap.add_argument("--lam-max", type=float, default=3e1)
    ap.add_argument("--n-lam", type=int, default=14)
    ap.add_argument("--out", default=os.path.join(_HERE, "pareto_tikhonov_aca"),
                    help="output basename (.json + .png)")
    args = ap.parse_args()

    P = FoldedPareto(plane_half=args.plane_half, maxh=args.maxh,
                     order=args.order, target=args.target)
    print(f"target={args.target}  M={P.M}  N={P.N}  ACA k={P.k}  "
          f"sigma in [{P.Sig.min():.2e},{P.Sig.max():.2e}]", flush=True)

    lams = np.concatenate([[0.0],
                           np.logspace(np.log10(args.lam_min),
                                       np.log10(args.lam_max),
                                       args.n_lam - 1)])

    result = {"target": args.target, "M": P.M, "N": P.N, "aca_k": int(P.k),
              "md": P.md, "lams": lams.tolist()}

    print(f"\n{'mode':6} {'lam':>9} {'misfit%':>9} {'peak':>11} {'energy':>11}",
          flush=True)
    fE = fP = None
    if args.front in ("energy", "both"):
        fE, fold0 = P.front_energy(lams)
        result["front_energy_h1"] = fE
        for p in fE:
            print(f"{'H1':6} {p['lam']:>9.1e} {p['misfit']*100:>9.2e} "
                  f"{p['peak']:>11.3e} {p['energy']:>11.3e}", flush=True)
        t_solve = P.bench_solve(fold0)
        result["t_solve_us"] = t_solve * 1e6
        print(f"   [folded Tikhonov SOLVE = {t_solve*1e6:.1f} us/alpha-point "
              f"after ONE ACA factorisation + ONE fold; whole "
              f"{len(lams)}-point L-curve = {t_solve*len(lams)*1e6:.0f} us]",
              flush=True)
        print("   [peak-objective (NGSolve vertex grad sampling) is a separate "
              "MEASUREMENT cost ~45 ms/point, vectorisable; not the solve]",
              flush=True)
    if args.front in ("peak", "both"):
        fP = P.front_peak(lams)
        result["front_peak_linf"] = fP
        for p in fP:
            print(f"{'Linf':6} {p['lam']:>9.1e} {p['misfit']*100:>9.2e} "
                  f"{p['peak']:>11.3e} {p['energy']:>11.3e}", flush=True)

    if args.front == "both":
        result["nondominated_h1"] = pareto_nondominated(fE)
        result["nondominated_linf"] = pareto_nondominated(fP)
        mE = np.array([p["misfit"] for p in fE])
        kE = np.array([p["peak"] for p in fE])
        o = np.argsort(mE)
        pushes = [100 * (1 - p["peak"]
                         / float(np.interp(p["misfit"], mE[o], kE[o])))
                  for p in fP]
        result["peak_push_pct_vs_h1"] = pushes
        print(f"\nL-inf push vs H1 at matched misfit: "
              f"median {np.median(pushes):+.0f}%, "
              f"range [{min(pushes):+.0f}%, {max(pushes):+.0f}%]", flush=True)

    with open(args.out + ".json", "w") as f:
        json.dump(result, f, indent=2)

    # ---- plot: non-dominated envelope (bold) over raw sweep (faint) ----
    import matplotlib
    matplotlib.use("Agg")
    fig, ax = lab_figure(embed_width_cm=8.0)

    def draw(pts, color, marker, label):
        m = np.array([p["misfit"] * 100 for p in pts])
        k = np.array([p["peak"] / 1e6 for p in pts])
        ax.scatter(m, k, s=18, color=color, alpha=0.35)
        nd = pareto_nondominated(pts)
        mn = np.array([p["misfit"] * 100 for p in nd])
        kn = np.array([p["peak"] / 1e6 for p in nd])
        ax.plot(mn, kn, marker + "-", color=color, label=label)

    if fE is not None:
        draw(fE, "#c0392b", "o", "Tikhonov L-curve, S=H1 (min energy)")
    if fP is not None:
        draw(fP, "#2471a3", "s", "Tikhonov + L-inf IRLS (min peak)")
    ax.set_xscale("log")
    ax.set_xlabel("field inhomogeneity  ||A psi - B|| / ||B||  (%)")
    ax.set_ylabel(r"peak $|\nabla\psi|$ ($\times10^6$ A/m)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    save_lab_figure(fig, args.out, embed_width_cm=8.0)
    print(f"\nsaved {args.out}.json / .png", flush=True)


if __name__ == "__main__":
    main()
