"""Sheet-metal (bankin-ho / 板金) deformation pushes the (homogeneity, peak-J)
Pareto front.

The conductor SURFACE is formed z = f(x,y) (entering the Biot-Savart kernel via
build_fem_matrix_deformed) and the shape is optimised, per homogeneity level, to
minimise the peak surface current density max|grad psi| -- with the folded
Tikhonov + ACA+TSVD routine (RegularizedTSVD) doing the inner psi solve.

HONEST decomposition.  A deformation lowers the peak two ways:
  * STANDOFF -- forming the sheet CLOSER to the target (smaller average gap):
    a large effect, but achievable by simply repositioning a FLAT former, so
    NOT genuine forming.
  * BENDING  -- reshaping at FIXED average standoff (zero-mean deformation):
    the genuine sheet-metal contribution.
This demo defaults to --zero-mean (pure bending; report the genuine number).
--allow-standoff lifts the constraint (bigger but standoff-dominated number).

Deformation basis (y-even, matching the Bz = Gx target symmetry), localised by a
Gaussian envelope so the footprint boundary is fixed:

  z(x,y) = env(x,y;s) * (c0 + c1 xn + c2 xn^2 + c3 xn^3 + c4 yn^2 + c5 xn yn^2)
           [- mean, if --zero-mean]            xn = x/ph,  yn = y/ph

Result (planar gradient hot-spot, order-3 psi, |z| <= 5 cm): genuine zero-mean
bending pushes the WHOLE front down ~ -5 to -18% (best -17% near exact
homogeneity); allowing standoff reaches ~ -53% but that is mostly repositioning.

Outputs JSON + PNG next to this script.
"""
import os
import json
import argparse

import numpy as np

from _validation_output import validation_json_for_basename

from radia_mcp.figure import lab_figure, save_lab_figure

_HERE = os.path.dirname(os.path.abspath(__file__))

from demo_planar_uniform_fem_psi_advanced import build_fem_matrix_deformed
from demo_planar_uniform_fem_psi import build_h1_stiffness
from radia.stream_function import aca_tsvd, RegularizedTSVD
from ngsolve import grad, GridFunction, exp, TaskManager


class DeformPareto:
    def __init__(self, plane_half=0.25, maxh=0.03, order=3, th=0.05, tz=0.10,
                 G=1.0, zero_mean=True, zcap=0.05):
        self.ph, self.maxh, self.order = plane_half, maxh, order
        self.zero_mean, self.zcap = zero_mean, zcap
        g = np.linspace(-th, th, 5)
        Xt, Yt = np.meshgrid(g, g, indexing="ij")
        self.targets = np.column_stack(
            [Xt.ravel(), Yt.ravel(), tz * np.ones(Xt.size)])
        self.B = G * self.targets[:, 0]
        gx = np.linspace(-plane_half, plane_half, 61)
        self._GX, self._GY = np.meshgrid(gx, gx)

    def _znp(self, c, s):
        env = np.exp(-(self._GX ** 2 + self._GY ** 2) / (2 * s * s))
        xn, yn = self._GX / self.ph, self._GY / self.ph
        z = (c[0] + c[1] * xn + c[2] * xn ** 2 + c[3] * xn ** 3
             + c[4] * yn ** 2 + c[5] * xn * yn ** 2) * env
        return z - z.mean() if self.zero_mean else z

    def zfn(self, c, s):
        if c is None:
            return None
        ph = self.ph
        if self.zero_mean:
            env = np.exp(-(self._GX ** 2 + self._GY ** 2) / (2 * s * s))
            xn, yn = self._GX / ph, self._GY / ph
            poly = (c[0] + c[1] * xn + c[2] * xn ** 2 + c[3] * xn ** 3
                    + c[4] * yn ** 2 + c[5] * xn * yn ** 2)
            zmean = float((env * poly).mean())
        else:
            zmean = 0.0

        def f(xx, yy):
            env = exp(-(xx * xx + yy * yy) / (2 * s * s))
            xn, yn = xx / ph, yy / ph
            z = (c[0] + c[1] * xn + c[2] * xn ** 2 + c[3] * xn ** 3
                 + c[4] * yn ** 2 + c[5] * xn * yn ** 2) * env
            return z - zmean
        return f

    def zmax(self, c, s):
        return float(np.max(np.abs(self._znp(c, s))))

    def engine(self, surface_z_fn):
        with TaskManager():
            A, fes, mesh, _ = build_fem_matrix_deformed(
                self.ph, self.maxh, self.order, "boundary", self.targets,
                surface_z_fn)
            fi = np.where(np.array(fes.FreeDofs()))[0]
            Af = A[:, fi]
            S = build_h1_stiffness(fes, fi)
        md = float(np.mean(np.diag(Af.T @ Af)))
        base = aca_tsvd(len(self.targets), len(fi),
                        lambda i, j: float(Af[i, j]), modes=len(self.targets),
                        kmax=min(len(self.targets), len(fi)), aca_eps=1e-10)
        reg = RegularizedTSVD.from_stiffness(base, S)
        vpts = np.array([list(v.point) for v in mesh.vertices])
        gf = GridFunction(fes)

        def sp(lam):
            psi = reg.solve(self.B, alpha=lam * md)
            gf.vec.FV().NumPy()[:] = 0.0
            gf.vec.FV().NumPy()[fi] = psi
            gv = np.array(grad(gf)(mesh(vpts[:, 0], vpts[:, 1])))
            pk = float(np.sqrt(gv[:, 0] ** 2 + gv[:, 1] ** 2).max())
            mf = float(np.linalg.norm(Af @ psi - self.B)
                       / (np.linalg.norm(self.B) + 1e-30))
            return mf, pk
        return sp

    def optimize_at(self, lam, warm, n_trials):
        import optuna
        from optuna.samplers import CmaEsSampler
        optuna.logging.set_verbosity(optuna.logging.WARNING)

        def obj(trial):
            c = [trial.suggest_float(f"c{i}", -0.06, 0.06) for i in range(6)]
            s = trial.suggest_float("s", 0.06, 0.30)
            zmax = self.zmax(c, s)
            mf, pk = self.engine(self.zfn(c, s))(lam)
            trial.set_user_attr("mf", mf)
            return pk + 1e6 * max(0.0, zmax - self.zcap) / self.zcap
        st = optuna.create_study(
            sampler=CmaEsSampler(seed=0, warn_independent_sampling=False))
        if warm is not None:
            st.enqueue_trial(warm)
        st.optimize(obj, n_trials=n_trials, show_progress_bar=False)
        bt = st.best_trial
        return bt.user_attrs["mf"], bt.value, bt.params


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--allow-standoff", action="store_true",
                    help="lift the zero-mean constraint (bigger but "
                         "standoff-dominated improvement)")
    ap.add_argument("--zcap", type=float, default=0.05, help="|z| cap [m]")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--maxh", type=float, default=0.03)
    ap.add_argument("--n-trials", type=int, default=110,
                    help="CMA-ES trials per homogeneity level")
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke: 3 alphas x 30 trials")
    ap.add_argument("--out", default=os.path.join(_HERE, "pareto_deform"))
    args = ap.parse_args()

    P = DeformPareto(maxh=args.maxh, order=args.order,
                     zero_mean=not args.allow_standoff, zcap=args.zcap)

    lams_flat = np.concatenate([[0.0], np.logspace(-4, 1.5, 9)])
    flat = P.engine(None)
    flat_front = [flat(l) for l in lams_flat]
    print(f"FLAT exact-homog peak = {flat_front[0][1]:.3e}  "
          f"(mode = {'zero-mean bending' if P.zero_mean else 'standoff-allowed'})\n",
          flush=True)

    # known good warm-start for the zero-mean planar case
    warm = dict(c0=-0.046, c1=-0.001, c2=-0.016, c3=0.0, c4=0.034, c5=0.0,
                s=0.105)
    lams = [0.0, 3e-2, 3e-1, 3e0] if args.quick else [0.0, 3e-3, 3e-2, 3e-1, 3e0]
    n_trials = 30 if args.quick else args.n_trials

    fm = np.array([p[0] for p in flat_front])
    fp = np.array([p[1] for p in flat_front])
    o = np.argsort(fm)
    opt_front = []
    print(f"{'misfit%':>9} {'flat_peak':>11} {'opt_peak':>11} {'vs_flat':>9}",
          flush=True)
    for lam in lams:
        mf, pk, bp = P.optimize_at(lam, warm, n_trials)
        warm = bp
        opt_front.append((mf, pk))
        pf = float(np.interp(mf, fm[o], fp[o]))
        print(f"{mf*100:9.3f} {pf:11.3e} {pk:11.3e} {100*(pk/pf-1):+8.1f}%",
              flush=True)

    json_path = validation_json_for_basename(args.out, "pareto_deform.json")
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump({"zero_mean": P.zero_mean, "zcap": args.zcap,
                   "flat_front": flat_front, "opt_front": opt_front},
                  stream, indent=2, allow_nan=False)

    import matplotlib
    matplotlib.use("Agg")
    fig, ax = lab_figure(embed_width_cm=8.0)
    ax.plot([p[0] * 100 for p in flat_front], [p[1] / 1e6 for p in flat_front],
            "o-", color="#c0392b", label="flat (baseline)")
    tag = "standoff-allowed" if not P.zero_mean else "zero-mean (genuine bending)"
    ax.plot([p[0] * 100 for p in opt_front], [p[1] / 1e6 for p in opt_front],
            "s-", color="#1e8449", label=f"sheet-metal deform, {tag}")
    ax.set_xscale("log")
    ax.set_xlabel("inhomogeneity ||A psi - B|| / ||B|| (%)")
    ax.set_ylabel(r"peak $|\nabla\psi|$ ($\times10^6$ A/m)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    save_lab_figure(fig, args.out, embed_width_cm=8.0)
    print(f"\nsaved {json_path} / {args.out}.png", flush=True)


if __name__ == "__main__":
    main()
