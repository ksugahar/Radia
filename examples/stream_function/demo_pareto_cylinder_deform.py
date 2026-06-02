"""Cylinder sheet-metal (bankin-ho / 板金) on the (homogeneity, peak-J) Pareto
front: IN-SURFACE axial bending -- the DOMINANT cylinder lever (opposite of the
planar out-of-surface bending).

A length-preserving axial reparametrisation  Z(z) = z + sum_k b_k sin(k pi
(z+L/2)/L)  redistributes the loop spacing ALONG the surface; the cylinder
radius is held EXACTLY fixed, so this is 100% genuine forming -- there is NO
standoff component to subtract (unlike the planar case, which needs a zero-mean
constraint).  Spreading loops where the current is hot lowers the local current
density.

Done correctly:
  * loop axial extent = local spacing dZ (np.gradient of Z),
  * peak K_phi = -dpsi/dz via the NON-UNIFORM axial derivative on Z,
  * seminorm S = WEIGHTED graph Laplacian (z-edge weight a*dphi/dZ, phi-edge
    weight dZ/(a*dphi)) so the folded-Tikhonov min-seminorm solve is the true
    minimum surface-current energy on the non-uniform grid.
Monotonicity dZ/dz>0 is enforced by a penalty.  CMA-ES optimises the shape per
homogeneity level (warm-started across alpha).  Inner psi solve = the library
RegularizedTSVD folded Tikhonov routine.

Result (a=15cm, L=50cm, Gx fingerprint, +-8cm DSV): in-surface axial bending
pushes the WHOLE front down ~ -10 to -25% (best ~ -25%), all genuine.  Compare
radial surface forming, which is WEAK on the cylinder (~ -3%) -- see notes in
docs/stream_function/regularization.md.

Outputs JSON + PNG next to this script.
"""
import os
import sys
import json
import argparse

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "src"))

from demo_sf_to_peec_gx import loop_corners, _loop_Hz, make_dsv
from radia.stream_function import aca_tsvd, RegularizedTSVD


class CylinderDeformPareto:
    def __init__(self, a=0.15, L=0.50, Gx=1.0, dsv=0.08, nphi=24, nz=17,
                 ndsv=5, n_basis=4, bcap=0.10):
        self.a, self.L, self.Gx = a, L, Gx
        self.nphi, self.nz, self.n_basis, self.bcap = nphi, nz, n_basis, bcap
        self.phi = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
        self.z0 = np.linspace(-L / 2, L / 2, nz)
        self.dphi = 2 * np.pi / nphi
        self.obs = make_dsv(dsv, ndsv)
        self.M = len(self.obs)
        self.B = Gx * self.obs[:, 0]

    def Zof(self, b):
        if b is None:
            return self.z0.copy()
        s = np.zeros(self.nz)
        for k, bk in enumerate(b, start=1):
            s += bk * np.sin(k * np.pi * (self.z0 + self.L / 2) / self.L)
        return self.z0 + s

    def Sweighted(self, Z):
        a, dphi, nz, nphi = self.a, self.dphi, self.nz, self.nphi
        dZ = np.gradient(Z)
        N = nz * nphi
        S = np.zeros((N, N))

        def idx(iz, ip):
            return iz * nphi + (ip % nphi)
        for iz in range(nz):
            for ip in range(nphi):
                c = idx(iz, ip)
                wphi = dZ[iz] / (a * dphi)
                for jp in (ip - 1, ip + 1):
                    S[c, c] += wphi
                    S[c, idx(iz, jp)] -= wphi
                for jz in (iz - 1, iz + 1):
                    if 0 <= jz < nz:
                        wz = a * dphi / abs(Z[jz] - Z[iz])
                        S[c, c] += wz
                        S[c, idx(jz, ip)] -= wz
        return S + 1e-6 * np.eye(N)

    def build(self, b):
        Z = self.Zof(b)
        dZ = np.gradient(Z)
        corners = []
        for q in range(self.nz):
            for pp in range(self.nphi):
                corners.append(loop_corners(self.a, self.phi[pp], Z[q],
                                            self.dphi, abs(dZ[q])))
        A = np.array([[_loop_Hz(self.obs[i], corners[j])
                       for j in range(len(corners))] for i in range(self.M)])
        S = self.Sweighted(Z)
        md = float(np.mean(np.diag(A.T @ A)))
        base = aca_tsvd(self.M, len(corners), lambda i, j: float(A[i, j]),
                        modes=self.M, kmax=min(self.M, len(corners)),
                        aca_eps=1e-10, method=3)
        reg = RegularizedTSVD.from_stiffness(base, S)

        def sp(lam):
            psi = reg.solve(self.B, alpha=lam * md)
            P = psi.reshape(self.nz, self.nphi)
            dpz = np.gradient(P, Z, axis=0)
            dpp = (np.roll(P, -1, 1) - np.roll(P, 1, 1)) / (2 * self.dphi)
            K = np.sqrt(dpz ** 2 + (dpp / self.a) ** 2)
            mf = float(np.linalg.norm(A @ psi - self.B)
                       / (np.linalg.norm(self.B) + 1e-30))
            return mf, float(K.max())
        return sp

    def _mono_penalty(self, b):
        dZ = np.gradient(self.Zof(b))
        return 1e6 * max(0.0, 0.01 - float(dZ.min()))

    def optimize_at(self, lam, warm, n_trials):
        import optuna
        from optuna.samplers import CmaEsSampler
        optuna.logging.set_verbosity(optuna.logging.WARNING)

        def obj(trial):
            b = [trial.suggest_float(f"b{k}", -self.bcap, self.bcap)
                 for k in range(self.n_basis)]
            mf, pk = self.build(b)(lam)
            trial.set_user_attr("mf", mf)
            return pk + self._mono_penalty(b)
        st = optuna.create_study(
            sampler=CmaEsSampler(seed=0, warn_independent_sampling=False))
        if warm:
            st.enqueue_trial(warm)
        st.optimize(obj, n_trials=n_trials, show_progress_bar=False)
        bt = st.best_trial
        return bt.user_attrs["mf"], bt.value, bt.params


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radius", type=float, default=0.15)
    ap.add_argument("--length", type=float, default=0.50)
    ap.add_argument("--n-trials", type=int, default=110,
                    help="CMA-ES trials per homogeneity level")
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke: 3 alphas x 30 trials")
    ap.add_argument("--out", default=os.path.join(_HERE,
                                                  "pareto_cylinder_deform"))
    args = ap.parse_args()

    P = CylinderDeformPareto(a=args.radius, L=args.length)
    flat = P.build(None)
    lams_flat = np.concatenate([[0.0], np.logspace(-4, 1.5, 9)])
    flat_front = [flat(l) for l in lams_flat]
    print(f"FLAT cylinder exact-homog peak = {flat_front[0][1]:.3e}  "
          f"(a={args.radius} L={args.length}, in-surface axial bending)\n",
          flush=True)

    fm = np.array([p[0] for p in flat_front])
    fp = np.array([p[1] for p in flat_front])
    o = np.argsort(fm)
    lams = [0.0, 3e-2, 3e0] if args.quick else [0.0, 3e-3, 3e-2, 3e-1, 3e0]
    n_trials = 30 if args.quick else args.n_trials

    warm, opt_front = None, []
    print(f"{'misfit%':>9} {'flat_pk':>10} {'opt_pk':>10} {'vs_flat':>9}",
          flush=True)
    for lam in lams:
        mf, pk, bp = P.optimize_at(lam, warm, n_trials)
        warm = bp
        opt_front.append((mf, pk))
        pf = float(np.interp(mf, fm[o], fp[o]))
        print(f"{mf*100:9.3f} {pf:10.3e} {pk:10.3e} {100*(pk/pf-1):+8.1f}%",
              flush=True)

    json.dump({"radius": args.radius, "length": args.length,
               "flat_front": flat_front, "opt_front": opt_front},
              open(args.out + ".json", "w"), indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.plot([p[0] * 100 for p in flat_front], [p[1] for p in flat_front],
            "o-", color="#c0392b", label="flat cylinder (baseline)")
    ax.plot([p[0] * 100 for p in opt_front], [p[1] for p in opt_front],
            "s-", color="#1e8449",
            label="in-surface axial bending (radius fixed)")
    ax.set_xscale("log")
    ax.set_xlabel("inhomogeneity ||A psi - B|| / ||B|| [%]")
    ax.set_ylabel("peak surface current density max|grad psi|")
    ax.set_title("Cylinder sheet-metal: IN-SURFACE axial bending pushes the\n"
                 "front (radius fixed = 100% genuine forming, no standoff)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=130)
    print(f"\nsaved {args.out}.json / .png", flush=True)


if __name__ == "__main__":
    main()
