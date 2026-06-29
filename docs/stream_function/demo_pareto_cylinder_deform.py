"""Cylinder sheet-metal (bankin-ho / 板金) on the (homogeneity, peak-J) Pareto
front: IN-SURFACE bending -- the DOMINANT cylinder lever (opposite of the
planar out-of-surface bending).

In-surface reparametrisation redistributes the loop spacing ALONG the surface
at FIXED radius, so it is 100% genuine forming with NO standoff component
(unlike the planar case which needs a zero-mean constraint).  Two directions:

  AXIAL    Z(z)   = z + sum_k b_k sin(k pi (z+L/2)/L)     (length-preserving)
  AZIMUTHAL Phi(p)= phi_p + sum_m c_m sin(m phi_p)        (periodic, --azimuthal)

The in-surface lever DIRECTION follows the TARGET symmetry (azimuthal order m):
  target          m    axial    azimuthal benefit (2D vs axial, exact homog)
  Gx = x          1   -16%      +0.0%   (smooth cos phi, no azimuthal hot spot)
  Z2 = 2z^2-r^2   0   -20%      -1.8%   (axisymmetric -> axial dominates)
  S2 = xy         2   -11%      -1.2%
  C2 = x^2-y^2    2    -5%      -5.4%   (2D ~DOUBLES the benefit to ~ -10%)
So AXIAL-only is the right default for Gx / Z2; turn on --azimuthal for the
high-m azimuthal shims (C2 ellipse / S2).

Done correctly: loop axial extent = local spacing dZ; peak K_phi = -dpsi/dz via
the NON-UNIFORM axial derivative and K_z = (1/a) dpsi/dPhi via the PERIODIC
non-uniform azimuthal derivative; seminorm S = spacing-WEIGHTED graph Laplacian
(z-edge a*dphi/dZ, phi-edge dZ/(a*dPhi)).  Monotone reparametrisation enforced.
Inner psi solve = the library RegularizedTSVD folded Tikhonov routine.

Outputs JSON + PNG next to this script.
"""
import os
import sys
import json
import argparse

import numpy as np

from radia_mcp.figure import lab_figure, save_lab_figure

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "src"))

from demo_sf_to_peec_gx import loop_corners, _loop_Hz, make_dsv
from radia.stream_function import aca_tsvd, RegularizedTSVD

TWO_PI = 2.0 * np.pi


def target_vector(name, obs):
    x, y, z = obs[:, 0], obs[:, 1], obs[:, 2]
    return {"gx": x, "c2": x ** 2 - y ** 2, "s2": x * y,
            "z2": 2 * z ** 2 - x ** 2 - y ** 2}[name]


class CylinderDeformPareto:
    def __init__(self, a=0.15, L=0.50, target="gx", dsv=0.08, nphi=24, nz=17,
                 ndsv=5, n_axial=4, n_azim=4, bcap=0.10, ccap=0.30):
        self.a, self.L = a, L
        self.nphi, self.nz = nphi, nz
        self.n_axial, self.n_azim = n_axial, n_azim
        self.bcap, self.ccap = bcap, ccap
        self.phi0 = np.linspace(0.0, TWO_PI, nphi, endpoint=False)
        self.z0 = np.linspace(-L / 2, L / 2, nz)
        self.obs = make_dsv(dsv, ndsv)
        self.M = len(self.obs)
        self.B = target_vector(target, self.obs)

    def Zof(self, b):
        s = np.zeros(self.nz)
        for k, bk in enumerate(b, start=1):
            s += bk * np.sin(k * np.pi * (self.z0 + self.L / 2) / self.L)
        return self.z0 + s

    def Phiof(self, c):
        if c is None:
            return self.phi0.copy()
        s = np.zeros(self.nphi)
        for m, cm in enumerate(c, start=1):
            s += cm * np.sin(m * self.phi0)
        return self.phi0 + s

    def build(self, b, c):
        a, dphi0 = self.a, TWO_PI / self.nphi
        Z = self.Zof(b) if b is not None else self.z0.copy()
        Phi = self.Phiof(c)
        dZ = np.gradient(Z)
        fwd = (np.roll(Phi, -1) - Phi) % TWO_PI
        bwd = (Phi - np.roll(Phi, 1)) % TWO_PI
        dPhi = 0.5 * (fwd + bwd)
        corners = []
        for q in range(self.nz):
            for p in range(self.nphi):
                corners.append(loop_corners(a, Phi[p], Z[q], dPhi[p],
                                            abs(dZ[q])))
        A = np.array([[_loop_Hz(self.obs[i], corners[j])
                       for j in range(len(corners))] for i in range(self.M)])
        N = self.nz * self.nphi
        S = np.zeros((N, N))

        def idx(iz, ip):
            return iz * self.nphi + (ip % self.nphi)
        for iz in range(self.nz):
            for ip in range(self.nphi):
                cc = idx(iz, ip)
                for jp, arc in ((ip + 1, fwd[ip]), (ip - 1, bwd[ip])):
                    wphi = abs(dZ[iz]) / (a * arc)
                    S[cc, cc] += wphi
                    S[cc, idx(iz, jp)] -= wphi
                for jz in (iz - 1, iz + 1):
                    if 0 <= jz < self.nz:
                        wz = a * dPhi[ip] / abs(Z[jz] - Z[iz])
                        S[cc, cc] += wz
                        S[cc, idx(jz, ip)] -= wz
        S += 1e-6 * np.eye(N)
        md = float(np.mean(np.diag(A.T @ A)))
        base = aca_tsvd(self.M, len(corners), lambda i, j: float(A[i, j]),
                        modes=self.M, kmax=min(self.M, len(corners)),
                        aca_eps=1e-10, method=3)
        reg = RegularizedTSVD.from_stiffness(base, S)
        den_phi = (np.roll(Phi, -1) - np.roll(Phi, 1)) % TWO_PI

        def sp(lam):
            psi = reg.solve(self.B, alpha=lam * md)
            P = psi.reshape(self.nz, self.nphi)
            dpz = np.gradient(P, Z, axis=0)
            dpp = (np.roll(P, -1, 1) - np.roll(P, 1, 1)) / den_phi[None, :]
            K = np.sqrt(dpz ** 2 + (dpp / a) ** 2)
            mf = float(np.linalg.norm(A @ psi - self.B)
                       / (np.linalg.norm(self.B) + 1e-30))
            return mf, float(K.max())
        return sp

    def _mono_penalty(self, b, c):
        p = 1e6 * max(0.0, 0.01 - float(np.gradient(self.Zof(b)).min()))
        if c is not None:
            Phi = self.Phiof(c)
            dmin = float(((np.roll(Phi, -1) - Phi) % TWO_PI).min())
            p += 1e6 * max(0.0, 0.01 - dmin)
        return p

    def optimize_at(self, lam, azimuthal, warm, n_trials):
        import optuna
        from optuna.samplers import CmaEsSampler
        optuna.logging.set_verbosity(optuna.logging.WARNING)

        def obj(trial):
            b = [trial.suggest_float(f"b{k}", -self.bcap, self.bcap)
                 for k in range(self.n_axial)]
            c = ([trial.suggest_float(f"c{m}", -self.ccap, self.ccap)
                  for m in range(self.n_azim)] if azimuthal else None)
            mf, pk = self.build(b, c)(lam)
            trial.set_user_attr("mf", mf)
            return pk + self._mono_penalty(b, c)
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
    ap.add_argument("--target", choices=["gx", "c2", "s2", "z2"], default="gx",
                    help="DSV target: gx=x (m=1), c2=x^2-y^2 (ellipse, m=2), "
                         "s2=xy (m=2), z2=2z^2-r^2 (axisymmetric, m=0)")
    ap.add_argument("--azimuthal", action="store_true",
                    help="enable the AZIMUTHAL reparametrisation (2D in-surface);"
                         " pays off for high-m targets (c2/s2), null for gx/z2")
    ap.add_argument("--radius", type=float, default=0.15)
    ap.add_argument("--length", type=float, default=0.50)
    ap.add_argument("--n-trials", type=int, default=110)
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke: 3 alphas x 30 trials")
    ap.add_argument("--out", default=os.path.join(_HERE,
                                                  "pareto_cylinder_deform"))
    args = ap.parse_args()

    P = CylinderDeformPareto(a=args.radius, L=args.length, target=args.target)
    flat = P.build(None, None)
    lams_flat = np.concatenate([[0.0], np.logspace(-4, 1.5, 9)])
    flat_front = [flat(l) for l in lams_flat]
    mode = "axial+azimuthal (2D)" if args.azimuthal else "axial-only"
    print(f"FLAT cylinder exact-homog peak = {flat_front[0][1]:.3e}  "
          f"(target={args.target}, {mode})\n", flush=True)

    fm = np.array([p[0] for p in flat_front])
    fp = np.array([p[1] for p in flat_front])
    o = np.argsort(fm)
    lams = [0.0, 3e-2, 3e0] if args.quick else [0.0, 3e-3, 3e-2, 3e-1, 3e0]
    n_trials = 30 if args.quick else args.n_trials

    warm, opt_front = None, []
    print(f"{'misfit%':>9} {'flat_pk':>10} {'opt_pk':>10} {'vs_flat':>9}",
          flush=True)
    for lam in lams:
        mf, pk, bp = P.optimize_at(lam, args.azimuthal, warm, n_trials)
        warm = bp
        opt_front.append((mf, pk))
        pf = float(np.interp(mf, fm[o], fp[o]))
        print(f"{mf*100:9.3f} {pf:10.3e} {pk:10.3e} {100*(pk/pf-1):+8.1f}%",
              flush=True)

    json.dump({"target": args.target, "azimuthal": args.azimuthal,
               "radius": args.radius, "length": args.length,
               "flat_front": flat_front, "opt_front": opt_front},
              open(args.out + ".json", "w"), indent=2)

    import matplotlib
    matplotlib.use("Agg")
    fig, ax = lab_figure(embed_width_cm=8.0)
    ax.plot([p[0] * 100 for p in flat_front], [p[1] for p in flat_front],
            "o-", color="#c0392b", label="flat cylinder (baseline)")
    ax.plot([p[0] * 100 for p in opt_front], [p[1] for p in opt_front],
            "s-", color="#1e8449",
            label=f"in-surface bending ({mode}, radius fixed)")
    ax.set_xscale("log")
    ax.set_xlabel("inhomogeneity ||A psi - B|| / ||B|| (%)")
    ax.set_ylabel("peak surface current density max|grad psi|")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    save_lab_figure(fig, args.out, embed_width_cm=8.0)
    print(f"\nsaved {args.out}.json / .png", flush=True)


if __name__ == "__main__":
    main()
