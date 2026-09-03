"""rotating_magnet_eddy.py -- moving-magnet eddy current: when is a per-step
FEM (or a CLN model-order reduction) actually needed?  The magnetic-Reynolds
crossover, on Yano's rotating-magnet-over-plate problem.

Geometry (Yano & Sugahara): a 1 mm cube permanent magnet (Br = 0.2 T,
M = Br/mu0) at y = 2 mm above a 12 x 0.5 x 12 mm copper plate
(sigma = 5.8e7 S/m), translating x = -6 -> +4 mm while rotating about z
(2 turns).  Radia computes the open-boundary external field A_s, B_ext
ANALYTICALLY; the conductor eddy reaction is a reduced-potential A-phi FEM on
the plate (HCurl nograds A_r + H1 phi, current continuity).  See
radia_mcp.maglev topic radia_iem_fem.

Three ways to get the eddy current J, Joule loss P and Lorentz force F:

  1. kinematic (source-only):  J = -sigma * d/dt A_s   -- NO eddy FEM at all.
     Exact in the low magnetic-Reynolds-number limit Rm = mu0 sigma omega L^2
     << 1, where the eddy reaction barely perturbs the field.
  2. full-FEM A-phi:           the reference (carries the reaction A_r, phi).
  3. CLN-reduced A-phi:        a constant-basis multiport Cauer Ladder Network
     (Tanimoto/Sugahara/Takahashi/Matsuo, IEEE TMag 2023, generalized to a
     2-parameter translating+rotating source by SVD-seeded block-Krylov over
     the transient iteration matrix A_sys^-1 M).  Per-step cost ~ M x M.

The crossover (this script's headline): the kinematic source-only error grows
with Rm.  Below Rm ~ 0.1 it is < 0.5 % -- the analytic source formula suffices
and NO per-step FEM is needed (Yano's actual case, Rm ~ 0.016, error 0.035 %).
Above Rm ~ 1 the reaction matters (error %-level), and the CLN reproduces the
full-FEM to ~1e-4..1e-6 at ~1000x less per-step cost -- THIS is where the CLN
earns its keep (e.g. faster motion, thicker / more conductive rails, kHz drive,
or the TEAM 28 aluminium disk).

Run:
    python rotating_magnet_eddy.py            # Rm sweep + Yano case -> json + png
    python rotating_magnet_eddy.py --maxh 0.001 --ncfg 24   # quick / coarse

Outputs durable JSON under validation_test/maglev/demos and keeps the public
figure next to this script.
"""
from __future__ import annotations
import argparse
import json
import os
import time

import numpy as np
import radia as rad
from ngsolve import (Mesh, HCurl, H1, VectorH1, BilinearForm, LinearForm,
                     GridFunction, curl, grad, InnerProduct, dx, Integrate,
                     TaskManager)
from netgen.occ import Box, OCCGeometry

from _validation_output import validation_output

MU0 = 4 * np.pi * 1e-7
SIGMA = 5.8e7              # copper conductivity [S/m]
BR = 0.2                  # magnet remanence [T]
MAG_MM = 1.0              # cube edge [mm]
GAP_Y = 0.002             # magnet centre height [m]
X0, X1 = -0.006, 0.004    # translation range [m]
NSTEP = 180               # nominal trajectory resolution
ROT_PER_STEP = 4.0        # deg/step (2 turns over NSTEP)
L_CHAR = 0.006            # plate half-size for Rm [m]
REG_MASS = 1e-4           # rel. mass regularization of the A-phi near-kernel


def _pose(step):
    xp = X0 + (X1 - X0) * step / NSTEP
    th = np.deg2rad(step * ROT_PER_STEP)
    c, sn = np.cos(th), np.sin(th)
    return dict(origin=[xp, GAP_Y, 0.0], u_axis=[c, sn, 0.0],
                v_axis=[-sn, c, 0.0], w_axis=[0.0, 0.0, 1.0])


def project_source(fesS, steps):
    """Project the analytic Radia A_s and B_ext onto VectorH1 once per config."""
    mag = rad.magnet_box([0, 0, 0], [MAG_MM * 1e-3] * 3, [0, BR / MU0, 0])
    As = np.zeros((fesS.ndof, len(steps)))
    Be = np.zeros((fesS.ndof, len(steps)))
    gA = GridFunction(fesS)
    gB = GridFunction(fesS)
    for j, st in enumerate(steps):
        p = _pose(int(st))
        gA.Set(rad.RadiaField(mag, 'a', **p)); As[:, j] = gA.vec.FV().NumPy().copy()
        gB.Set(rad.RadiaField(mag, 'b', **p)); Be[:, j] = gB.vec.FV().NumPy().copy()
    rad.UtiDelAll()
    return As, Be


def build(maxh):
    plate = Box((-0.006, -0.0005, -0.006), (0.006, 0.0, 0.006)).mat("copper")
    mesh = Mesh(OCCGeometry(plate).GenerateMesh(maxh=maxh))
    fesA = HCurl(mesh, order=2, nograds=True, dirichlet=".*")
    fesP = H1(mesh, order=2, dirichlet=".*")
    fes = fesA * fesP
    fesS = VectorH1(mesh, order=2)
    (A, phi), (W, psi) = fes.TnT()
    # sigma-mass of the COMBINED field (A + grad phi).  This M is near-singular
    # (its kernel is A + grad phi = 0); a tiny full-mass regularization
    # (REG_MASS rel.) lifts that near-kernel so the same M / A_sys = M + dt K is
    # robustly factorable on ANY mesh + dt (used consistently in BOTH the system
    # matrix and the RHS history term -> an O(REG_MASS) well-posed perturbation,
    # no dt-dependent bias).  ~0.01% shift, well under the example's accuracy.
    msig = SIGMA * (InnerProduct(A, W) + InnerProduct(grad(phi), W)
                    + InnerProduct(A, grad(psi)) + InnerProduct(grad(phi), grad(psi))
                    + REG_MASS * (InnerProduct(A, W) + InnerProduct(grad(phi), grad(psi))))
    Mbf = BilinearForm(fes, symmetric=True); Mbf += msig * dx; Mbf.Assemble()
    Kterm = (1.0 / MU0) * InnerProduct(curl(A), curl(W))
    return mesh, fes, fesA, fesP, fesS, (A, phi, W, psi), Mbf, msig, Kterm


def run(maxh=0.0005, ncfg=48, T_list=(2.0, 0.2, 0.02, 2e-3), verbose=True):
    t0 = time.time()
    with TaskManager():
        mesh, fes, fesA, fesP, fesS, tnt, Mbf, msig, Kterm = build(maxh)
        A, phi, W, psi = tnt
        ndof = fes.ndof
    Vc = Integrate(1.0, mesh)
    steps = np.linspace(0, NSTEP - 1, ncfg).round().astype(int)
    if verbose:
        print(f"plate ndof={ndof} (A {fesA.ndof}+phi {fesP.ndof}) Vc={Vc:.3e} m^3"
              f"  projecting source... [{time.time()-t0:.1f}s]", flush=True)
    As, Bext = project_source(fesS, steps)

    # source loads g(config) = int sigma A_s . (W, grad psi)
    gAs = GridFunction(fesS); G = np.zeros((ndof, ncfg))
    for j in range(ncfg):
        gAs.vec.FV().NumPy()[:] = As[:, j]
        lf = LinearForm(fes); lf += SIGMA * (InnerProduct(gAs, W) + InnerProduct(gAs, grad(psi))) * dx
        lf.Assemble(); G[:, j] = lf.vec.FV().NumPy().copy()

    _t = GridFunction(fes); _r = GridFunction(fes)

    def Mv(x):
        _t.vec.FV().NumPy()[:] = x; _r.vec.data = Mbf.mat * _t.vec
        return _r.vec.FV().NumPy().copy()

    gfu = GridFunction(fes); gfo = GridFunction(fes)
    gAn = GridFunction(fesS); gAo = GridFunction(fesS); gBn = GridFunction(fesS)

    def outputs(u, u_o, n, dt, src_only=False):
        """eddy J_rms, Joule P, Lorentz |F| at step n."""
        gAn.vec.FV().NumPy()[:] = As[:, n]; gAo.vec.FV().NumPy()[:] = As[:, n - 1]
        gBn.vec.FV().NumPy()[:] = Bext[:, n]
        if src_only:
            d = (gAn - gAo)
        else:
            gfu.vec.FV().NumPy()[:] = u; gfo.vec.FV().NumPy()[:] = u_o
            Ar, Ph = gfu.components; Ar_o, Ph_o = gfo.components
            d = (gAn - gAo) + (Ar - Ar_o) + (grad(Ph) - grad(Ph_o))
        J = (-SIGMA / dt) * d
        Jr = np.sqrt(Integrate(InnerProduct(J, J), mesh) / Vc)
        P = Integrate(InnerProduct(J, J) / SIGMA, mesh)
        fx = Integrate(J[1] * gBn[2] - J[2] * gBn[1], mesh, order=6)
        fy = Integrate(J[2] * gBn[0] - J[0] * gBn[2], mesh, order=6)
        fz = Integrate(J[0] * gBn[1] - J[1] * gBn[0], mesh, order=6)
        return Jr, P, float(np.sqrt(fx * fx + fy * fy + fz * fz))

    rows = []
    yano = None
    for T in T_list:
        dt = T / ncfg
        Rm = MU0 * SIGMA * (2 * np.pi * 2 / T) * L_CHAR ** 2
        with TaskManager():
            Abf = BilinearForm(fes, symmetric=True); Abf += (msig + dt * Kterm) * dx; Abf.Assemble()
            Ainv = Abf.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

        def Av(b):
            _t.vec.FV().NumPy()[:] = b; _r.vec.data = Ainv * _t.vec
            return _r.vec.FV().NumPy().copy()

        def Asv(x):
            _t.vec.FV().NumPy()[:] = x; _r.vec.data = Abf.mat * _t.vec
            return _r.vec.FV().NumPy().copy()

        # --- full-FEM reference + kinematic (source-only) ---
        u = np.zeros(ndof); U = np.zeros((ndof, ncfg))
        Jf = np.zeros(ncfg); Pf = np.zeros(ncfg); Ff = np.zeros(ncfg)
        Js = np.zeros(ncfg); Ps = np.zeros(ncfg); Fs = np.zeros(ncfg)
        tfull = 0.0
        for n in range(1, ncfg):
            rhs = Mv(u) - (G[:, n] - G[:, n - 1])
            t1 = time.perf_counter(); u = Av(rhs); tfull += time.perf_counter() - t1
            U[:, n] = u
            Jf[n], Pf[n], Ff[n] = outputs(U[:, n], U[:, n - 1], n, dt)
            Js[n], Ps[n], Fs[n] = outputs(None, None, n, dt, src_only=True)

        # --- CLN-reduced (M=20) ---
        X0_ = np.column_stack([Av(G[:, j]) for j in range(ncfg)])
        X1_ = np.column_stack([Av(Mv(X0_[:, j])) for j in range(ncfg)])
        X2_ = np.column_stack([Av(Mv(X1_[:, j])) for j in range(ncfg)])
        pod = lambda Aa, k: np.linalg.svd(Aa, full_matrices=False)[0][:, :min(k, Aa.shape[1])]
        Vb, _ = np.linalg.qr(np.hstack([pod(X0_, 14), pod(X1_, 8), pod(X2_, 6)]))
        M = min(20, Vb.shape[1]); V = Vb[:, :M]
        Mr = V.T @ np.column_stack([Mv(V[:, k]) for k in range(M)])
        Ar_ = V.T @ np.column_stack([Asv(V[:, k]) for k in range(M)]); Arin = np.linalg.inv(Ar_)
        Gr = V.T @ G
        c = np.zeros(M); Jc = np.zeros(ncfg); Pc = np.zeros(ncfg); Fc = np.zeros(ncfg); rx = np.zeros(ncfg)
        u_o = np.zeros(ndof); tred = 0.0
        for n in range(1, ncfg):
            rhsr = (Mr @ c) - (Gr[:, n] - Gr[:, n - 1])
            t1 = time.perf_counter()
            for _ in range(50):
                cnew = Arin @ rhsr
            tred += (time.perf_counter() - t1) / 50
            ur = V @ cnew
            rx[n] = np.linalg.norm(ur - U[:, n]) / (np.linalg.norm(U[:, n]) + 1e-300)
            Jc[n], Pc[n], Fc[n] = outputs(ur, u_o, n, dt); u_o = ur; c = cnew

        def relmax(a, b):
            return float(np.abs(a[1:] - b[1:]).max() / (np.abs(b[1:]).max() + 1e-300))

        ms_full = tfull / (ncfg - 1) * 1e3
        ms_red = tred / (ncfg - 1) * 1e3
        row = {
            "T_s": T, "Rm": float(Rm), "Jrms_peak": float(Jf[1:].max()),
            "F_peak_N": float(Ff[1:].max()),
            "source_only_Jrms_err": relmax(Js, Jf), "source_only_F_err": relmax(Fs, Ff),
            "cln_Jrms_err": relmax(Jc, Jf), "cln_F_err": relmax(Fc, Ff),
            "cln_reaction_err": float(rx[1:].max()),
            "full_ms_per_step": ms_full, "cln_ms_per_step": ms_red,
            "speedup": ms_full / ms_red if ms_red > 0 else None,
        }
        rows.append(row)
        if abs(T - 2.0) < 1e-9:
            yano = {"Jrms_t": Jf.tolist(), "F_t": Ff.tolist(),
                    "Jrms_src_t": Js.tolist(), "F_src_t": Fs.tolist(), "dt_s": dt}
        if verbose:
            print(f"  Rm={Rm:8.3f} | J_rms peak={Jf[1:].max():9.1f} | source-only J err={relmax(Js,Jf):8.2e}"
                  f" F err={relmax(Fs,Ff):8.2e} | CLN J err={relmax(Jc,Jf):8.2e}"
                  f" | reaction captured={1-rx[1:].max():6.1%} | speedup {row['speedup']:.0f}x"
                  f" [{time.time()-t0:.0f}s]", flush=True)
    return {"ndof": int(ndof), "ncfg": int(ncfg), "maxh": maxh, "Vc": float(Vc),
            "rows": rows, "yano": yano}


def make_figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = res["rows"]
    Rm = [r["Rm"] for r in rows]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax1.loglog(Rm, [r["source_only_F_err"] for r in rows], 'o-', label="source-only (no FEM)")
    ax1.loglog(Rm, [r["cln_F_err"] for r in rows], 's-', label="CLN (M=20)")
    ax1.axhline(0.01, ls=":", color="0.5"); ax1.axvline(1.0, ls=":", color="0.5")
    ax1.set_xlabel("magnetic Reynolds number  Rm"); ax1.set_ylabel("Lorentz force rel. error")
    ax1.legend(fontsize=8, loc="upper left"); ax1.grid(True, which="both", alpha=0.3)
    y = res["yano"]
    tt = np.arange(len(y["F_t"])) * y["dt_s"]
    ax2.plot(tt, np.array(y["F_t"]) * 1e9, '-', label="full-FEM A-phi")
    ax2.plot(tt, np.array(y["F_src_t"]) * 1e9, '--', label="source-only")
    ax2.set_xlabel("time (s)"); ax2.set_ylabel("|F| (nN)")
    ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)
    fig.tight_layout(pad=0.4)
    fig.savefig(path, dpi=150); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maxh", type=float, default=0.0005)
    ap.add_argument("--ncfg", type=int, default=48)
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    res = run(maxh=args.maxh, ncfg=args.ncfg)
    output = validation_output("rotating_magnet_eddy_results.json", here)
    with open(output, "w") as f:
        json.dump(res, f, indent=2, allow_nan=False)
    if not args.no_figure:
        try:
            make_figure(res, os.path.join(here, "rotating_magnet_eddy.png"))
        except Exception as e:                      # matplotlib optional
            print(f"(figure skipped: {e!r})")
    print("\n=== Rm crossover ===")
    print("  Rm < ~0.1  : source-only J = -sigma dA_s/dt suffices (no per-step FEM) -- Yano's case")
    print("  Rm > ~1    : eddy reaction matters; CLN(M=20) reproduces full-FEM to ~1e-4 at ~1000x")
    print(f"  saved {output}")


if __name__ == "__main__":
    main()
