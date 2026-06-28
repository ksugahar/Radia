"""bench_multipole_moment_scaling.py -- scale benchmark for the parameter-free multipole-moment MMM NONLINEAR solver.

Runs the validated (no-deflation, Anderson-accelerated secant) moment formulation on the voxelized hex
C-yoke at increasing resolution and records the Benchmark-Policy fields (peak_memory_mb, t_setup, t_solve,
iterations, converged) + correctness self-checks (div(B)=0, per-element constitutive residual) + the dense
operator footprint dof^2 (the O(N^2) wall that motivates the HACApK phase).

mdx-ready: imports `radia` DIRECTLY (PyPI on mdx / editable on LAB) -- NO src-path hack that would shadow the
PyPI install -- and is fully self-contained (no sibling-module import), so

    python bench_multipole_moment_scaling.py --sizes 20 30 40 50

runs on mdx the moment a radia release carrying rad.GetCentroidFieldGrad / rad.GetFaceGeom lands.  The solver
here is byte-identical to the validated examples/vim/multipole_moment_nonlinear.py moment_nonlinear (split into a
SETUP phase = objs + C++ accessors + per-element local geometry, and a SOLVE phase = the Anderson iteration,
so the two costs are timed separately).  Kept inline so the benchmark has zero import-path dependency.

THE WALL this benchmark measures: the moment system A is assembled DENSE (dof x dof) and re-solved each
Anderson iteration, so memory is O(dof^2) and the per-iteration solve is O(dof^3).  That caps the dense path
at a few thousand dof -- exactly the regime the next phase (HACApK ACA-compression of the centroid
field+gradient operator F0/Ginv) is meant to break.  This script is the dense BASELINE the HACApK path must
reproduce (correctness) and beat (memory + time) at scale.
"""
import argparse
import json
import math
import os
import platform
import time
from datetime import datetime

import numpy as np
import psutil

import radia as rad   # PyPI (mdx) or editable (LAB); NO sys.path hack -> mdx-clean

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * math.pi
MU_R0 = 1000.0; CHI0 = MU_R0 - 1.0; MSAT = 1.0e6


def chi_secant(Hmag):
    return CHI0 / (1.0 + CHI0 * Hmag / MSAT)


def M_BH_vec(H):
    Hm = np.linalg.norm(H, axis=1, keepdims=True)
    Mm = CHI0 * Hm / (1.0 + CHI0 * Hm / MSAT)
    return np.where(Hm > 1e-30, Mm * H / Hm, 0.0)


def _inside_cyoke(cx, cy):
    # outer square [-0.06,0.06]^2, minus inner cavity [-0.035,0.035]^2, minus gap x>0.018 (C opens +x)
    return (-0.06 <= cx <= 0.06) and (-0.06 <= cy <= 0.06) and not (-0.035 < cx < 0.035 and -0.035 < cy < 0.035) and not (cx > 0.018)


def build_cyoke_hexes(nxy, nz):
    xs = np.linspace(-0.06, 0.06, nxy + 1); zs = np.linspace(-0.02, 0.02, nz + 1)
    hexes = []
    for k in range(nz):
        for j in range(nxy):
            for i in range(nxy):
                if not _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
                    continue
                x0, x1, y0, y1, z0, z1 = xs[i], xs[i + 1], xs[j], xs[j + 1], zs[k], zs[k + 1]
                hexes.append([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                              [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    return hexes


def _norm(row, rhs):
    nn = np.linalg.norm(row)
    return (row / nn, rhs / nn) if nn > 1e-300 else (row, rhs)


def _peak_mb():
    mi = psutil.Process(os.getpid()).memory_info()
    return (mi.peak_wset if hasattr(mi, "peak_wset") else mi.rss) / (1024.0 * 1024.0)


def moment_setup(hexes):
    """SETUP phase: objs + C++ accessors (rad.GetFaceGeom / rad.GetCentroidFieldGrad) + per-element local
    geometry (dipole/quadrupole functionals).  Returns the assembled-once per-element data the SOLVE reuses."""
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in hexes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MU_R0))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    G = np.asarray(rad.GetFaceGeom(handle), float); C = np.asarray(rad.GetCentroidFieldGrad(handle), float)
    dof = G.shape[0]; elem = G[:, 0].astype(int); area = G[:, 1]; fc = G[:, 2:5]; nrm = G[:, 5:8]; ecen = G[:, 8:11]
    n_el = int(elem.max()) + 1; dofs_of = [np.where(elem == e)[0] for e in range(n_el)]
    EL = []
    for e in range(n_el):
        fs = dofs_of[e]; Ae = area[fs]; ne = nrm[fs]; cc = ecen[fs[0]]; d = fc[fs] - cc
        Ve = (1 / 3) * np.sum(Ae * np.sum(fc[fs] * ne, axis=1)); F0 = C[e, 0:3, :]; Ginv = C[e, 3:9, :]
        dip = np.zeros((3, dof))
        for k in range(3):
            dip[k, fs] = Ae * d[:, k]
        quads = []
        for Bm in (d[:, 0]**2 - d[:, 1]**2, d[:, 1]**2 - d[:, 2]**2):
            qf = np.zeros(dof); qf[fs] = Ae * Bm
            cm = np.array([np.sum(Ae * ne[:, kk] * Bm) for kk in range(3)]); qf -= (cm @ dip) / Ve
            Dm = np.array([[np.sum(Ae * d[:, jj] * ne[:, ii] * Bm) for jj in range(3)] for ii in range(3)])
            Dvec = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], Dm[0, 1]+Dm[1, 0], Dm[0, 2]+Dm[2, 0], Dm[1, 2]+Dm[2, 1]])
            quads.append((qf, Dvec @ Ginv))
        EL.append((fs, F0, dip, Ve, quads))
    return dict(dof=dof, n_el=n_el, area=area, fc=fc, ecen=ecen, dofs_of=dofs_of, EL=EL)


def moment_solve(S, Happ, maxit=80, tol=1e-10, m_depth=6):
    """SOLVE phase: Anderson-accelerated secant fixed-point on the DENSE moment system, NO loop deflation
    (A is non-singular -- solve directly).  Same map as multipole_moment_nonlinear.py."""
    dof = S["dof"]; n_el = S["n_el"]; EL = S["EL"]; area = S["area"]

    def Gmap(sig):
        Hs = np.array([Happ + EL[e][1] @ sig for e in range(n_el)])
        chi = chi_secant(np.linalg.norm(Hs, axis=1))
        A = np.zeros((dof, dof)); b = np.zeros(dof); r = 0
        for e in range(n_el):
            ce = chi[e]; fs, F0, dip, Ve, quads = EL[e]
            for k in range(3):
                row, rhs = _norm(dip[k, :] / Ve - ce * F0[k, :], ce * Happ[k]); A[r, :] = row; b[r] = rhs; r += 1
            mono = np.zeros(dof); mono[fs] = area[fs]
            row, rhs = _norm(mono, 0.0); A[r, :] = row; b[r] = rhs; r += 1
            for (qf, gradrow) in quads:
                row, rhs = _norm(qf - ce * gradrow, 0.0); A[r, :] = row; b[r] = rhs; r += 1
        return np.linalg.solve(A, b)

    sig = np.zeros(dof); Sh = []; Fh = []; nit = 0; conv = False
    for it in range(maxit):
        nit = it + 1
        g = Gmap(sig); f = g - sig
        if np.linalg.norm(f) / (np.linalg.norm(g) + 1e-30) < tol:
            sig = g; conv = True; break
        Sh.append(g); Fh.append(f)
        if len(Fh) > m_depth:
            Sh.pop(0); Fh.pop(0)
        if len(Fh) == 1:
            sig = g
        else:
            dF = np.array([Fh[i+1] - Fh[i] for i in range(len(Fh)-1)]).T
            dG = np.array([Sh[i+1] - Sh[i] for i in range(len(Sh)-1)]).T
            gamma, *_ = np.linalg.lstsq(dF, Fh[-1], rcond=None)
            sig = Sh[-1] - dG @ gamma
    return sig, nit, conv


def run_case(nxy, nz, H0, maxit):
    hexes = build_cyoke_hexes(nxy, nz)
    t0 = time.perf_counter()
    S = moment_setup(hexes)
    t_setup = time.perf_counter() - t0
    Happ = np.array([0.0, H0, 0.0])
    t1 = time.perf_counter()
    sig, nit, conv = moment_solve(S, Happ, maxit=maxit)
    t_solve = time.perf_counter() - t1
    EL = S["EL"]; n_el = S["n_el"]; area = S["area"]; fc = S["fc"]; ecen = S["ecen"]; dofs_of = S["dofs_of"]; dof = S["dof"]
    Hs = np.array([Happ + EL[e][1] @ sig for e in range(n_el)])
    divB = max(abs(np.sum(area[dofs_of[e]] * sig[dofs_of[e]])) / (np.sum(area[dofs_of[e]] * np.abs(sig[dofs_of[e]])) + 1e-30) for e in range(n_el))
    constit = float(np.median([np.linalg.norm(EL[e][2] @ sig / EL[e][3] - M_BH_vec(Hs[e:e+1])[0]) / (np.linalg.norm(EL[e][2] @ sig / EL[e][3]) + 1e-30) for e in range(n_el)]))
    my = float(np.sum(sig * area * (fc[:, 1] - ecen[:, 1])))
    rad.UtiDelAll()
    return dict(nxy=nxy, nz=nz, nhex=int(n_el), dof=int(dof), H0=float(H0), moment_my=my,
                iterations=int(nit), converged=bool(conv),
                t_setup=float(t_setup), t_solve=float(t_solve), peak_memory_mb=float(_peak_mb()),
                dense_matrix_mb=float(dof * dof * 8 / 1.0e6), divB_resid=float(divB), constit_resid=constit)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[8, 12, 16, 20],
                    help="C-yoke in-plane resolutions nxy to sweep (mdx: push higher, e.g. 20 30 40 50)")
    ap.add_argument("--nz", type=int, default=2, help="through-thickness layers")
    ap.add_argument("--field", type=float, default=5.0e4, help="applied H_y in A/m (default 5e4 = BH knee)")
    ap.add_argument("--maxit", type=int, default=80, help="max Anderson iterations")
    ap.add_argument("--out", default=os.path.join(HERE, "bench_multipole_moment_scaling.json"),
                    help="output JSON path (committed next to the script per the Data Persistence Policy)")
    args = ap.parse_args()

    print(f"\nmultipole-moment MMM NONLINEAR scale benchmark (dense baseline): chi0={CHI0:.0f}, Msat={MSAT:.0e}, "
          f"H_app={args.field:.0e} A/m\n")
    print(f"  {'mesh':>10} {'nhex':>6} {'dof':>6} | {'iters':>5} {'conv':>5} | {'t_setup':>8} {'t_solve':>8} "
          f"{'peak_MB':>8} {'denseMB':>8} | {'div(B)':>8} {'constit':>8}")
    print("  " + "-" * 104)
    results = []
    for nxy in args.sizes:
        rc = run_case(nxy, args.nz, args.field, args.maxit)
        results.append(rc)
        print(f"  {nxy}x{nxy}x{args.nz:<2} {rc['nhex']:>6} {rc['dof']:>6} | {rc['iterations']:>5} "
              f"{str(rc['converged']):>5} | {rc['t_setup']:>7.2f}s {rc['t_solve']:>7.2f}s "
              f"{rc['peak_memory_mb']:>7.0f} {rc['dense_matrix_mb']:>7.1f} | {rc['divB_resid']:>8.1e} "
              f"{rc['constit_resid']:>8.1e}")

    data = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": "multipole_moment_scaling",
        "problem": dict(geometry="voxelized hex C-yoke", mu_r0=MU_R0, Msat=MSAT, H_applied=args.field,
                        nz=args.nz, sizes=args.sizes, solver="Anderson-accelerated secant, dense, no deflation"),
        "results": results,
        "note": ("Dense baseline for the parameter-free multipole-moment MMM nonlinear solver. Memory is O(dof^2) "
                 "(dense_matrix_mb) and the per-iteration solve is O(dof^3), so this path caps at a few "
                 "thousand dof -- the HACApK phase (ACA-compression of the centroid field+gradient operator "
                 "F0/Ginv) must reproduce these moment_my / div(B) / constit numbers and beat peak_memory_mb "
                 "+ t_solve at scale. mdx-ready: imports radia directly (no src-path hack)."),
    }
    with open(args.out, "w") as f:
        json.dump(data, f, indent=2, default=float)
    print(f"\n  saved {args.out}")
    print("  dense baseline -- HACApK phase must reproduce moment_my/div(B)/constit and beat peak_MB/t_solve.")


if __name__ == "__main__":
    main()
