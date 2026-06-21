"""Consideration: what does the element-center cancellation charge DO, and can div(B)=0-via-Lagrange
replace it?  (user, 2026-06-21: "中心磁荷を外した場合も考察しよう".)

The yano-MSC kernel adds a point charge -FaceArea at the element center per source face, so each element
is net-neutral by construction.  rad.SolverConfig(yano_no_center_charge=True) drops it (RAW collocation).

Measured on a sheared hex cube (mu_r=200), three formulations vs an independent MMM (tet) reference:
  * CENTER (Yano)            : the historical kernel.
  * RAW                      : no center charge, solved as-is.
  * RAW + div(B)=0 + loops   : no center charge, but deflate the loop space AND the monopole (net-charge)
                               space -- i.e. impose div(B)=0 by Lagrange/deflation instead of the charge.

FINDING:
  * The center charge is STRUCTURALLY ESSENTIAL, not a near-0 nicety.  Removing it ENLARGES the field-null
    from the cycle space (cycles) to the FULL internal-face space (null = n_internal_faces = cycles +
    (n_cells-1)); the RAW solve then has ~63% per-element net charge and a >100% moment error (unusable).
  * The (n_cells-1) extra null modes are monopole(net-charge)-type: deflating [loops + div(B)=0] removes
    them, restores bounded conditioning, and recovers the SAME accuracy as the center-charge kernel
    (~1% vs MMM).  So div(B)=0-by-Lagrange CAN replace the center charge -- it must handle BOTH the loops
    AND the monopole modes (div=0 alone or loops alone is not enough), and it costs the loop + monopole
    deflation that the single center charge does implicitly for free.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src", "radia"))
import radia as rad  # noqa: E402

MU0 = 4e-7 * math.pi
MU_R = 200.0; H0 = 1000.0; L = 0.020; CHI = MU_R - 1.0
R1, R2 = 0.40, 0.80
FREUD = [(0, 1, 2, 6), (0, 1, 5, 6), (0, 3, 2, 6), (0, 3, 7, 6), (0, 4, 5, 6), (0, 4, 7, 6)]


def phi(P, s):
    x, y, z = P; zr = z / L
    return [x + s * L * (zr + 0.5 * zr * zr), y + 0.4 * s * L * (zr - 0.3 * zr * zr), z]


def box_corners(ax, i, j, k):
    return [[ax[i], ax[j], ax[k]], [ax[i+1], ax[j], ax[k]], [ax[i+1], ax[j+1], ax[k]], [ax[i], ax[j+1], ax[k]],
            [ax[i], ax[j], ax[k+1]], [ax[i+1], ax[j], ax[k+1]], [ax[i+1], ax[j+1], ax[k+1]], [ax[i], ax[j+1], ax[k+1]]]


def build_hex(n, s):
    ax = np.linspace(-L, L, n + 1)
    return [[phi(P, s) for P in box_corners(ax, i, j, k)] for k in range(n) for j in range(n) for i in range(n)]


def build_tet(n, s):
    ax = np.linspace(-L, L, n + 1); cells = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = [phi(P, s) for P in box_corners(ax, i, j, k)]
                cells += [[c[a], c[b], c[cc], c[d]] for (a, b, cc, d) in FREUD]
    return cells


def external_moment(cont):
    b1 = rad.Fld(cont, "b", [0, 0, R1])[2] - MU0 * H0
    b2 = rad.Fld(cont, "b", [0, 0, R2])[2] - MU0 * H0
    A = np.array([[1 / R1**3, 1 / R1**5], [1 / R2**3, 1 / R2**5]])
    a, _ = np.linalg.solve(A, [b1, b2]); return 2 * math.pi * a / MU0


def mmm_moment(cells):
    rad.UtiDelAll(); rad.set_demag_backend("auto"); rad.SolverConfig(yano_pyramid_cloud=False, yano_no_center_charge=False)
    objs = [rad.ObjTetrahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for t in objs:
        rad.MatApl(t, rad.MatLin(MU_R))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-10, 8000, 0); m = external_moment(cont); rad.UtiDelAll(); return m


def matgeom(cells, no_center):
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(yano_pyramid_cloud=False, yano_no_center_charge=no_center)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MU_R))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    N, dof = rad.GetInteractMatrix(handle); G = rad.GetFaceGeom(handle); Lb, nLoop = rad.GetLoopBasis(handle)
    rad.SolverConfig(yano_no_center_charge=False); rad.UtiDelAll()
    return np.asarray(N, float), np.asarray(G, float), np.asarray(Lb, float), nLoop, dof


def comp(W):
    Wn = W / (np.linalg.norm(W, axis=0, keepdims=True) + 1e-300)
    U, S, _ = np.linalg.svd(Wn, full_matrices=True)
    return U[:, int(np.sum(S > 1e-9 * S[0])):]


def metrics(N, G, Lb, nLoop, n_el, m_mmm, P=None):
    elem = G[:, 0].astype(int); area = G[:, 1]; cen = G[:, 2:5]; ecen = G[:, 8:11]; nrm = G[:, 5:8]
    dof = N.shape[0]
    A = (1.0 / CHI) * np.eye(dof) - N
    b = H0 * nrm[:, 2]
    sig = np.linalg.solve(A, b) if P is None else P @ np.linalg.solve(P.T @ A @ P, P.T @ b)
    m = float(np.sum(sig * area * (cen[:, 2] - ecen[:, 2])))
    net = np.zeros(n_el); den = np.zeros(n_el)
    for d in range(dof):
        net[elem[d]] += area[d] * sig[d]; den[elem[d]] += area[d] * abs(sig[d])
    sv = np.linalg.svd(N, compute_uv=False); nulld = int(np.sum(sv < 1e-9 * sv[0]))
    Ad = A if P is None else P.T @ A @ P
    return dict(null=nulld, net=float((np.abs(net) / (den + 1e-300)).max()),
                err=m / m_mmm - 1, cond=float(np.linalg.cond(Ad)))


def main():
    print(f"\nCenter cancellation charge: what it does, and can div(B)=0 replace it?  mu_r={MU_R:.0f}\n")
    out = []
    for n in (4, 6):
        s = 0.6
        cells = build_hex(n, s); n_el = len(cells)
        m_mmm = mmm_moment(build_tet(n, s))
        Nc, G, Lb, nLoop, dof = matgeom(cells, False)
        Nr, _, _, _, _ = matgeom(cells, True)
        area = G[:, 1]; elem = G[:, 0].astype(int)
        C = np.zeros((n_el, dof))
        for d in range(dof):
            C[elem[d], d] = area[d]
        P_lm = comp(np.hstack([Lb, C.T]))
        rc = metrics(Nc, G, Lb, nLoop, n_el, m_mmm)
        rr = metrics(Nr, G, Lb, nLoop, n_el, m_mmm)
        rd = metrics(Nr, G, Lb, nLoop, n_el, m_mmm, P=P_lm)
        print(f"=== n={n} s={s}: dof={dof}, n_cells={n_el}, cycles={nLoop}, n_internal_faces={nLoop + n_el - 1} ===")
        print(f"  {'':<26}{'null(N)':>9}{'net charge':>12}{'err vs MMM':>12}{'cond':>11}")
        print(f"  {'CENTER (Yano)':<26}{rc['null']:>9d}{rc['net']:>12.1e}{rc['err']:>+11.2%}{rc['cond']:>11.1e}")
        print(f"  {'RAW (no center charge)':<26}{rr['null']:>9d}{rr['net']:>12.1e}{rr['err']:>+11.1%}{rr['cond']:>11.1e}")
        print(f"  {'RAW + div(B)=0 + loops':<26}{'-':>9}{rd['net']:>12.1e}{rd['err']:>+11.2%}{rd['cond']:>11.1e}")
        print(f"  (RAW null {rr['null']} == cycles {nLoop} + (n_cells-1) {n_el-1};  MMM moment {m_mmm:.4e})\n")
        out.append(dict(n=n, n_cells=n_el, cycles=nLoop, n_internal_faces=nLoop + n_el - 1,
                        center=rc, raw=rr, raw_div0_loops=rd, mmm_moment=m_mmm))
    with open(os.path.join(HERE, "yano_center_charge_study.json"), "w") as f:
        json.dump({"mu_r": MU_R, "cases": out,
                   "conclusion": ("The center cancellation charge is STRUCTURALLY ESSENTIAL: removing it "
                                  "enlarges the field-null from the cycle space to the full internal-face "
                                  "space (null = n_internal_faces = cycles + (n_cells-1)); the RAW solve has "
                                  "~63% net charge and >100% moment error. The (n_cells-1) extra null modes "
                                  "are monopole-type: RAW + div(B)=0 + loop deflation removes them, restores "
                                  "bounded conditioning AND recovers ~1% accuracy vs MMM. So div(B)=0-by-"
                                  "Lagrange CAN replace the center charge, but must deflate BOTH loops and "
                                  "monopoles -- what the single center charge does implicitly for free.")},
                  f, indent=2, default=float)
    print("saved", os.path.join(HERE, "yano_center_charge_study.json"))


if __name__ == "__main__":
    main()
