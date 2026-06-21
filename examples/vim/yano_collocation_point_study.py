"""The collocation-point problem: yano-MSC accuracy depends on WHERE we evaluate.

A collocation method enforces the constitutive relation at discrete EVALUATION POINTS.  yano-MSC uses,
per face, EvalPt = alpha*FaceCenter + (1-alpha)*ElementCenter, with alpha an UN-DERIVED choice
(EIEM2 = 0.5, the pyramid kernel = 0.75).  This sweeps alpha (rad.SolverConfig(yano_eval_alpha=...),
center charge on, pyramid off) and measures accuracy vs an independent MMM reference on a sheared cube.

FINDING: accuracy varies STRONGLY with alpha -- a ~4% swing from alpha=0.3 (-0.6% vs MMM) to alpha=0.95
(-4.7%) -- with a clean ACCURACY <-> CONDITIONING trade-off:
  * eval near the element CENTER (small alpha): the field there is the more representative element-
    average-like value -> MORE ACCURATE, but it is neighbour-influenced / less self-dominant ->
    WORSE conditioned (and alpha->0 makes all 6 face eval points coincide -> degenerate).
  * eval near the FACE (large alpha): the field is the self-dominant near-face value -> WELL conditioned,
    but it is the near-singular surface-charge field -> LESS representative -> LESS ACCURATE.
So EIEM2 (0.5) and the pyramid point (0.75) are heuristic spots on this curve, NEITHER optimal; EIEM2 is
the more accurate of the two (which is why the pyramid kernel measured "worse" on the bulk moment).

CONSEQUENCE: the collocation point is a real, un-principled accuracy knob -- the ~1% accuracy of EIEM2 is
a point on a trade-off, not a robust property.  Loop deflation bounds the conditioning for ALL alpha (the
cond column is mu_r-independent, only its constant moves), so alpha can be tuned for accuracy; but there is
no principled choice.  The principled fix that REMOVES the eval-point arbitrariness is INTEGRATION
(Galerkin) instead of point collocation -- exactly what the HDiv-VIM analytic Gram does (no eval point).
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
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    rad.SolverConfig(yano_pyramid_cloud=False, yano_no_center_charge=False, yano_eval_alpha=-1.0)
    objs = [rad.ObjTetrahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for t in objs:
        rad.MatApl(t, rad.MatLin(MU_R))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-10, 8000, 0); m = external_moment(cont); rad.UtiDelAll(); return m


def matgeom(cells, alpha):
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(yano_pyramid_cloud=False, yano_no_center_charge=False, yano_eval_alpha=alpha)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MU_R))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    N, dof = rad.GetInteractMatrix(handle); G = rad.GetFaceGeom(handle); Lb, nLoop = rad.GetLoopBasis(handle)
    rad.SolverConfig(yano_eval_alpha=-1.0); rad.UtiDelAll()
    return np.asarray(N, float), np.asarray(G, float), np.asarray(Lb, float), nLoop, dof


def main():
    n, s = 6, 0.6
    cells = build_hex(n, s); m_mmm = mmm_moment(build_tet(n, s))
    print(f"\nyano-MSC collocation-point sweep, sheared cube n={n} s={s}, center charge ON, vs MMM {m_mmm:.4e}\n")
    print(f"  {'alpha':>6} | {'moment':>11} | {'err vs MMM':>11} | {'cond (loop-defl)':>16}")
    print("  " + "-" * 56)
    rows = []
    for alpha in (0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.85, 0.95):
        N, G, Lb, nLoop, dof = matgeom(cells, alpha)
        area = G[:, 1]; cen = G[:, 2:5]; ecen = G[:, 8:11]; nrm = G[:, 5:8]
        A = (1.0 / CHI) * np.eye(dof) - N
        sig = np.linalg.solve(A, H0 * nrm[:, 2])
        m = float(np.sum(sig * area * (cen[:, 2] - ecen[:, 2])))
        U, _S, _ = np.linalg.svd(Lb, full_matrices=True); P = U[:, nLoop:]
        cond = float(np.linalg.cond(P.T @ A @ P))
        tag = "  <- EIEM2" if alpha == 0.5 else ("  <- pyramid pt" if alpha == 0.75 else "")
        print(f"  {alpha:6.2f} | {m:11.5e} | {m/m_mmm-1:+10.2%} | {cond:16.1f}{tag}")
        rows.append(dict(alpha=alpha, moment=m, err_vs_mmm=m / m_mmm - 1, cond_loop_deflated=cond))
    with open(os.path.join(HERE, "yano_collocation_point_study.json"), "w") as f:
        json.dump({"n": n, "s": s, "mmm_moment": m_mmm, "sweep": rows,
                   "conclusion": ("The yano-MSC collocation point (EvalPt = alpha*FaceCenter + "
                                  "(1-alpha)*center) is an un-derived accuracy knob: err vs MMM swings ~4% "
                                  "from alpha=0.3 (-0.6%) to alpha=0.95 (-4.7%), a clean accuracy<->conditioning "
                                  "trade-off (small alpha=center: accurate, ill-conditioned; large alpha=face: "
                                  "well-conditioned, inaccurate). EIEM2 (0.5) and pyramid (0.75) are heuristic "
                                  "points, neither optimal. Loop deflation bounds cond for all alpha so alpha can "
                                  "be tuned for accuracy, but there is no principled choice; the principled fix "
                                  "is INTEGRATION (Galerkin / analytic Gram, as in HDiv-VIM) -- no eval point.")},
                  f, indent=2, default=float)
    print("\n  ~4% accuracy swing with alpha -> the collocation-point choice IS a real accuracy problem.")
    print("  No principled alpha; loop deflation bounds cond for all of them. Principled fix = integration")
    print("  (Galerkin / analytic Gram, as in HDiv-VIM), which has no eval point.")
    print("saved", os.path.join(HERE, "yano_collocation_point_study.json"))


if __name__ == "__main__":
    main()
