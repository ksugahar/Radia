"""Parameter-free MOMENT formulation of the yano-MSC -- the principled fix for the collocation-point
(eval-point alpha) problem.  Designed from physics (the multipole structure of the 6 face charges):

  6 face charges  =  monopole(1) + dipole(3) + quadrupole(2),

and the constitutive relation M(r) = chi H(r), expanded about the element CENTROID, fixes each block by
its own field MOMENT -- no eval-point alpha:

  * dipole (3):     <M> = dipole(sigma)/V = chi * H(centroid)        [0th moment: field at centroid]
  * quadrupole (2): grad(M) = chi * grad(H(centroid))                [1st moment: field gradient]
  * monopole (1):   sum_f area_f sigma_f = 0                         [neutrality / div(B)=0]
  * loops:          deflated (the cell-graph cycle basis, field-null)

Evaluating at the centroid is the parameter-free choice (the centroid is where the point field best
represents the volume average -- patch-test exact for the dipole); the 6 face eval-points only coincide
there (degenerate) for POINT collocation, but the MOMENTS (field + gradient) are independent conditions,
so there is no degeneracy and no alpha to tune.

KEY DERIVATION (what makes the quadrupole correct on DISTORTED hexes): the quadrupole moment of sigma
MIXES a dipole(M_0) part and a gradient(grad M) part.  On a symmetric (cube) element the M_0 part
vanishes, but on a sheared hex it does not -- so the quadrupole condition must SUBTRACT the M_0 mixing
(coefficient c_m = sum_f area_f n_f B_m) before equating the gradient part to chi*gradH.

Historical implementation note: the original prototype estimated the centroid field and gradient by
finite-differencing the now-retired EIEM2 eval-point operator. Current Radia no longer exposes those
research knobs; the production moment system is built directly by `rad.BuildMomentSystem` /
`BuildMomentSystemCore`.

RESULT (vs an independent MMM reference): the parameter-free moment formulation is near-exact on regular
meshes (~0.05%) and ~0.3% on strongly sheared meshes -- 5x to 50x better than the (alpha-tuned) EIEM2
collocation, with NO eval-point parameter.  (The conditioning on sheared meshes is higher than EIEM2 and
is the remaining refinement; the accuracy is the point.)
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

MU0 = 4e-7 * math.pi
MU_R = 200.0; H0 = 1000.0; L = 0.020; CHI = MU_R - 1.0
R1, R2 = 0.40, 0.80
FREUD = [(0, 1, 2, 6), (0, 1, 5, 6), (0, 3, 2, 6), (0, 3, 7, 6), (0, 4, 5, 6), (0, 4, 7, 6)]


def _radia():
    """Lazy import so the archive-printing path works without a local built _radia_pybind.pyd."""
    if "radia" not in sys.modules:
        sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
    import radia as rad
    return rad


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
    rad = _radia()
    b1 = rad.Fld(cont, "b", [0, 0, R1])[2] - MU0 * H0
    b2 = rad.Fld(cont, "b", [0, 0, R2])[2] - MU0 * H0
    A = np.array([[1 / R1**3, 1 / R1**5], [1 / R2**3, 1 / R2**5]])
    a, _ = np.linalg.solve(A, [b1, b2]); return 2 * math.pi * a / MU0


def mmm_moment(cells):
    rad = _radia()
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    objs = [rad.ObjTetrahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for t in objs:
        rad.MatApl(t, rad.MatLin(MU_R))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-10, 8000, 0); m = external_moment(cont); rad.UtiDelAll(); return m


def matgeom(cells, alpha):
    raise RuntimeError(
        "This historical EIEM2 alpha-sweep prototype cannot be recomputed on current Radia. "
        "Use yano_moment_formulation.json for the archived numbers, and use rad.BuildMomentSystem "
        "or GetCentroidFieldGrad for the canonical moment-yano path."
    )


def eiem2_err(cells, m_mmm):
    """the alpha=0.5 EIEM2 collocation, loop-deflated (the tuned baseline)."""
    N, G, Lb, nLoop, dof = matgeom(cells, 0.5)
    area = G[:, 1]; cen = G[:, 2:5]; ecen = G[:, 8:11]; nrm = G[:, 5:8]
    U, _S, _ = np.linalg.svd(Lb, full_matrices=True); P = U[:, nLoop:]
    A = (1.0 / CHI) * np.eye(dof) - N
    sig = P @ np.linalg.solve(P.T @ A @ P, P.T @ (H0 * nrm[:, 2]))
    return float(np.sum(sig * area * (cen[:, 2] - ecen[:, 2]))) / m_mmm - 1


def moment_err(cells, m_mmm):
    """parameter-free moment formulation: dipole<-H(centroid), quadrupole<-gradH, monopole<-neutrality."""
    a1, a2 = 0.10, 0.30
    Na, G, Lb, nLoop, dof = matgeom(cells, a1)
    Nb, _, _, _, _ = matgeom(cells, a2)
    N0 = (a2 * Na - a1 * Nb) / (a2 - a1)          # field at centroid (alpha -> 0)
    N1 = (Nb - Na) / (a2 - a1)                    # field gradient (dN/dalpha)
    elem = G[:, 0].astype(int); area = G[:, 1]; fc = G[:, 2:5]; nrm = G[:, 5:8]; ecen = G[:, 8:11]
    n_el = int(elem.max()) + 1
    dofs_of = [np.where(elem == e)[0] for e in range(n_el)]
    A = np.zeros((dof, dof)); b = np.zeros(dof); r = 0

    def norm(row, rhs):
        nn = np.linalg.norm(row)
        return (row / nn, rhs / nn) if nn > 1e-300 else (row, rhs)

    for e in range(n_el):
        fs = dofs_of[e]; Ae = area[fs]; ne = nrm[fs]; ce = ecen[fs[0]]
        d = fc[fs] - ce
        Ve = (1.0 / 3.0) * np.sum(Ae * np.sum(fc[fs] * ne, axis=1))
        F0 = np.linalg.pinv(ne) @ N0[fs, :]                              # field vector at centroid
        Mproj = np.array([[ne[fi, 0]*d[fi, 0], ne[fi, 1]*d[fi, 1], ne[fi, 2]*d[fi, 2],
                           ne[fi, 0]*d[fi, 1] + ne[fi, 1]*d[fi, 0], ne[fi, 0]*d[fi, 2] + ne[fi, 2]*d[fi, 0],
                           ne[fi, 1]*d[fi, 2] + ne[fi, 2]*d[fi, 1]] for fi in range(6)])
        Ginv = np.linalg.pinv(Mproj) @ N1[fs, :]                          # gradH (sym 6) at centroid
        dip = np.zeros((3, dof))
        for k in range(3):
            dip[k, fs] = Ae * d[:, k]
        for k in range(3):                                               # dipole: <M> = chi H(centroid)
            row, rhs = norm(dip[k, :] / Ve - CHI * F0[k, :], CHI * (H0 if k == 2 else 0.0))
            A[r, :] = row; b[r] = rhs; r += 1
        mono = np.zeros(dof); mono[fs] = Ae                             # monopole: neutrality
        row, rhs = norm(mono, 0.0); A[r, :] = row; b[r] = rhs; r += 1
        for Bm in (d[:, 0]**2 - d[:, 1]**2, d[:, 1]**2 - d[:, 2]**2):     # quadrupole (mixing-corrected)
            row = np.zeros(dof); row[fs] = Ae * Bm
            cm = np.array([np.sum(Ae * ne[:, k] * Bm) for k in range(3)])
            row -= (cm @ dip) / Ve                                       # subtract the dipole(M_0) mixing
            Dm = np.array([[np.sum(Ae * d[:, jj] * ne[:, ii] * Bm) for jj in range(3)] for ii in range(3)])
            Dvec = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], Dm[0, 1]+Dm[1, 0], Dm[0, 2]+Dm[2, 0], Dm[1, 2]+Dm[2, 1]])
            row -= CHI * (Dvec @ Ginv)                                   # = chi * gradient part
            row, rhs = norm(row, 0.0); A[r, :] = row; b[r] = rhs; r += 1
    U, _S, _ = np.linalg.svd(Lb, full_matrices=True); P = U[:, nLoop:]
    sig = P @ np.linalg.solve(P.T @ A @ P, P.T @ b)
    m = float(np.sum(sig * area * (fc[:, 2] - ecen[:, 2])))
    return m / m_mmm - 1, float(np.linalg.cond(P.T @ A @ P))


def main():
    json_path = os.path.join(HERE, "yano_moment_formulation.json")
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    print(f"\nArchived parameter-free MOMENT formulation results (mu_r={data['mu_r']:.0f}, vs MMM)\n")
    print(f"  {'mesh':>11} | {'EIEM2 (tuned a)':>15} | {'MOMENT (param-free)':>20}")
    print("  " + "-" * 54)
    for row in data["rows"]:
        tag = f"{row['n']}^3 {'reg' if row['s'] == 0 else 'shear'}"
        print(f"  {tag:>11} | {row['eiem2_err']:>+14.2%} | {row['moment_err']:>+13.2%} (cond {row['moment_cond']:.0e})")
    print("\n  Historical EIEM2 alpha-sweep knobs were removed in Phase 3.")
    print("  Current production recomputation goes through rad.BuildMomentSystem / BuildMomentSystemCore.")
    print("source", json_path)


if __name__ == "__main__":
    main()
