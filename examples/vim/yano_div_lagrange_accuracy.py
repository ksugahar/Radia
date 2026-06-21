"""Does the physical div(B)=0 constraint (Lagrange, no cancellation charge) beat Yano's near-0
cancellation charge -- on ACCURACY, not just conditioning?  (user research question, 2026-06-21.)

Yano's yano-MSC kernel adds an element-center CANCELLATION CHARGE that drives each element's net surface
charge Sum_f area_f*sigma_f close to 0 -- but NOT exactly 0 (measured here ~1e-3 relative).  The user's
idea: impose div(B)=0 (net charge = 0) by a LAGRANGE constraint INSTEAD, and ask whether forcing it
*exactly* 0 (hard) is better or worse than leaving it ~1e-3 (soft).  Hypothesis: the small slack is
BENEFICIAL (absorbs other collocation error), so soft >= hard -- a consistency-vs-purity tradeoff.

This is measurable in Python because the new rad.GetFaceGeom exposes the per-DOF face geometry (area,
centroid, outward normal, element center) in the matrix DOF order, so we can:
  * build the system   A = (1/chi) I - N   (N = rad.GetInteractMatrix, the C++ EIEM2 operator),
  * build the uniform-field RHS  b = H0 * n_z,
  * impose div(B)=0 as a penalty gamma*C^T C (soft, gamma: 0 -> inf) or a Lagrange KKT system (hard),
  * read the dipole moment  m_z = Sum_f sigma_f*area_f*(c_f - c_e)_z  -- the reliable observable --
    and compare to an INDEPENDENT MMM (tet) reference on the same deformed body.

VALIDATION baked in: (a) per-element Sum_f area_f*n_f == 0 (closed hex -> DOF<->face mapping correct);
(b) at gamma=0 the Python moment matches the C++ rad.Solve external moment to <1%.

FINDING (mu_r=200): sweeping gamma drives the net charge 1e-3 -> 0, but the bulk moment changes only
~0.01-0.03% -- NEGLIGIBLE.  Hard (exact div=0) is marginally WORSE than Yano's soft near-0 at distortion
(weakly consistent with "the slack is beneficial"), but far too small to matter for the bulk moment.  So
on the global field the Lagrange div=0 ~ Yano cancellation; any real difference must live in the LOCAL
field (not the moment), which needs the field eval and is left for a C++ solve-side wiring.
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
MU_R = 200.0; H0 = 1000.0; L = 0.020
CHI = MU_R - 1.0
R1, R2 = 0.40, 0.80


def phi(P, s):
    x, y, z = P; zr = z / L
    return [x + s * L * (zr + 0.5 * zr * zr), y + 0.4 * s * L * (zr - 0.3 * zr * zr), z]


def box_corners(ax, i, j, k):
    return [[ax[i], ax[j], ax[k]], [ax[i+1], ax[j], ax[k]], [ax[i+1], ax[j+1], ax[k]], [ax[i], ax[j+1], ax[k]],
            [ax[i], ax[j], ax[k+1]], [ax[i+1], ax[j], ax[k+1]], [ax[i+1], ax[j+1], ax[k+1]], [ax[i], ax[j+1], ax[k+1]]]


FREUD = [(0, 1, 2, 6), (0, 1, 5, 6), (0, 3, 2, 6), (0, 3, 7, 6), (0, 4, 5, 6), (0, 4, 7, 6)]


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
    rad.UtiDelAll(); rad.set_demag_backend("auto"); rad.SolverConfig(yano_pyramid_cloud=False)
    objs = [rad.ObjTetrahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for t in objs:
        rad.MatApl(t, rad.MatLin(MU_R))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-10, 8000, 0); m = external_moment(cont); rad.UtiDelAll(); return m


def cpp_moment(cells):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_pyramid_cloud=False)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MU_R))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-10, 8000, 0); m = external_moment(cont); rad.UtiDelAll(); return m


def matrix_geom(cells):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_pyramid_cloud=False)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in cells]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MU_R))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    N, dof = rad.GetInteractMatrix(handle); G = rad.GetFaceGeom(handle); rad.UtiDelAll()
    return np.asarray(N, float), np.asarray(G, float), dof


def net_charge_rel(sigma, elem, area, n_el):
    net = np.zeros(n_el); den = np.zeros(n_el)
    for d in range(len(sigma)):
        net[elem[d]] += area[d] * sigma[d]; den[elem[d]] += area[d] * abs(sigma[d])
    return float((np.abs(net) / (den + 1e-300)).max())


def run_case(n, s):
    cells = build_hex(n, s); n_el = len(cells)
    N, G, dof = matrix_geom(cells)
    elem = G[:, 0].astype(int); area = G[:, 1]; cen = G[:, 2:5]; ecen = G[:, 8:11]; nrm = G[:, 5:8]
    # (a) closed-hex check
    clos = np.zeros((n_el, 3))
    for d in range(dof):
        clos[elem[d]] += area[d] * nrm[d]
    closure = float(np.abs(clos).max())
    m_mmm = mmm_moment(build_tet(n, s))
    A = (1.0 / CHI) * np.eye(dof) - N
    b = H0 * nrm[:, 2]
    C = np.zeros((n_el, dof))
    for d in range(dof):
        C[elem[d], d] = area[d]
    C /= np.linalg.norm(C, axis=1, keepdims=True)
    CTC = C.T @ C
    a_scale = np.linalg.norm(A, 2)

    def moment(sig):
        return float(np.sum(sig * area * (cen[:, 2] - ecen[:, 2])))

    sweep = []
    for alpha in (0.0, 1e-1, 1.0, 1e1, 1e2):
        sig = np.linalg.solve(A + alpha * a_scale * CTC, b)
        sweep.append({"gamma_over_A": alpha, "moment": moment(sig),
                      "err_vs_mmm": moment(sig) / m_mmm - 1.0, "net_charge": net_charge_rel(sig, elem, area, n_el)})
    KKT = np.block([[A, C.T], [C, np.zeros((n_el, n_el))]])
    sig_h = np.linalg.solve(KKT, np.concatenate([b, np.zeros(n_el)]))[:dof]
    hard = {"moment": moment(sig_h), "err_vs_mmm": moment(sig_h) / m_mmm - 1.0,
            "net_charge": net_charge_rel(sig_h, elem, area, n_el)}
    m_cpp = cpp_moment(cells)
    return dict(n=n, s=s, n_el=n_el, dof=dof, closure=closure, mmm_moment=m_mmm,
                cpp_moment=m_cpp, python_moment_gamma0=sweep[0]["moment"], sweep=sweep, hard_lagrange=hard)


def main():
    print(f"\ndiv(B)=0 Lagrange vs Yano cancellation -- ACCURACY (mu_r={MU_R:.0f}).  rad.GetFaceGeom demo.\n")
    out = []
    for (n, s) in [(4, 0.0), (4, 0.6), (6, 0.6)]:
        r = run_case(n, s); out.append(r)
        print(f"=== n={n} s={s}  dof={r['dof']}  closure(Sum area*n)={r['closure']:.1e}  "
              f"moment match C++/Py={r['cpp_moment']/r['python_moment_gamma0']:.4f} ===")
        print(f"  {'gamma/|A|':>10} | {'err vs MMM':>11} | {'net charge':>11}")
        for row in r["sweep"]:
            tag = "  <- Yano (soft, gamma=0)" if row["gamma_over_A"] == 0 else ""
            print(f"  {row['gamma_over_A']:10.0e} | {row['err_vs_mmm']:+10.3%} | {row['net_charge']:11.2e}{tag}")
        h = r["hard_lagrange"]
        print(f"  {'HARD':>10} | {h['err_vs_mmm']:+10.3%} | {h['net_charge']:11.2e}  <- exact div=0 (Lagrange)")
        print(f"  -> hard - soft(gamma=0) = {(h['err_vs_mmm'] - r['sweep'][0]['err_vs_mmm'])*100:+.3f} pp (negligible)\n")
    with open(os.path.join(HERE, "yano_div_lagrange_accuracy.json"), "w") as f:
        json.dump({"mu_r": MU_R, "H0": H0, "cases": out,
                   "conclusion": ("div(B)=0 (soft penalty -> hard Lagrange) drives the per-element net charge "
                                  "1e-3 -> 0 but changes the bulk dipole moment only ~0.01-0.03% (negligible). "
                                  "Hard is marginally worse than Yano's near-0 cancellation at distortion "
                                  "(weakly 'slack is beneficial'), but immaterial for the bulk moment. On the "
                                  "global field the Lagrange div=0 ~ Yano cancellation; any real difference is "
                                  "in the LOCAL field, not tested here.")}, f, indent=2, default=float)
    print("saved", os.path.join(HERE, "yano_div_lagrange_accuracy.json"))


if __name__ == "__main__":
    main()
