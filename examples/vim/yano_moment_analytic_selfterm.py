"""Analytic self-term for the parameter-free yano-MSC MOMENT formulation -- removes BOTH the eval-point
alpha AND the finite-difference conditioning noise.

Background (see yano_moment_formulation.py): the moment formulation closes the 6 face charges of each hex by
their field MOMENTS about the centroid -- parameter-free, no eval-point alpha:

    dipole(3)     : <M> = dipole(sigma)/V = chi * H(centroid)
    quadrupole(2) : quad(sigma) (mixing-corrected) = chi * gradH(centroid)
    monopole(1)   : sum_f A_f sigma_f = 0           (neutrality == discrete div(B)=0)
    loops         : NOT deflated -- the moment matrix A is non-singular, so the (field-null) cell-graph
                    cycle modes are solved, not projected out.  Deflating them to zero (an earlier prototype
                    did) breaks the per-element constitutive (the loops carry per-element dipole); that is
                    invisible in the external moment but FATAL in nonlinear -- see yano_moment_nonlinear.py.

That prototype obtained H(centroid) and gradH(centroid) by finite-differencing the yano operator over two
small alphas (N0 = N(a->0), N1 = dN/da).  The differencing injects noise that drives the conditioning up on
sheared meshes (cond ~ 3e5 at 6^3 shear).  This file replaces it with the CLOSED-FORM centroid field and
gradient:

    H(r)     = (1/4pi) \\int_face sigma (r-r')/|r-r'|^3 dA'
    gradH(r) = (1/4pi) \\int_face sigma [ I/|r-r'|^3 - 3 (r-r')(r-r')^T/|r-r'|^5 ] dA'   (symmetric, traceless)

evaluated at the element centroid (interior -> smooth integrand -> exact by Gauss quadrature here; closed-form
Wilton/Graglia face integrals in the C++ kernel).  Two physics points the validation pins down:

  * SELF term uses the BARE face charge (no center charge).  The single-cube patch test (exact demag N=1/3)
    is reproduced to 0.000% this way, and the centroid is interior so it is finite -- no singularity, no
    center point needed.

  * MUTUAL term uses the yano DIPOLE LAYER = (face charge) - area*(point charge @ source-element center).
    With the deflation removed, the bare variant (no mutual center charge, the 'bare' column) is ALSO
    accurate -- it matches the analytic; the per-face center charge gives only a modest conditioning
    improvement (e.g. 6^3 shear: analytic 6e2 vs bare 2e3).  (An earlier deflated prototype reported the
    center charge as REQUIRED -- bare cond 2e6, -20% at 4^3 regular -- but that breakage was a DEFLATION
    artifact, not intrinsic.)  The source-element center is at finite distance from the self centroid, so
    the mutual center charge is finite (it is excluded exactly where it would be singular: self).

RESULT: the analytic self-term matches the finite-difference moment accuracy (~0.1-0.6%, vs EIEM2's 1.5-4.5%)
while dropping the conditioning by ~300x on sheared meshes (6^3 shear: analytic 6e2 vs finite-diff 2e5 -- the
finite-diff alpha-extrapolation noise is real and survives without deflation).  This is the parameter-free,
well-conditioned closed-form that the C++ GetCentroidFieldGrad kernel implements.
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
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (0, 4, 7, 3)]


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
    rad.UtiDelAll(); rad.set_demag_backend("auto"); rad.SolverConfig(yano_eval_alpha=-1.0, yano_no_center_charge=False)
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


# ----- analytic charged-quad-face field/gradient kernel -----------------------------------------------

def _gauss01(ng):
    x, w = np.polynomial.legendre.leggauss(ng)
    return 0.5 * (x + 1.0), 0.5 * w


def _face_quadrature(verts, gp, gw):
    V0, V1, V2, V3 = [np.asarray(v, float) for v in verts]
    u, v = np.meshgrid(gp, gp, indexing="ij"); wu, wv = np.meshgrid(gw, gw, indexing="ij")
    u = u.ravel(); v = v.ravel(); w = (wu * wv).ravel()
    P = (np.outer((1-u)*(1-v), V0) + np.outer(u*(1-v), V1) + np.outer(u*v, V2) + np.outer((1-u)*v, V3))
    dPdu = np.outer(1-v, V1-V0) + np.outer(v, V2-V3)
    dPdv = np.outer(1-u, V3-V0) + np.outer(u, V2-V1)
    dA = np.linalg.norm(np.cross(dPdu, dPdv), axis=1) * w
    return P, dA


def recon_face_verts(cells, G):
    """4 verts per DOF, aligned to the C++ DOF order (match face center)."""
    elem = G[:, 0].astype(int); fc = G[:, 2:5]; out = []
    for f in range(G.shape[0]):
        C = np.asarray(cells[elem[f]], float); tgt = fc[f]; best = None; bd = 1e18
        for q in HEX_FACES:
            quad = C[list(q)]; dd = np.sum((quad.mean(0) - tgt) ** 2)
            if dd < bd:
                bd = dd; best = quad
        out.append(best)
    return out


def precompute_quadrature(verts_list, ng):
    gp, gw = _gauss01(ng)
    P, dA = [], []
    for V in verts_list:
        p, a = _face_quadrature(V, gp, gw); P.append(p); dA.append(a)
    return np.array(P), np.array(dA)


def field_grad_at(obs, allP, alldA):
    """bare face-charge H (dof,3), gradH (dof,3,3) at obs."""
    d = obs[None, None, :] - allP
    r2 = np.einsum("fnk,fnk->fn", d, d)
    w3 = alldA * r2 ** -1.5
    w5 = alldA * r2 ** -2.5
    H = np.einsum("fnk,fn->fk", d, w3) / (4 * np.pi)
    goff = -3.0 * np.einsum("fni,fnj,fn->fij", d, d, w5)
    gH = (goff + np.eye(3)[None] * w3.sum(1)[:, None, None]) / (4 * np.pi)
    return H, gH


def field_grad_dipole_layer(obs, allP, alldA, src_elem, src_area, src_ecenter, self_e):
    """dipole-layer H, gradH: bare face MINUS area*(point @ source center).  SELF center excluded (would sit
    at obs); MUTUAL center is finite-distance -> singularity-free."""
    H, gH = field_grad_at(obs, allP, alldA)
    mutual = (src_elem != self_e)
    d = obs[None, :] - src_ecenter
    r2 = np.where(mutual, np.einsum("fk,fk->f", d, d), 1.0)
    inv_r3 = np.where(mutual, r2 ** -1.5, 0.0); inv_r5 = np.where(mutual, r2 ** -2.5, 0.0)
    q = -src_area
    H = H + (q * inv_r3)[:, None] * d / (4 * np.pi)
    gc = q[:, None, None] * (np.eye(3)[None] * inv_r3[:, None, None]
                             - 3.0 * np.einsum("fi,fj->fij", d, d) * inv_r5[:, None, None]) / (4 * np.pi)
    return H, gH + gc


# ----- moment assembly (shared row logic; field source = finite-diff OR analytic) ---------------------

def _norm(row, rhs):
    nn = np.linalg.norm(row)
    return (row / nn, rhs / nn) if nn > 1e-300 else (row, rhs)


def _assemble(G, dof, get_F0_Ginv):
    elem = G[:, 0].astype(int); area = G[:, 1]; fc = G[:, 2:5]; nrm = G[:, 5:8]; ecen = G[:, 8:11]
    n_el = int(elem.max()) + 1; dofs_of = [np.where(elem == e)[0] for e in range(n_el)]
    A = np.zeros((dof, dof)); b = np.zeros(dof); r = 0
    for e in range(n_el):
        fs = dofs_of[e]; Ae = area[fs]; ne = nrm[fs]; ce = ecen[fs[0]]; d = fc[fs] - ce
        Ve = (1.0 / 3.0) * np.sum(Ae * np.sum(fc[fs] * ne, axis=1))
        F0, Ginv = get_F0_Ginv(e, fs, ne, d, ce)
        dip = np.zeros((3, dof))
        for k in range(3):
            dip[k, fs] = Ae * d[:, k]
        for k in range(3):                                # dipole: <M> = chi H(centroid)
            row, rhs = _norm(dip[k, :] / Ve - CHI * F0[k, :], CHI * (H0 if k == 2 else 0.0))
            A[r, :] = row; b[r] = rhs; r += 1
        mono = np.zeros(dof); mono[fs] = Ae               # monopole: neutrality == div(B)=0
        row, rhs = _norm(mono, 0.0); A[r, :] = row; b[r] = rhs; r += 1
        for Bm in (d[:, 0]**2 - d[:, 1]**2, d[:, 1]**2 - d[:, 2]**2):   # quadrupole (mixing-corrected)
            row = np.zeros(dof); row[fs] = Ae * Bm
            cm = np.array([np.sum(Ae * ne[:, k] * Bm) for k in range(3)])
            row -= (cm @ dip) / Ve
            Dm = np.array([[np.sum(Ae * d[:, jj] * ne[:, ii] * Bm) for jj in range(3)] for ii in range(3)])
            Dvec = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], Dm[0, 1]+Dm[1, 0], Dm[0, 2]+Dm[2, 0], Dm[1, 2]+Dm[2, 1]])
            row -= CHI * (Dvec @ Ginv)
            row, rhs = _norm(row, 0.0); A[r, :] = row; b[r] = rhs; r += 1
    return A, b


def _solve(A, b, Lb, nLoop, G):
    # The moment matrix A is NON-SINGULAR -> solve directly.  Deflating the (field-null) loop modes to zero
    # (an earlier prototype did) would break the per-element constitutive -- the loops carry per-element
    # dipole -- which is invisible in the external moment but FATAL in nonlinear (see yano_moment_nonlinear.py).
    # Reporting cond(A) of the FULL system also corrects the earlier conditioning claim, which was measured on
    # the deflated P^T A P.  (Lb, nLoop kept in the signature for call-site compatibility; unused here.)
    area = G[:, 1]; fc = G[:, 2:5]; ecen = G[:, 8:11]
    sig = np.linalg.solve(A, b)
    m = float(np.sum(sig * area * (fc[:, 2] - ecen[:, 2])))
    return m, float(np.linalg.cond(A))


def eiem2_err(cells, m_mmm):
    N, G, Lb, nLoop, dof = matgeom(cells, 0.5)
    area = G[:, 1]; cen = G[:, 2:5]; ecen = G[:, 8:11]; nrm = G[:, 5:8]
    A = (1.0 / CHI) * np.eye(dof) - N                              # non-singular -> no deflation needed
    sig = np.linalg.solve(A, H0 * nrm[:, 2])
    return float(np.sum(sig * area * (cen[:, 2] - ecen[:, 2]))) / m_mmm - 1


def moment_finite_diff(cells, m_mmm):
    a1, a2 = 0.10, 0.30
    Na, G, Lb, nLoop, dof = matgeom(cells, a1)
    Nb, _, _, _, _ = matgeom(cells, a2)
    N0 = (a2 * Na - a1 * Nb) / (a2 - a1); N1 = (Nb - Na) / (a2 - a1)

    def get(e, fs, ne, d, ce):
        F0 = np.linalg.pinv(ne) @ N0[fs, :]
        Mproj = np.array([[ne[fi, 0]*d[fi, 0], ne[fi, 1]*d[fi, 1], ne[fi, 2]*d[fi, 2],
                           ne[fi, 0]*d[fi, 1]+ne[fi, 1]*d[fi, 0], ne[fi, 0]*d[fi, 2]+ne[fi, 2]*d[fi, 0],
                           ne[fi, 1]*d[fi, 2]+ne[fi, 2]*d[fi, 1]] for fi in range(6)])
        Ginv = np.linalg.pinv(Mproj) @ N1[fs, :]
        return F0, Ginv

    A, b = _assemble(G, dof, get)
    m, c = _solve(A, b, Lb, nLoop, G)
    return m / m_mmm - 1, c


def moment_analytic(cells, m_mmm, ng=8, dipole_layer=True):
    _, G, Lb, nLoop, dof = matgeom(cells, 0.5)
    elem = G[:, 0].astype(int); area = G[:, 1]; ecen = G[:, 8:11]
    allP, alldA = precompute_quadrature(recon_face_verts(cells, G), ng)

    def get(e, fs, ne, d, ce):
        if dipole_layer:
            H, gH = field_grad_dipole_layer(ce, allP, alldA, elem, area, ecen, e)
        else:
            H, gH = field_grad_at(ce, allP, alldA)
        return H.T, np.array([gH[:, 0, 0], gH[:, 1, 1], gH[:, 2, 2], gH[:, 0, 1], gH[:, 0, 2], gH[:, 1, 2]])

    A, b = _assemble(G, dof, get)
    m, c = _solve(A, b, Lb, nLoop, G)
    return m / m_mmm - 1, c


def patch_test():
    """single cube: exact volume-average demag N=1/3 -> M_z = chi H0/(1+chi/3); analytic SELF must hit it."""
    cells = build_hex(1, 0.0)
    _, G, Lb, nLoop, dof = matgeom(cells, 0.5)
    allP, alldA = precompute_quadrature(recon_face_verts(cells, G), 10)

    def get(e, fs, ne, d, ce):
        H, gH = field_grad_at(ce, allP, alldA)
        return H.T, np.array([gH[:, 0, 0], gH[:, 1, 1], gH[:, 2, 2], gH[:, 0, 1], gH[:, 0, 2], gH[:, 1, 2]])

    A, b = _assemble(G, dof, get)
    area = G[:, 1]; fc = G[:, 2:5]; ecen = G[:, 8:11]
    sig = np.linalg.solve(A, b)                                    # non-singular -> no deflation
    Mz = float(np.sum(sig * area * (fc[:, 2] - ecen[:, 2]))) / (2 * L) ** 3
    return Mz, CHI * H0 / (1.0 + CHI / 3.0)


def main():
    Mz, Mz_exact = patch_test()
    print(f"\nPATCH TEST (single cube, exact N=1/3): M_z={Mz:.3f} vs {Mz_exact:.3f}  err={Mz/Mz_exact-1:+.3%}\n")
    print("  accuracy (vs MMM) and CONDITIONING -- EIEM2(tuned) vs moment(finite-diff) vs moment(ANALYTIC)")
    print("  ANALYTIC = closed-form centroid field+grad (SELF=bare face, MUTUAL=dipole layer, mono=div(B)=0);")
    print("  'bare' = no mutual center charge (now accurate too -- a modest conditioning aid, not required).\n")
    print(f"  {'mesh':>11} | {'EIEM2':>8} | {'moment FD':>10} {'cond':>7} | {'ANALYTIC':>10} {'cond':>7} | {'bare':>9} {'cond':>7}")
    print("  " + "-" * 86)
    rows = []
    for (n, s) in [(3, 0.0), (3, 0.6), (4, 0.0), (4, 0.6), (6, 0.0), (6, 0.6)]:
        cells = build_hex(n, s); m_mmm = mmm_moment(build_tet(n, s))
        e_e = eiem2_err(cells, m_mmm)
        e_fd, c_fd = moment_finite_diff(cells, m_mmm)
        e_an, c_an = moment_analytic(cells, m_mmm, dipole_layer=True)
        e_ba, c_ba = moment_analytic(cells, m_mmm, dipole_layer=False)
        tag = f"{n}^3 {'reg' if s == 0 else 'shear'}"
        print(f"  {tag:>11} | {e_e:>+7.2%} | {e_fd:>+9.2%} {c_fd:>7.0e} | {e_an:>+9.2%} {c_an:>7.0e} | {e_ba:>+8.2%} {c_ba:>7.0e}")
        rows.append(dict(n=n, s=s, eiem2_err=e_e, fd_err=e_fd, fd_cond=c_fd,
                         an_err=e_an, an_cond=c_an, bare_err=e_ba, bare_cond=c_ba))

    out = dict(mu_r=MU_R, patch_Mz=Mz, patch_Mz_exact=Mz_exact, rows=rows,
               conclusion=("Analytic closed-form centroid field+gradient (SELF=bare face charge, patch-test "
                           "exact; MUTUAL=yano dipole layer; monopole=div(B)=0) solved with NO loop deflation "
                           "(A is non-singular). Reproduces the finite-difference moment accuracy (~0.1-0.6%, "
                           "vs EIEM2 1.5-4.5%) AND resolves the finite-difference alpha-extrapolation "
                           "conditioning noise: ~300x lower cond on sheared meshes (6^3 shear analytic ~6e2 vs "
                           "FD ~2e5) -- a REAL benefit that survives without deflation. The mutual center charge "
                           "is a modest conditioning aid, NOT required (the 'bare' column now matches analytic; "
                           "the earlier 'bare breaks regular grids / cond 2e6' was a DEFLATION artifact). This is "
                           "the kernel for C++ GetCentroidFieldGrad; deflation removed per yano_moment_nonlinear.py."))
    with open(os.path.join(HERE, "yano_moment_analytic_selfterm.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\n  Analytic self-term: parameter-free AND well-conditioned. saved",
          os.path.join(HERE, "yano_moment_analytic_selfterm.json"))


if __name__ == "__main__":
    main()
