"""Golden: the C++ HIGH-ORDER charge-Gram (radia._radia_pybind._ChargeGramHMatrix high-order mode) -- the
SCALABLE production path for the order-p HDiv-type VIM demag operator (the order-0 analytic charge Gram
extended to polynomial charges).

The C++ Gram uses a MONOMIAL charge basis on each host element's reference coords; this test builds the
matching monomial charge-density map B_mono (a universal-per-element-type change-of-basis from the NGSolve
L2/SurfaceL2 density map -- the affine |J| cancels) + the per-charge data, hands the reference Gauss-Duffy
rules to the C++ constructor, and forms the demag with a single H-matvec:
    demag = <c, G c> / <m, M_mass m>,   c = B_mono @ m.
N = B_mono^T G B_mono is basis-invariant, so the demag is checked against the ANALYTIC gates (no dense
Python reference needed in the repo):
  * uniform M_z -> demag ~ 1/3, ORDER-INVARIANT p=0,1,2 (the polynomial charge map + monomial Gram are
    consistent; the high-order DOFs carry the constant surface charge without corrupting it);
  * non-uniform M=(0,0,z) -> demag p-CONVERGES (p=1 == p=2 because the surface charge M.n = z*n_z is linear
    and SurfaceL2(>=1) captures it exactly; p=0 distinct).  Exercises the volume Gram blocks (rho != 0).

NGSolve + Netgen required (importorskip); meshing is non-deterministic, so bands + the convergence
structure are checked (the values match the C:/temp dense reference to ~1e-4: 0.3333 / 0.374->0.410->0.410).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402
import radia._radia_pybind as rp  # noqa: E402


def _g01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def _tet_ref(o):
    s, ws = _g01(o)
    P, W = [], []
    for a, wa in zip(s, ws):
        for b, wb in zip(s, ws):
            for c, wc in zip(s, ws):
                P.append((a, b * (1 - a), c * (1 - a) * (1 - b)))
                W.append(wa * wb * wc * (1 - a) ** 2 * (1 - b))
    return np.array(P), np.array(W)


def _tri_ref(o):
    s, ws = _g01(o)
    P, W = [], []
    for u, wu in zip(s, ws):
        for v, wv in zip(s, ws):
            P.append((u, v * (1 - u)))
            W.append(wu * wv * (1 - u))
    return np.array(P), np.array(W)


def _monos_vol(pv):
    return [(i, j, k) for i in range(pv + 1) for j in range(pv + 1 - i) for k in range(pv + 1 - i - j)]


def _monos_surf(p):
    return [(i, j) for i in range(p + 1) for j in range(p + 1 - i)]


def _S(fe, mons, refP, refW, dim):
    """universal L2->monomial change-of-basis on the reference element: S = Mmono^{-1} C (|J| cancels)."""
    nm, nsh = len(mons), fe.ndof
    M = np.zeros((nm, nm))
    C = np.zeros((nm, nsh))
    for pt, w in zip(refP, refW):
        if dim == 3:
            l = pt
            mv = np.array([l[0] ** i * l[1] ** j * l[2] ** k for (i, j, k) in mons])
            sh = np.array(fe.CalcShape(l[0], l[1], l[2]))
        else:
            mv = np.array([pt[0] ** i * pt[1] ** j for (i, j) in mons])
            sh = np.array(fe.CalcShape(pt[0], pt[1]))
        M += w * np.outer(mv, mv)
        C += w * np.outer(mv, sh)
    return np.linalg.solve(M, C)


def _cube(h):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=h))


def _cpp_demag(mesh, p, Mcf, quad=6, eps=1e-7):
    """demag via the C++ high-order charge-Gram H-matrix (one matvec)."""
    import scipy.sparse as sp
    pv = max(p - 1, 0)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=p)
        nn = ng.specialcf.normal(mesh.dim)
        L2v, L2b = ng.L2(mesh, order=pv), ng.SurfaceL2(mesh, order=p)
        u = fes.TrialFunction()
        bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
        bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
        mv = ng.BilinearForm(L2v); mv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; mv.Assemble()
        mb = ng.BilinearForm(L2b); mb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; mb.Assemble()
        mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()

        def dense(B):
            r, c, v = B.mat.COO()
            return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(B.mat.height, B.mat.width)).toarray()
        Bv_d = np.linalg.solve(dense(mv), dense(bv))
        Bb_d = np.linalg.solve(dense(mb), dense(bb))
        Mmass = dense(mh)
        vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
        bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
        vV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in vels]
        bV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in bels]
        vdof = [list(L2v.GetDofNrs(e)) for e in vels]
        bdof = [list(L2b.GetDofNrs(e)) for e in bels]
        mons_v, mons_s = _monos_vol(pv), _monos_surf(p)
        Sv = _S(L2v.GetFE(vels[0]), mons_v, *_tet_ref(quad), dim=3)
        Ss = _S(L2b.GetFE(bels[0]), mons_s, *_tri_ref(quad), dim=2)
        gfu = ng.GridFunction(fes); gfu.Set(Mcf); m = np.array(gfu.vec)

    Brows, host, kind, expo = [], [], [], []
    for c in range(len(vels)):
        blk = Sv @ Bv_d[vdof[c], :]
        for a, (i, j, k) in enumerate(mons_v):
            Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    for f in range(len(bels)):
        blk = Ss @ Bb_d[bdof[f], :]
        for a, (i, j) in enumerate(mons_s):
            Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B_mono = np.array(Brows)

    rtp, rtw = _tet_ref(quad)
    rsp, rsw = _tri_ref(quad)
    G = rp._ChargeGramHMatrix(
        cell_verts=np.concatenate([V.ravel() for V in vV]).tolist(),
        face_verts=np.concatenate([V.ravel() for V in bV]).tolist(),
        n_el=len(vels), charge_host=host, charge_kind=kind, charge_expo=expo,
        ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
        ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(), eps=eps)
    c = B_mono @ m
    Gc = np.array(G.matvec(c.tolist()))
    return float(c @ Gc) / float(m @ Mmass @ m)


def test_cpp_ho_uniform_demag_order_invariant():
    """C++ high-order: uniform M_z on a cube -> demag ~ 1/3, IDENTICAL across p=0,1,2."""
    mesh = _cube(0.8)
    D = {p: _cpp_demag(mesh, p, ng.CoefficientFunction((0, 0, 1))) for p in (0, 1, 2)}
    for p, v in D.items():
        assert 0.31 < v < 0.345, f"p={p} C++ uniform demag {v:.5f} not ~1/3"
    assert max(D.values()) - min(D.values()) < 1e-3, f"C++ demag not order-invariant: {D}"


def test_cpp_ho_nonuniform_p_convergence():
    """C++ high-order: non-uniform M=(0,0,z) -> demag p-converges (p=1 == p=2, p=0 distinct)."""
    mesh = _cube(0.8)
    Mcf = ng.CoefficientFunction((0, 0, ng.z))
    D = {p: _cpp_demag(mesh, p, Mcf) for p in (0, 1, 2)}
    assert abs(D[1] - D[2]) < 1e-3, f"C++ p=1,p=2 should agree (linear M.n exact at p>=1): {D}"
    assert abs(D[0] - D[1]) > 1e-3, f"C++ p=0 should differ (SurfaceL2(0) under-resolves linear M.n): {D}"
