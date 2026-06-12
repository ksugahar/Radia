"""radia.hdiv_vim._vim -- an ngsolve.bem-STYLE API for the HDiv-type VIM demag operator.

Mirrors the ngsolve.bem design (SingleLayerPotentialOperator etc.): construct the operator from an NGSolve
FESpace -- the polynomial ORDER comes from the fes, exactly like `HDiv(mesh, order=p)` -- and expose
`.mat`, an H-matrix-backed NGSolve `BaseMatrix` that composes with NGSolve's solvers / BlockMatrix just
like `SingleLayerPotentialOperator(fes, ...).mat`.  RT0 (order=0) and order=p go through ONE call (order=0
is the degenerate constant-monomial case):

    from ngsolve import *
    from ngsolve.krylovspace import GMRes
    from radia.hdiv_vim import DemagOperator

    mesh = Mesh(...)
    fes  = HDiv(mesh, order=p)                       # order from the fes (NGSolve idiom)
    with TaskManager():
        N = DemagOperator(fes, intorder=3*p+6, eps=1e-7)   # like SingleLayerPotentialOperator(fes, ...)
        # N.mat : BaseMatrix == the demag operator B^T G B (G = the C++ charge-Gram H-matrix)
        u, v = fes.TnT()
        M = BilinearForm(u*v*dx).Assemble()          # the HDiv mass
        A = (1.0/chi)*M.mat - N.mat                  # BaseMatrix composition (NGSolve)
        gfm = GridFunction(fes)
        gfm.vec.data = GMRes(A=A, b=rhs.vec, tol=1e-8, maxsteps=400)

Convenience: ``N.DemagFactor(CF((0,0,1)))`` -> the Rayleigh quotient (the demag factor, ~1/3 for a sphere).

Backend: the C++ charge-Gram H-matrix (radia._radia_pybind._ChargeGramHMatrix, the order-p mode merged in
a27d1a5c).  The charge basis is element-local monomials (host reference coords); N = B^T G B is
basis-invariant, so the demag matches the NGSolve-L2-basis dense reference.  Pure Python glue -- no dense
O(N^2) operator is ever formed (the H-matrix gives the O(N log N) matvec).

TaskManager: per the caller-wraps policy, this module does NOT open a TaskManager; the CALLER wraps the
DemagOperator construction + DemagFactor / solve in `with TaskManager():` (the ngsolve.bem idiom).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import ngsolve as ng

import radia._radia_pybind as _rp


# ------------------------------------------------------------------ reference Gauss-Duffy quadrature
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
    return np.array(P), np.array(W)            # ref-tet pts (lam1,lam2,lam3); weights sum 1/6


def _tri_ref(o):
    s, ws = _g01(o)
    P, W = [], []
    for u, wu in zip(s, ws):
        for v, wv in zip(s, ws):
            P.append((u, v * (1 - u)))
            W.append(wu * wv * (1 - u))
    return np.array(P), np.array(W)            # ref-tri pts (lam1,lam2); weights sum 1/2


def _monos_vol(pv):
    return [(i, j, k) for i in range(pv + 1) for j in range(pv + 1 - i) for k in range(pv + 1 - i - j)]


def _monos_surf(p):
    return [(i, j) for i in range(p + 1) for j in range(p + 1 - i)]


def _ngsolve_affine(trafo, eltype, ndim):
    """Reconstruct the (affine) NGSolve element map P = P0 + Jng @ pt from its ElementTransformation, by fitting
    a few IntegrationRule evaluations.  pt is in NGSolve's reference frame."""
    ir = ng.IntegrationRule(eltype, 2)
    rp = np.array([list(p.point)[:ndim] for p in ir])
    Pp = np.array([list(trafo(p).point) for p in ir])
    A = np.hstack([rp, np.ones((len(rp), 1))])
    X, *_ = np.linalg.lstsq(A, Pp, rcond=None)             # [Jng^T (ndim rows) ; P0]
    return X[ndim, :], X[:ndim, :].T                        # (P0, Jng 3xndim)


def _change_of_basis(fe, mons, refP, refW, dim, trafo, Vmesh):
    """L2/SurfaceL2 -> monomial change-of-basis, in the GRAM's cell_verts geometry frame.

    CRITICAL (root cause of the high-order demag bug, 2026-06-13): NGSolve's L2/SurfaceL2 `CalcShape` uses its
    OWN reference-element frame (ref(0,0,0)->the LAST mesh vertex, the standard Netgen ordering), but the C++
    charge-Gram interprets the resulting monomials via `cell_verts` in MESH-VERTEX order (ref(0,0,0)->V0).  If
    the monomials are built in NGSolve's frame (the old code evaluated m_a and CalcShape at the same pt), the
    charge B feeds the Gram is geometrically scrambled by a fixed vertex permutation -- INVISIBLE to every
    uniform-M / demag-factor test (uniform M has div M = 0 => no volume charge) but it makes non-uniform
    (high-order) solves diverge: the demag operator under-counts divM-heavy modes, which the chi*(...) solve
    then amplifies.  So we evaluate the MONOMIAL at the cell_verts-frame coord `g` that corresponds to the
    same physical point as the NGSolve-ref point `pt` (via GetTrafo), keeping CalcShape at `pt`.  This lands
    the monomial coefficients in the Gram's frame.  The map is a fixed reference permutation (NGSolve orders
    local vertices consistently), so computing it from one element of each type is exact for all.

    C[a][k] = INT m_a(g(pt)) phi_k(pt), Mmono[a][b] = INT m_a(g(pt)) m_b(g(pt)); S = Mmono^{-1} C.  The
    quadrature is exact for the polynomial degree, so integrating over pt (vs g) is immaterial."""
    P0, Jng = _ngsolve_affine(trafo, ng.TET if dim == 3 else ng.TRIG, dim)
    if dim == 3:
        Jm = np.array([Vmesh[1] - Vmesh[0], Vmesh[2] - Vmesh[0], Vmesh[3] - Vmesh[0]]).T
    else:
        Jm = np.array([Vmesh[1] - Vmesh[0], Vmesh[2] - Vmesh[0]]).T
    Jm_pinv = np.linalg.pinv(Jm)                            # 3x3 (tet) or 2x3 (tri) -> maps phys offset to ref
    nm, nsh = len(mons), fe.ndof
    M = np.zeros((nm, nm))
    C = np.zeros((nm, nsh))
    for pt, w in zip(refP, refW):
        P = P0 + Jng @ np.array(pt[:dim])                   # physical point NGSolve maps pt to
        g = Jm_pinv @ (P - Vmesh[0])                        # same point's coord in the cell_verts (Gram) frame
        if dim == 3:
            mv = np.array([g[0] ** i * g[1] ** j * g[2] ** k for (i, j, k) in mons])
            sh = np.array(fe.CalcShape(pt[0], pt[1], pt[2]))
        else:
            mv = np.array([g[0] ** i * g[1] ** j for (i, j) in mons])
            sh = np.array(fe.CalcShape(pt[0], pt[1]))
        M += w * np.outer(mv, mv)
        C += w * np.outer(mv, sh)
    return np.linalg.solve(M, C)


def _csr(bf):
    r, c, v = bf.mat.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(bf.mat.height, bf.mat.width))


def build_charge_gram(fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0, far_quad=3, ho_far_factor=2.0,
                      inner_quad=None):
    """From an HDiv FESpace (order p, the order from the fes), build the monomial charge-density map
    B (scipy CSR, n_charge x ndof), the C++ charge-Gram H-matrix G, and the HDiv mass M_mass (CSR).
    order=0 is the degenerate constant-monomial case (== RT0).  The CALLER wraps in TaskManager.

    NEAR/FAR adaptive quadrature (order>0, the DEFAULT build speedup -- accuracy-preserving + golden-locked):
    far charge pairs use a cheap LOW-quad plain double-Gauss in the C++ Gram (the kernel is smooth there),
    near/self pairs keep the full high-quad subtraction.  far_quad = the LOW Gauss pts/dim (<< the near
    `quad`); ho_far_factor = the separation threshold (FAR if |c_a-c_b| > ho_far_factor*(size_a+size_b), the
    size = each host's BOUNDING RADIUS, so touching high-aspect-ratio needle/sliver elements are never
    misclassified FAR).  DEFAULT ho_far_factor=2.0 => the split is ON; it is validated to match the exact
    all-high-quad build to <1e-3 by tests/feec/test_hdiv_vim_nearfar_highorder.py -- so this is NOT a silent
    approximation, it is a TESTED accuracy-preserving quadrature-order choice (a default param == the user's
    contract).  Pass ho_far_factor=inf to FORCE the exact all-high-quad build (e.g. a golden reference)."""
    mesh = fes.mesh
    p = fes.globalorder
    pv = max(p - 1, 0)
    # Gauss pts/dim for the NEAR/SELF singular entries (the far/smooth pairs use the cheaper far_quad).  N =
    # B^T G B is a demag SELF-ENERGY and MUST be positive-semidefinite; UNDER-integrating the near/self pairs
    # makes N NON-PSD (a spurious negative eigenvalue) -- harmless for the demag FACTOR / linear solve (the
    # negative eig is ~1e-3 relative and M_mass dominates ((1/chi)M_mass+N) anyway), but it SILENTLY CORRUPTS
    # any ENERGY / EIGENVALUE use of N (the nonlinear energy-min solve slides down the negative mode -> garbage
    # high-order M).  Measured PSD floor (sphere, min eig of N): quad >= 3*p -- p=1 PSD at quad=3, p=2 needs
    # quad>=6 (quad=4,5 give min eig -3.2e-4 NON-PSD), p=3 needs quad>=9 (quad=8 still NON-PSD).  So the floor
    # is 3*p, NOT the old max(4,p+2) (a prior "floor 6->4" speed opt traded away PSD at p=2).  near cost is
    # O(quad^6) but only on NEAR pairs (far pairs keep far_quad); intorder OVERRIDES (intorder < 2*(3p)-1 may be
    # NON-PSD -- use only for a fast demag-factor estimate, NEVER for an energy/eigenvalue solve).  See
    # tests/feec/test_hdiv_vim_psd.py.
    quad = max(p + 2, (intorder + 1) // 2) if intorder is not None else max(3 * p, 4)
    # INNER subtraction quad (B2 speedup): the subtraction remainder (m_src(y)-m_src(p)) is SMOOTH (the
    # singular part is carried EXACTLY by the analytic PhiTet/TriPotential base), so the inner sum uses a
    # COARSER rule than the outer -> another ~1.5-2x on the O(quad_out^3 * quad_in^3) near entries.  Floor at
    # max(quad-2, p+1); only passed to C++ when iq < quad (else inner = outer).  Validated to hold the same
    # demag accuracy as inner=outer by the nearfar/operator goldens + the uniform-1/3 metric.
    iq = inner_quad if inner_quad is not None else max(quad - 2, p + 1)
    nn = ng.specialcf.normal(mesh.dim)
    L2v, L2b = ng.L2(mesh, order=pv), ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mv = ng.BilinearForm(L2v); mv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; mv.Assemble()
    mb = ng.BilinearForm(L2b); mb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; mb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    # block-diagonal change-of-basis: mv (L2) and mb (SurfaceL2) are DG mass matrices -> BLOCK-DIAGONAL, so
    # sparse spsolve on the CSC is O(N).  The old np.linalg.solve(mv.toarray(), bv.toarray()) was DENSE
    # O(n_cells^3) AND materialized an (n_cells x ndof) dense array -> >300s + GBs at ~5000 tets (the build
    # did NOT scale, even though the matvec does).  spsolve keeps Bv_d/Bb_d SPARSE and the build scalable.
    Bv_d = spla.spsolve(sp.csc_matrix(_csr(mv)), sp.csc_matrix(_csr(bv))).tocsr()   # sparse (vol density map)
    Bb_d = spla.spsolve(sp.csc_matrix(_csr(mb)), sp.csc_matrix(_csr(bb))).tocsr()   # sparse (surf)
    M_mass = _csr(mh)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    vV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in vels]
    bV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in bels]
    vdof = [list(L2v.GetDofNrs(e)) for e in vels]
    bdof = [list(L2b.GetDofNrs(e)) for e in bels]
    mons_v, mons_s = _monos_vol(pv), _monos_surf(p)
    # change-of-basis built in the GRAM's cell_verts frame (see _change_of_basis docstring): pass element 0's
    # NGSolve transform + its mesh-order vertices so the monomial coeffs land in the same frame the C++ Gram
    # interprets `charge_expo` in.  The ref-corner->local-vertex permutation is constant per element type, so
    # one element of each type is exact for all.
    Sv = _change_of_basis(L2v.GetFE(vels[0]), mons_v, *_tet_ref(quad), dim=3,
                          trafo=mesh.GetTrafo(vels[0]), Vmesh=vV[0])
    Ss = _change_of_basis(L2b.GetFE(bels[0]), mons_s, *_tri_ref(quad), dim=2,
                          trafo=mesh.GetTrafo(bels[0]), Vmesh=bV[0])

    Sv_sp, Ss_sp = sp.csr_matrix(Sv), sp.csr_matrix(Ss)         # change-of-basis kept sparse (block x block)
    Brows, host, kind, expo = [], [], [], []
    for c in range(len(vels)):
        blk = Sv_sp @ Bv_d[vdof[c], :]                          # sparse (nmons_v x ndof)
        for a, (i, j, k) in enumerate(mons_v):
            Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    for f in range(len(bels)):
        blk = Ss_sp @ Bb_d[bdof[f], :]                          # sparse (nmons_s x ndof)
        for a, (i, j) in enumerate(mons_s):
            Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B = sp.vstack(Brows).tocsr()                                # (n_charge, ndof)

    rtp, rtw = _tet_ref(quad)
    rsp, rsw = _tri_ref(quad)
    cell_verts = np.concatenate([V.ravel() for V in vV]).tolist()
    face_verts = np.concatenate([V.ravel() for V in bV]).tolist()
    if p == 0:
        # order-0 = CONSTANT charges -> use the FAST ANALYTIC Gram (Wilton/PhiTet, exact), NOT quadrature.
        # The high-order QUADRATURE constructor (Sauter-Schwab over charge pairs) is ~100x slower per
        # H-matrix entry and is pure waste for constant charges -- it was THE DemagOperator build bottleneck
        # (>195s @ 5310 tets; the change-of-basis is only ~4s).  Charge order [cells..., faces...] matches B's
        # row order, so this is the same B^T G B as the validated solve_nonlinear_newton_scalable order-0 path.
        G = _rp._ChargeGramHMatrix(cell_verts=cell_verts, face_verts=face_verts, n_el=len(vels),
                                   eps=eps, leaf=leafsize, eta=eta, near_factor=1e30)
    else:
        # order p>0: monomial-charge quadrature Gram.  By DEFAULT (ho_far_factor=2.0) the near/far split is
        # ON -- well-separated pairs use the cheap LOW-quad QuadDotFar, near/self keep the full high-quad
        # subtraction (validated accuracy-preserving by the nearfar golden).  Pass ho_far_factor=inf to
        # DISABLE the split (every pair high-quad = the exact reference build; the _lo rules are then not
        # built or passed, binding the same way the pre-far-split overload did).
        kw = dict(cell_verts=cell_verts, face_verts=face_verts,
                  n_el=len(vels), charge_host=host, charge_kind=kind, charge_expo=expo,
                  ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
                  ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(),
                  eps=eps, leaf=leafsize, eta=eta)
        if np.isfinite(ho_far_factor):
            rtp_lo, rtw_lo = _tet_ref(far_quad)
            rsp_lo, rsw_lo = _tri_ref(far_quad)
            kw.update(ref_tet_pts_lo=rtp_lo.ravel().tolist(), ref_tet_w_lo=rtw_lo.tolist(),
                      ref_tri_pts_lo=rsp_lo.ravel().tolist(), ref_tri_w_lo=rsw_lo.tolist(),
                      ho_far_factor=ho_far_factor)
        if iq < quad:
            rtp_in, rtw_in = _tet_ref(iq)
            rsp_in, rsw_in = _tri_ref(iq)
            kw.update(ref_tet_pts_in=rtp_in.ravel().tolist(), ref_tet_w_in=rtw_in.tolist(),
                      ref_tri_pts_in=rsp_in.ravel().tolist(), ref_tri_w_in=rsw_in.tolist())
        G = _rp._ChargeGramHMatrix(**kw)
    return B, G, M_mass


class _DemagMat(ng.BaseMatrix):
    """NGSolve BaseMatrix for N = B^T G B (HDiv-VIM demag): N x = B^T (G (B x)), G the C++ charge-Gram
    H-matvec.  Symmetric (Mult == its own transpose).  Composes with NGSolve solvers / BlockMatrix."""
    def __init__(self, fes, B, G):
        super().__init__()
        self._fes = fes
        self._gf = ng.GridFunction(fes)          # template for CreateColVector/RowVector (HDiv has none)
        self._B = B.tocsr()
        self._BT = B.T.tocsr()
        self._G = G

    def IsComplex(self):
        return False

    def Height(self):
        return self._fes.ndof

    def Width(self):
        return self._fes.ndof

    def CreateRowVector(self):
        return self._gf.vec.CreateVector()

    def CreateColVector(self):
        return self._gf.vec.CreateVector()

    def _apply(self, xv):
        c = self._B @ xv
        Gc = np.asarray(self._G.matvec(c.tolist()))
        return self._BT @ Gc

    def Mult(self, x, y):
        y.FV().NumPy()[:] = self._apply(x.FV().NumPy())

    def MultAdd(self, scal, x, y):
        y.FV().NumPy()[:] += scal * self._apply(x.FV().NumPy())

    def MultTransAdd(self, scal, x, y):          # N symmetric -> transpose == itself
        self.MultAdd(scal, x, y)


class DemagOperator:
    """ngsolve.bem-style HDiv-type VIM demag operator.  Construct from an HDiv FESpace; `.mat` is the
    H-matrix-backed NGSolve BaseMatrix N = B^T G B.  See the module docstring for the idiom.  The CALLER
    wraps construction + DemagFactor in `with TaskManager():`."""

    def __init__(self, fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0,
                 far_quad=3, ho_far_factor=2.0, inner_quad=None):
        self.space = fes
        self._B, self._G, self._Mmass = build_charge_gram(
            fes, intorder=intorder, eps=eps, leafsize=leafsize, eta=eta,
            far_quad=far_quad, ho_far_factor=ho_far_factor, inner_quad=inner_quad)
        self.mat = _DemagMat(fes, self._B, self._G)

    @property
    def ndof(self):
        return self.space.ndof

    def DemagFactor(self, M_cf):
        """Rayleigh quotient (the demag factor) for a magnetization CoefficientFunction M_cf:
        <c, G c> / <m, M_mass m>, c = B m, m = the HDiv projection of M_cf.  ~1/3 for a sphere/cube."""
        gfu = ng.GridFunction(self.space)
        gfu.Set(M_cf)
        m = np.asarray(gfu.vec)
        c = self._B @ m
        Gc = np.asarray(self._G.matvec(c.tolist()))
        return float(c @ Gc) / float(m @ (self._Mmass @ m))
