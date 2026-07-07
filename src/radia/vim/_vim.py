"""radia.vim._vim -- an ngsolve.bem-STYLE API for the HDiv-type VIM demag operator.

Mirrors the ngsolve.bem design (SingleLayerPotentialOperator etc.): construct the operator from an NGSolve
FESpace -- the polynomial ORDER comes from the fes, exactly like `HDiv(mesh, order=p)` -- and expose
`.mat`, an H-matrix-backed NGSolve `BaseMatrix` that composes with NGSolve's solvers / BlockMatrix just
like `SingleLayerPotentialOperator(fes, ...).mat`.  RT0 (order=0) and order=p go through ONE call (order=0
is the degenerate constant-monomial case):

    from ngsolve import *
    from ngsolve.krylovspace import GMRes
    from radia.vim import DemagOperator

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


# ---------------- SYMMETRIC degree-5 simplex rules (Keast 1986 15-pt tet / Dunavant 1985 7-pt tri)
# These REPLACE the PRODUCT Gauss-Duffy OUTER charge-Gram quadrature at quad==3 (linear near + default far)
# AND quad==4 (the nonlinear energy-Newton near rule).  KEY FACT (measured): the Duffy-collapse Jacobian
# (1-a)^2(1-b) LOWERS the effective simplex degree of _tet_ref(o) to 2o-3 (NOT the 1D 2o-1) -- so product
# _tet_ref(3)=27pts is only degree 3, and _tet_ref(4)=64pts is degree 5.  Keast-15 is degree 5, so it
# (a) UPGRADES the linear near rule (degree 3 -> 5) at 27->15 pts, and (b) MATCHES the nonlinear rule's
# degree 5 at 64->15 pts (4.3x fewer).  Valid because the INNER integral is carried by the exact analytic
# PhiTet/TriPotential, so the outer integrand is C^{1,alpha} (smooth) even on self/face/edge/vertex-adjacent
# pairs -- a symmetric rule does NOT need the Duffy point-clustering.  Validated (locks below): degree-5
# exact to 1e-17; on the charge-Gram it reproduces demag to <=7e-6, leaves the transverse leak identical,
# and preserves PSD (min eig ~0); the NONLINEAR energy-Newton converges in the SAME or FEWER iterations
# (deep saturation H0=3e6: 15 vs product-64's 17) with M->Msat.  Build ~1.5-2.8x faster (grows with N; the
# nonlinear win is larger since product-64 has 4.3x the near points).  The fully-double-ANALYTIC route was
# surveyed and rejected: no tractable closed form for the dominant tet-tet Galerkin double integral -- the
# symmetric outer rule is the real, cheap lever.  Only the degree-5 pair (quad in {3,4}: linear near +
# default far_quad=3 + nonlinear near) is tabulated; curved, inner-subtraction (iq=2), and any other order
# fall back to the product rule.  The change-of-basis quadrature stays on _tet_ref/_tri_ref (S is exact at
# either rule -> bit-identical).
def _sym_orbit(bary, ncoord):
    """Expand a barycentric orbit to all distinct permutations, dropping the first coord (x,y[,z] = the last
    ncoord barycentric coords); any assignment is equivalent for a symmetric rule."""
    from itertools import permutations
    seen = set()
    for p in permutations(bary):
        if p not in seen:
            seen.add(p)
            yield p[1:1 + ncoord]


def _tet_ref_sym5():
    b = 1.0 / 3.0
    c = 8.0 / 11.0; d = 1.0 / 11.0
    e = 0.0665501535736643; f = 0.4334498464263357
    orbits = [((0.25, 0.25, 0.25, 0.25), 0.1817020685825351),
              ((0.0, b, b, b),            0.0361607142857143),
              ((c, d, d, d),              0.0698714945161738),
              ((e, e, f, f),              0.0656948493683187)]
    P, W = [], []
    for bary, w in orbits:
        for pt in _sym_orbit(bary, 3):
            P.append(pt); W.append(w / 6.0)    # normalize vol-1 -> ref-tet vol 1/6
    return np.array(P), np.array(W)            # (15,3), weights sum 1/6


def _tri_ref_sym5():
    a = 0.4701420641051151; b = 0.1012865073234563
    orbits = [((1.0 / 3, 1.0 / 3, 1.0 / 3), 0.2250000000000000),
              ((a, a, 1 - 2 * a),           0.1323941527885062),
              ((b, b, 1 - 2 * b),           0.1259391805448271)]
    P, W = [], []
    for bary, w in orbits:
        for pt in _sym_orbit(bary, 2):
            P.append(pt); W.append(w * 0.5)    # normalize area-1 -> ref-tri area 1/2
    return np.array(P), np.array(W)            # (7,2), weights sum 1/2


_SYM5_TET = _tet_ref_sym5()
_SYM5_TRI = _tri_ref_sym5()


def _outer_tet(quad):
    """OUTER Gram tet quadrature: symmetric degree-5 (Keast-15) for quad in {3,4} -- quad==3 is the linear
    near + default far rule (product _tet_ref(3) is only degree 3), quad==4 is the nonlinear energy-Newton
    near rule (product _tet_ref(4)=64pts is degree 5, matched by the 15-pt symmetric rule).  Any other order
    (inner subtraction iq=2, intorder overrides, curved) falls back to product Gauss-Duffy."""
    return _SYM5_TET if quad in (3, 4) else _tet_ref(quad)


def _outer_tri(quad):
    """OUTER Gram tri quadrature: symmetric degree-5 (Dunavant-7) for quad in {3,4} (linear + nonlinear near +
    default far); else product Gauss-Duffy."""
    return _SYM5_TRI if quad in (3, 4) else _tri_ref(quad)


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

    CRITICAL (root cause of the high-order demag bug, 2026-06-13, corrected 2026-06-28): NGSolve's L2/SurfaceL2 `CalcShape` uses its
    OWN reference-element frame (ref(0,0,0)->the LAST mesh vertex, the standard Netgen ordering), but the C++
    charge-Gram interprets the resulting monomials via `cell_verts` in MESH-VERTEX order (ref(0,0,0)->V0).  If
    the monomials are built in NGSolve's frame (the old code evaluated m_a and CalcShape at the same pt), the
    charge B feeds the Gram is geometrically scrambled by a fixed vertex permutation -- INVISIBLE to every
    uniform-M / demag-factor test (uniform M has div M = 0 => no volume charge) but it makes non-uniform
    (high-order) solves diverge.  A later M4 audit found the remaining trap: the transform cannot be reused
    from one representative element.  NGSolve orients high-order shapes and the element map by each element's
    global vertex order, so `_charge_basis` must call this helper per element with that element's `trafo` and
    mesh-vertex frame.

    We evaluate the MONOMIAL at the cell_verts-frame coord `g` that corresponds to the same physical point as
    the NGSolve-ref point `pt` (via GetTrafo), keeping CalcShape at `pt`.  This lands the monomial coefficients
    in the Gram's frame for the specific element being processed.

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


def _blockdiag_density_map(M, Bx, fes_dg, vorb, mesh):
    """Solve M @ X = Bx where M is a BLOCK-DIAGONAL DG mass matrix (one block per element = that element's
    DOF list), returning X (CSR).  This is the charge-density map (vol: -div u / mass; surf: u.n / mass).

    M is DG (L2 / SurfaceL2) so it is block-diagonal by construction: a general scipy `spsolve` does a full
    SuperLU of M + one back-substitution per RHS column (n_face columns), which MEASURED at ~1.97 s / ~47%
    of the serial change-of-basis -- pure waste.  We invert the small per-element blocks directly (batched
    np.linalg.inv): L2 order-0 -> 1x1 (M diagonal); SurfaceL2 order-1 -> 3x3.  VERIFIED machine-identical to
    spsolve (max|.| ~6e-14) and ~36x faster (1.97 s -> 0.054 s).  Fail-loud (No-Fallbacks) on non-uniform
    block sizes -- a single-order DG space on one element type has uniform blocks, so this must hold."""
    els = [ng.ElementId(vorb, i) for i in range(mesh.GetNE(vorb))]
    dof_lists = [list(fes_dg.GetDofNrs(e)) for e in els]
    M = sp.csr_matrix(M)
    bs = len(dof_lists[0])
    if not all(len(d) == bs for d in dof_lists):
        raise ValueError("_blockdiag_density_map: non-uniform DG block sizes %s (expected one element type, "
                         "one FE order)" % sorted(set(len(d) for d in dof_lists)))
    if bs == 1:                                            # L2 order 0 -> M is diagonal: X = diag(1/M) Bx
        return (sp.diags(1.0 / M.diagonal()) @ Bx).tocsr()
    idx = np.asarray(dof_lists)                            # (nblk, bs)
    blocks = np.array([M[np.ix_(d, d)].toarray() for d in dof_lists])   # (nblk, bs, bs)
    invb = np.linalg.inv(blocks)
    rows = np.repeat(idx, bs, axis=1).ravel()
    cols = np.tile(idx, (1, bs)).ravel()
    Minv = sp.csr_matrix((invb.ravel(), (rows, cols)), shape=M.shape)
    return (Minv @ Bx).tocsr()


def _charge_basis(fes, quad):
    """Shared geometry + monomial charge-density map for the order-p HDiv-VIM charge operators
    (vim.ChargeGram's analytic Gram AND vim.ChargeGramGauss's point operator): returns B (CSR
    n_charge x ndof), M_mass (CSR), the per-charge (host, kind, flat expo) in the cell_verts reference
    frame, and the host vertex geometry.  `quad` only sets the change-of-basis quadrature (exact for the
    polynomial degree).  CALLER wraps in TaskManager.  Charge order = [cell monomials..., face monomials...]
    (matches B's rows); host[a] is the PER-KIND index (cell idx for kind 0, face idx for kind 1)."""
    mesh = fes.mesh
    p = fes.globalorder
    pv = max(p - 1, 0)
    nn = ng.specialcf.normal(mesh.dim)
    L2v, L2b = ng.L2(mesh, order=pv), ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mv = ng.BilinearForm(L2v); mv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; mv.Assemble()
    mb = ng.BilinearForm(L2b); mb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; mb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    # block-diagonal change-of-basis: mv (L2) and mb (SurfaceL2) are DG mass matrices -> BLOCK-DIAGONAL, so
    # sparse spsolve on the CSC is O(N) (a dense solve was the >300s @ ~5000 tets build bottleneck).
    Bv_d = _blockdiag_density_map(_csr(mv), _csr(bv), L2v, ng.VOL, mesh)   # vol density map (block-diag inv)
    Bb_d = _blockdiag_density_map(_csr(mb), _csr(bb), L2b, ng.BND, mesh)   # surf density map (block-diag inv)
    M_mass = _csr(mh)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    vV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in vels]
    bV = [np.array([mesh[v].point for v in mesh[e].vertices]) for e in bels]
    vdof = [list(L2v.GetDofNrs(e)) for e in vels]
    bdof = [list(L2b.GetDofNrs(e)) for e in bels]
    mons_v, mons_s = _monos_vol(pv), _monos_surf(p)
    # PER-ELEMENT change-of-basis into the GRAM's cell_verts frame.  The map from the L2/SurfaceL2 shape
    # coefficients to the monomial charge the C++ Gram interprets is NOT the same for every element: NGSolve
    # builds the element map AND orients the (high-order) shape functions by the element's GLOBAL vertex order,
    # so both the ref-frame affine map and the shape orientation vary per element.  The old code computed ONE
    # change-of-basis S from element 0 and reused it for all elements -- this geometrically scrambles the
    # monomial charge of every differently-oriented element, injecting a SPURIOUS NET CHARGE (broken neutrality:
    # a real M has INT rho + INT_bnd sigma = 0).  It is INVISIBLE to every uniform-M / demag-FACTOR test
    # (constant charge -> the monopole error cancels), but the spurious monopole has huge Coulomb self-energy,
    # so the high-order demag operator N = B^T G B gets UNPHYSICAL eigenvalues > 1 (a true demag spectrum is in
    # [0,1]: p=1 max 1.22 / 3 modes, p=2 max 8.0 / 15 modes) and the high-order MATERIAL solve over-magnetizes
    # ~2x (p=1) / ~4x (p=2).  Computing S PER ELEMENT restores neutrality, eig in [0,1], and the correct
    # material M (p=1/p=2 then agree with order-0 instead of blowing up).  No cache: S depends on the per-
    # element shape orientation (NOT just the affine map), so a geometry-only key is unsafe; the per-element
    # _change_of_basis is a tiny dense solve, negligible vs the O(N^2) Gram build.  See memory
    # hdiv-highorder-material-solve-wrong round 4.
    rtq, rsq = _tet_ref(quad), _tri_ref(quad)
    Brows, host, kind, expo = [], [], [], []
    if pv == 0:
        # RT1 fast path: the VOLUME change-of-basis is the IDENTITY (Sv == [[1.0]], geometry-INDEPENDENT -- the
        # single constant monomial IS the L2-order-0 constant shape; verified max spread 0.0).  Take the Bv_d rows
        # directly, skipping ~1056 redundant per-cell _change_of_basis calls (~0.8s, the bulk of the loop).
        for c in range(len(vels)):
            Brows.append(Bv_d[vdof[c], :]); host.append(c); kind.append(0); expo += [0, 0, 0]
    else:
        for c in range(len(vels)):
            Sv = sp.csr_matrix(_change_of_basis(L2v.GetFE(vels[c]), mons_v, *rtq, dim=3,
                                                trafo=mesh.GetTrafo(vels[c]), Vmesh=vV[c]))
            blk = Sv @ Bv_d[vdof[c], :]                         # sparse (nmons_v x ndof)
            for a, (i, j, k) in enumerate(mons_v):
                Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    for f in range(len(bels)):
        Ss = sp.csr_matrix(_change_of_basis(L2b.GetFE(bels[f]), mons_s, *rsq, dim=2,
                                            trafo=mesh.GetTrafo(bels[f]), Vmesh=bV[f]))
        blk = Ss @ Bb_d[bdof[f], :]                             # sparse (nmons_s x ndof)
        for a, (i, j) in enumerate(mons_s):
            Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B = sp.vstack(Brows).tocsr()                                # (n_charge, ndof)
    return dict(B=B, M_mass=M_mass, host=host, kind=kind, expo=expo, vV=vV, bV=bV,
                cell_verts=np.concatenate([V.ravel() for V in vV]).tolist(),
                face_verts=np.concatenate([V.ravel() for V in bV]).tolist(),
                mons_v=mons_v, mons_s=mons_s, n_el=len(vels))


# ------------------------------------------------------------------ CURVED (isoparametric P2) charge basis
# P2 reference node positions in the C++ CurvedTet/TriPotential convention (corners then mid-edges).  Extracting
# the P2 nodes via GetTrafo AT these reference positions makes the C++ curved map X(xi) reproduce NGSolve's
# element map, so the monomial xi^e is in NGSolve's reference frame -- the SAME frame CalcShape uses, so the
# curved change-of-basis below is a PURE reference-frame projection (g = pt, no physical round-trip).
_TET_REFNODES = [(0., 0, 0), (1., 0, 0), (0., 1, 0), (0., 0, 1),
                 (.5, 0, 0), (.5, .5, 0), (0., .5, 0), (0., 0, .5), (.5, 0, .5), (0., .5, .5)]
_TRI_REFNODES = [(0., 0), (1., 0), (0., 1), (.5, 0), (.5, .5), (0., .5)]
_IR_TET_NODES = ng.IntegrationRule([tuple(r) for r in _TET_REFNODES], [1.0] * 10)
_IR_TRI_NODES = ng.IntegrationRule([(r[0], r[1]) for r in _TRI_REFNODES], [1.0] * 6)


def _change_of_basis_ref(fe, mons, refP, refW, dim):
    """Reference-frame L2/SurfaceL2 -> monomial change-of-basis for the CURVED Gram: g = the NGSolve ref pt
    directly (the curved C++ map is aligned to NGSolve's reference frame via GetTrafo node extraction, so no
    physical round-trip is needed -- unlike the flat `_change_of_basis`).  Still PER ELEMENT (NGSolve orients
    high-order shapes by global vertex order, picked up by CalcShape)."""
    nm, nsh = len(mons), fe.ndof
    M = np.zeros((nm, nm)); C = np.zeros((nm, nsh))
    for pt, w in zip(refP, refW):
        if dim == 3:
            mv = np.array([pt[0] ** i * pt[1] ** j * pt[2] ** k for (i, j, k) in mons])
            sh = np.array(fe.CalcShape(pt[0], pt[1], pt[2]))
        else:
            mv = np.array([pt[0] ** i * pt[1] ** j for (i, j) in mons])
            sh = np.array(fe.CalcShape(pt[0], pt[1]))
        M += w * np.outer(mv, mv); C += w * np.outer(mv, sh)
    return np.linalg.solve(M, C)


def _charge_basis_curved(fes, quad):
    """CURVED (mesh.Curve(2)) analogue of `_charge_basis`: the charge map B is curved-correct (NGSolve
    integrates -div M / M.n on the curved mesh), the change-of-basis is reference-frame (g=pt), and the
    per-element P2 high-order nodes (10/tet, 6/tri, in the C++ convention) are extracted via GetTrafo.  CALLER
    wraps in TaskManager.  Returns cell_nodes [n_cell*30] / face_nodes [n_bf*18] (flat) + B / M_mass / host /
    kind / expo (same layout as `_charge_basis`)."""
    mesh = fes.mesh
    p = fes.globalorder
    pv = max(p - 1, 0)
    nn = ng.specialcf.normal(mesh.dim)
    L2v, L2b = ng.L2(mesh, order=pv), ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mv = ng.BilinearForm(L2v); mv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; mv.Assemble()
    mb = ng.BilinearForm(L2b); mb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; mb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv_d = _blockdiag_density_map(_csr(mv), _csr(bv), L2v, ng.VOL, mesh)
    Bb_d = _blockdiag_density_map(_csr(mb), _csr(bb), L2b, ng.BND, mesh)
    M_mass = _csr(mh)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    vdof = [list(L2v.GetDofNrs(e)) for e in vels]
    bdof = [list(L2b.GetDofNrs(e)) for e in bels]
    mons_v, mons_s = _monos_vol(pv), _monos_surf(p)
    rtp, rtw = _tet_ref(quad); rsp, rsw = _tri_ref(quad)

    Brows, host, kind, expo = [], [], [], []
    cell_nodes, face_nodes = [], []
    for c, e in enumerate(vels):
        cell_nodes.append(_trafo_lattice_nodes(mesh, e, _IR_TET_NODES))   # P2 nodes (curved geom, kept)
        if pv == 0:                                            # RT1: volume Sv == [[1]] (identity) -> Bv_d row direct
            Brows.append(Bv_d[vdof[c], :]); host.append(c); kind.append(0); expo += [0, 0, 0]
        else:
            Sv = sp.csr_matrix(_change_of_basis_ref(L2v.GetFE(e), mons_v, rtp, rtw, dim=3))
            blk = Sv @ Bv_d[vdof[c], :]
            for a, (i, j, k) in enumerate(mons_v):
                Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    for f, e in enumerate(bels):
        face_nodes.append(_trafo_lattice_nodes(mesh, e, _IR_TRI_NODES))
        Ss = sp.csr_matrix(_change_of_basis_ref(L2b.GetFE(e), mons_s, rsp, rsw, dim=2))
        blk = Ss @ Bb_d[bdof[f], :]
        for a, (i, j) in enumerate(mons_s):
            Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B = sp.vstack(Brows).tocsr()
    return dict(B=B, M_mass=M_mass, host=host, kind=kind, expo=expo,
                cell_nodes=np.concatenate([V.ravel() for V in cell_nodes]).tolist() if cell_nodes else [],
                face_nodes=np.concatenate([V.ravel() for V in face_nodes]).tolist() if face_nodes else [],
                n_el=len(vels))


# ----------------------------------------------------------- HEX RT1 (pure-hex HDiv order-1) charge basis ---
# The hex volume charge = -div(HDiv-hex order-1) lives in L2(order=1) = Q1 trilinear (8 modes/hex) -- the "hex
# gotcha" (NOT the tet's P0; testing div against only P0 pollutes ker(N) -> the demag solve blows up ~24x).
# Geometry = 27-node triquadratic (Q2) lattice via GetTrafo -> ONE code path for FLAT (trilinear subset of Q2,
# exact) AND CURVED (mesh.Curve(2)) hexes.  Face charge = SurfaceL2(order=1) bilinear (4/quad-face).  The C++
# RadHACApKChargeGram hex mode decomposes each hex into 6 sub-tets (quad face -> 2 sub-tris) as an INTEGRATION
# device and does the both-domains-graded Duffy singular quadrature that keeps eig(M_mass^-1 N) <= 1.  See
# the de-risk + block-memo perf fix in memory hdiv-tet-hex-coupling-pyramid-gated.
_MONS_HEX = [(i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)]      # 8 Q1 monomials, x-fastest
_MONS_QUAD = [(i, j) for j in (0, 1) for i in (0, 1)]                        # 4 Q1 monomials
_Q2_LATTICE_3D = [(ix / 2.0, iy / 2.0, iz / 2.0)
                  for iz in range(3) for iy in range(3) for ix in range(3)]  # n = ix + 3*iy + 9*iz
_Q2_LATTICE_2D = [(iu / 2.0, iv / 2.0) for iv in range(3) for iu in range(3)]  # n = iu + 3*iv


def _ref_prod_gauss(n, dim):
    """Product Gauss rule on the reference hex [0,1]^dim (for the Q1 change-of-basis projection)."""
    g, gw = _g01(n)
    P, W = [], []
    if dim == 3:
        for a, wa in zip(g, gw):
            for b, wb in zip(g, gw):
                for c, wc in zip(g, gw):
                    P.append((a, b, c)); W.append(wa * wb * wc)
    else:
        for a, wa in zip(g, gw):
            for b, wb in zip(g, gw):
                P.append((a, b)); W.append(wa * wb)
    return np.array(P), np.array(W)


def _trafo_lattice_nodes(mesh, e, ir, max_tries=16):
    """GetTrafo lattice extraction with a DETERMINISM CONTRACT (2026-07-03): evaluate the element's
    lattice until two CONSECUTIVE evaluations agree bit-for-bit (and are finite), fail LOUD otherwise.

    Why: the scalar `tr(ip).point` path occasionally returns UNINITIALIZED memory (~1e-310 denormals in
    single coordinates of a few lattice points) on the FIRST process-wide basis extraction -- caught in
    the wild at ~2-8%% of fresh processes (load-sensitive) on the affine wiring golden, poisoning the
    Gram geometry while B/M_mass (assembled via the SIMD mapped-rule path) stay bit-exact.  The
    corruption window can persist across several re-evaluations AND is bursty under heavy machine
    load (one loaded-machine run saw an element stay unstable for 8 consecutive back-to-back tries),
    so the contract retries generously WITH a short decorrelation sleep after each mismatch; two
    consecutive bit-identical finite evaluations are accepted, anything else raises -- nothing is
    silently accepted.  See memory ngsolve-gettrafo-first-touch-garbage."""
    import time as _time
    prev = None
    for k in range(max_tries):
        tr = mesh.GetTrafo(e)
        cur = np.array([list(tr(ip).point) for ip in ir])
        if prev is not None and np.array_equal(prev, cur) and np.all(np.isfinite(cur)):
            return cur
        prev = cur
        if k >= 1:
            _time.sleep(0.002)      # decorrelate the bursty corruption window
    raise RuntimeError(
        f"GetTrafo lattice evaluation unstable for element {e} after {max_tries} tries "
        "(NGSolve scalar trafo path returned differing node coordinates -- do not trust this mesh "
        "extraction; rerun, and report the incident).")


def _charge_basis_hex(fes, cob_quad=3):
    """HEX analogue of `_charge_basis_curved`: charge map B + 27/9-node Q2 geometry nodes (via GetTrafo ->
    flat + curved ONE path).  fes = HDiv(hexmesh, order=1).  CALLER wraps TaskManager.

    PIOLA-EXACT charge model (the warped-hex correctness fix): on a mapped hex the TRUE volume charge is
    rho(x) = q_ref(xi)/J(xi) with q_ref = -div_ref(u_ref) in Q1(xi) (and sigma = (u.n)_ref(u,v)/Js on faces)
    -- the 1/J is the Piola transform, NOT representable by ref-frame polynomials.  Projecting rho onto
    ref-Q1 (the old model) is not Coulomb-orthogonal, so the projected charge's self-energy can EXCEED the
    true one -> the demag spectrum leaked above 1 by O(warp^2) on strongly distorted hexes (real cylinder
    mesh eig ~1.01 converged; affine boxes were exact because J = const there, and the tet path is immune
    because tets are affine).  The fix represents the charge EXACTLY: B extracts the ref-frame Q1
    coefficients of q_ref via  INT (-div u) phi dx == INT (-div_ref u_ref) phi dxi  (the J's cancel), solved
    with the REF-measure L2 mass; the C++ Gram then integrates  INT INT m_a(xi) m_b(eta)/|X(xi)-X(eta)|
    dxi deta  with NO Jacobian factors (they cancel against the two 1/J densities)."""
    mesh = fes.mesh
    assert fes.globalorder == 1, "RT1-hex is HDiv order-1 only"
    nn = ng.specialcf.normal(mesh.dim)
    L2v = ng.L2(mesh, order=1)               # HEX GOTCHA: div_ref(HDiv-hex order1) is Q1 -> L2 order 1
    L2b = ng.SurfaceL2(mesh, order=1)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv = _csr(bv)                            # REF-measure moments of q_ref (the physical J cancels)
    Bb = _csr(bb)
    M_mass = _csr(mh)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    rhp, rhw = _ref_prod_gauss(cob_quad, 3)
    rqp, rqw = _ref_prod_gauss(cob_quad, 2)
    ir_hex = ng.IntegrationRule(_Q2_LATTICE_3D, [1.0] * 27)
    ir_quad = ng.IntegrationRule(_Q2_LATTICE_2D, [1.0] * 9)

    Brows, host, kind, expo = [], [], [], []
    cell_nodes, face_nodes = [], []
    for c, e in enumerate(vels):
        cell_nodes.append(_trafo_lattice_nodes(mesh, e, ir_hex))
        fe = L2v.GetFE(e)
        Phi = np.array([fe.CalcShape(*pt) for pt in rhp])            # (nq, nphi) ref shapes
        Mref = (Phi * rhw[:, None]).T @ Phi                          # INT_ref phi phi dxi (exact, Q1 shapes)
        rows = np.linalg.solve(Mref, Bv[list(L2v.GetDofNrs(e)), :].toarray())
        Sv = _change_of_basis_ref(fe, _MONS_HEX, rhp, rhw, dim=3)
        blk = sp.csr_matrix(Sv @ rows)
        for a, (i, j, k) in enumerate(_MONS_HEX):
            Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    n_el = len(vels)
    for f, e in enumerate(bels):
        face_nodes.append(_trafo_lattice_nodes(mesh, e, ir_quad))
        fe = L2b.GetFE(e)
        Phi = np.array([fe.CalcShape(pt[0], pt[1]) for pt in rqp])
        Mref = (Phi * rqw[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bb[list(L2b.GetDofNrs(e)), :].toarray())
        Ss = _change_of_basis_ref(fe, _MONS_QUAD, rqp, rqw, dim=2)
        blk = sp.csr_matrix(Ss @ rows)
        for a, (i, j) in enumerate(_MONS_QUAD):
            Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B = sp.vstack(Brows).tocsr()
    return dict(B=B, M_mass=M_mass, host=host, kind=kind, expo=expo, n_el=n_el, n_bf=len(bels),
                cell_nodes=np.concatenate([n.ravel() for n in cell_nodes]).tolist(),
                face_nodes=(np.concatenate([n.ravel() for n in face_nodes]).tolist() if face_nodes else []))


def _build_charge_gram_hex(fes, glout_n=6, glin_n=5, near_grade=0.6, far_inner=1.5,
                           eps=1e-12, leafsize=64, eta=2.0, image_masks=None, image_signs=None):
    """Pure-hex RT1 charge Gram via the hex-mode C++ _ChargeGramHMatrix.  FLAT and CURVED (mesh.Curve(2))
    share ONE path (the 27-node Q2 lattice is extracted via GetTrafo either way -- the caller Curve(2)'s the
    mesh for curved).  glin_n = the 1D rule of the REF-frame RADIAL near/self inner (the PhiAtHO_Duffy port
    -- robust on distorted/curved hexes, where graded clouds left eig > 1 on the real cylinder mesh);
    eig(M_mass^-1 N) <= 1 gated on box AND cylinder meshes; block-memo build (~59x vs naive per-entry).
    far_inner = the PER-OUTER-POINT radial reach: an outer point farther than far_inner*size from a source
    sub-simplex integrates it with the cheap CACHED far cloud instead of the per-point radial rule.  1.5
    keeps the radial exactly where the kernel peaks (self + the facing side of touching neighbours); the
    2026-07-03 cost attribution measured the 4.0 shell as ~60%% of the whole build on the real cylinder
    for <=1e-4 entry drift.  The far TET cloud is the SAME Keast-15 degree-5 rule as the outer (a degree-3
    WV rule at reach 1.5 ate the flat-cylinder eig margin, 1.0005 -> 1.0044; Keast-15 restores 1.0006 for
    +4%% build).  The build also skips the strictly-lower H-matrix leaves (symmetric fill -- every apply of
    the Gram routes through the exactly-symmetric matvec, so they are never read)."""
    cb = _charge_basis_hex(fes)
    glo, gwo = _g01(glout_n)
    gli, gwi = _g01(glin_n)
    ftp = np.asarray(_SYM5_TET[0]); ftw = np.asarray(_SYM5_TET[1])
    G = _rp._ChargeGramHMatrix(
        hex_cell_nodes=cb["cell_nodes"], quad_face_nodes=cb["face_nodes"],
        n_el=int(cb["n_el"]), n_bf=int(cb["n_bf"]),
        charge_host=list(cb["host"]), charge_kind=list(cb["kind"]), charge_expo=list(cb["expo"]),
        sym_tet_pts=np.asarray(_SYM5_TET[0]).ravel().tolist(), sym_tet_w=np.asarray(_SYM5_TET[1]).tolist(),
        sym_tri_pts=np.asarray(_SYM5_TRI[0]).ravel().tolist(), sym_tri_w=np.asarray(_SYM5_TRI[1]).tolist(),
        gl_out=glo.tolist(), gw_out=gwo.tolist(), gl_in=gli.tolist(), gw_in=gwi.tolist(),
        far_tet_pts=ftp.ravel().tolist(), far_tet_w=ftw.tolist(),
        far_tri_pts=np.asarray(_SYM5_TRI[0]).ravel().tolist(), far_tri_w=np.asarray(_SYM5_TRI[1]).tolist(),
        near_grade=near_grade, far_inner_factor=far_inner,
        image_masks=([] if image_masks is None else list(image_masks)),
        image_signs=([] if image_signs is None else list(image_signs)),
        eps=eps, leaf=leafsize, eta=eta)
    chk = G.hex_state_check()
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "hex charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption "
            "(0xc0000374 class) -- do NOT trust this Gram; rerun, and report the incident.")
    return cb["B"], G, cb["M_mass"]


# 2D PLANAR (motor cross-section) layer -- ref lattices handed to GetTrafo (order matches the C++ maps:
# Tri6Map quadratic barycentric v0,v1,v2,m01,m12,m20; Quad9Map lag3 x-fastest; Edge3Map ends+mid).
_TRI6_LAT = [(1, 0), (0, 1), (0, 0), (0.5, 0.5), (0, 0.5), (0.5, 0)]
_QUAD9_LAT = [(i/2, j/2) for j in range(3) for i in range(3)]
_EDGE3_LAT = [0.0, 0.5, 1.0]
_MONS_TRI2D = [(0, 0)]                                   # div(HDiv order-1 tri) is P0 (probed 2026-07-03)
_MONS_QUAD2D = [(0, 0), (1, 0), (0, 1), (1, 1)]          # div(HDiv order-1 quad) is Q1 (the hex-gotcha twin)
_MONS_EDGE2D = [(0,), (1,)]                              # (u.n)_ref on an edge is P1


def _charge_basis_2d(fes, cob_quad=3):
    """2D analogue of `_charge_basis_hex` (motor cross-sections; memory hdiv-vim-tri-quad-motor): charge
    map B + P2 lattice geometry for tri/quad cells and boundary edges, all in the NGSolve REF frame with
    the Piola-exact extraction (the dimension-independent J-cancellation identity).  Kernel side is the
    2D log Gram (C++ dim2 mode).  CALLER wraps TaskManager."""
    mesh = fes.mesh
    assert fes.globalorder == 1, "2D HDiv-VIM is HDiv order-1 only"
    nn2 = ng.specialcf.normal(2)
    L2v = ng.L2(mesh, order=1)                # Q1-capable: covers quad Q1; the tri P0 monomial sits inside
    Sb2 = ng.SurfaceL2(mesh, order=1)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=Sb2); bb += (u.Trace() * nn2) * Sb2.TestFunction() * ng.ds; bb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv = _csr(bv); Bb = _csr(bb); M_mass = _csr(mh)

    g, gw = _g01(cob_quad)
    tp = np.array([[uu, vv*(1 - uu)] for uu in g for vv in g])            # Duffy on the NGSolve tri ref
    tw = np.array([wu*wv*(1 - uu) for uu, wu in zip(g, gw) for vv, wv in zip(g, gw)])
    qp = np.array([[uu, vv] for uu in g for vv in g])
    qw = np.array([wu*wv for _, wu in zip(g, gw) for _, wv in zip(g, gw)])
    ep = g.reshape(-1, 1); ew = gw

    ir_tri6 = ng.IntegrationRule([(p[0], p[1], 0) for p in _TRI6_LAT], [1.0]*6)
    ir_quad9 = ng.IntegrationRule([(p[0], p[1], 0) for p in _QUAD9_LAT], [1.0]*9)
    ir_edge3 = ng.IntegrationRule([(t, 0, 0) for t in _EDGE3_LAT], [1.0]*3)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    Brows, host, kind, expo = [], [], [], []
    cell_nodes9, cell_type, edge_nodes3 = [], [], []
    for c, e in enumerate(vels):
        quad = len(mesh[e].vertices) == 4
        nodes = _trafo_lattice_nodes(mesh, e, ir_quad9 if quad else ir_tri6)[:, :2]
        slot = np.zeros((9, 2))
        slot[:nodes.shape[0]] = nodes
        cell_nodes9.append(slot); cell_type.append(1 if quad else 0)
        fe = L2v.GetFE(e)
        rp, rw = (qp, qw) if quad else (tp, tw)
        Phi = np.array([fe.CalcShape(pt[0], pt[1], 0.0) for pt in rp])
        Mref = (Phi * rw[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bv[list(L2v.GetDofNrs(e)), :].toarray())
        mons = _MONS_QUAD2D if quad else _MONS_TRI2D
        Sv = _change_of_basis_ref(fe, mons, rp, rw, dim=2)
        blk = Sv @ rows
        for a, (mi, mj) in enumerate(mons):
            Brows.append(sp.csr_matrix(blk[a])); host.append(c); kind.append(0); expo += [mi, mj, 0]
    for f, e in enumerate(bels):
        edge_nodes3.append(_trafo_lattice_nodes(mesh, e, ir_edge3)[:, :2])
        fe = Sb2.GetFE(e)
        Phi = np.array([fe.CalcShape(t[0], 0.0, 0.0) for t in ep])
        Mref = (Phi * ew[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bb[list(Sb2.GetDofNrs(e)), :].toarray())
        # dim=1 change of basis (S = M_mono^-1 C), mirroring _change_of_basis_ref
        Mm = np.array([[np.sum(ew * ep[:, 0]**(mi + mj)) for (mj,) in _MONS_EDGE2D] for (mi,) in _MONS_EDGE2D])
        Cm = np.array([[np.sum(ew * (ep[:, 0]**mi) * Phi[:, k]) for k in range(Phi.shape[1])]
                       for (mi,) in _MONS_EDGE2D])
        Se = np.linalg.solve(Mm, Cm)
        blk = Se @ rows
        for a, (mi,) in enumerate(_MONS_EDGE2D):
            Brows.append(sp.csr_matrix(blk[a])); host.append(f); kind.append(1); expo += [mi, 0, 0]
    B = sp.vstack(Brows).tocsr()
    return dict(B=B, M_mass=M_mass, host=host, kind=kind, expo=expo,
                n_el=len(vels), n_be=len(bels),
                cell_nodes9=np.concatenate([n.ravel() for n in cell_nodes9]).tolist(),
                cell_type=cell_type,
                edge_nodes3=np.concatenate([n.ravel() for n in edge_nodes3]).tolist())


def _prod_tri01(n):
    """product-Gauss bary rule on the unit simplex (lam1 = u, lam2 = v(1-u); W sums 1/2)."""
    g, gw = _g01(n)
    P, W = [], []
    for u, wu in zip(g, gw):
        for v, wv in zip(g, gw):
            P.append((u, v*(1 - u))); W.append(wu*wv*(1 - u))
    return np.array(P), np.array(W)


def _build_charge_gram_2d(fes, outer_n=4, glin_n=8, gledge_n=12, near_grade=0.6, far_inner=1.5,
                          eps=1e-12, leafsize=64, eta=2.0):
    """2D planar charge Gram via the C++ dim2 _ChargeGramHMatrix (kernel -ln(r)/(2pi)).  Regular
    (ungraded) outer everywhere -- but it MUST be the PRODUCT-GAUSS rule (outer_n^2/sub-tri), NOT a
    sparse symmetric rule: the outer integrand m_a(xi)*Phi_b(X(xi)) has C1 kinks where the source charge
    lives, and Dunavant-7's coherent misintegration of those kinks -- amplified through the div-scaled
    charge map (|q| ~ 1/h, components cancel) -- leaked the quad-mesh demag spectrum to eig 1.072 while
    every INDIVIDUAL entry looked fine at ~3e-5 (measured 2026-07-03; product-Gauss 4/6/8 all give
    0.9993).  Inner: radial cones split-graded at the kernel peak (glin_n per dim; edges split-grade at
    the projection parameter -- endpoint grading had the same coherent-overestimate disease), far cloud
    otherwise.  Gates: eig in [0,1] on tri/quad/distorted/disk/ellipse; disk demag 0.50000; ellipse
    0.3344/0.6656 vs 1/3, 2/3; 2D Clausius-Mossotti solve to 2-3e-4."""
    cb = _charge_basis_2d(fes)
    otp, otw = _prod_tri01(outer_n)
    gli, gwi = _g01(glin_n)
    gle, gwe = _g01(gledge_n)
    G = _rp._ChargeGramHMatrix(
        dim2=2,
        cell_nodes9=cb["cell_nodes9"], cell_type=list(cb["cell_type"]), edge_nodes3=cb["edge_nodes3"],
        n_el=int(cb["n_el"]), n_be=int(cb["n_be"]),
        charge_host=list(cb["host"]), charge_kind=list(cb["kind"]), charge_expo=list(cb["expo"]),
        sym_tri_pts=otp.ravel().tolist(), sym_tri_w=otw.tolist(),
        gl_edge=gle.tolist(), gw_edge=gwe.tolist(), gl_in=gli.tolist(), gw_in=gwi.tolist(),
        far_tri_pts=np.asarray(_SYM5_TRI[0]).ravel().tolist(), far_tri_w=np.asarray(_SYM5_TRI[1]).tolist(),
        near_grade=near_grade, far_inner_factor=far_inner, eps=eps, leaf=leafsize, eta=eta)
    chk = G.hex_state_check()
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "2D charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption "
            "(0xc0000374 class) -- do NOT trust this Gram; rerun, and report the incident.")
    return cb["B"], G, cb["M_mass"]


# WEDGE (PRISM) RT1 (2026-07-04, memory hdiv-tet-hex-coupling-pyramid-gated): the prism div-image is
# L2(prism,order=1) = tri-P1 (x) z-P1 = {1,x,y,z,xz,yz} (6/prism, a SUBSET of the hex's 8 Q1 monomials);
# boundary faces are MIXED tri (SurfaceL2 P1, 3) + quad (SurfaceL2 Q1, 4).  Geometry = the 18-node tri-P2
# (x) z-P2 lattice (n = t + 6*iz, t = _TRI6_LAT node) via GetTrafo -> flat + curved ONE path.
_MONS_WEDGE = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]   # tri-P1 (x) z-P1
_MONS_TRI3D = [(0, 0), (1, 0), (0, 1)]                                             # SurfaceL2(order1) on a tri = P1
_WEDGE_Q2_LATTICE = [(_TRI6_LAT[t][0], _TRI6_LAT[t][1], iz / 2.0) for iz in range(3) for t in range(6)]


def _prism_cob_quad(nz=3):
    """Ref-prism change-of-basis quadrature: SYM5 tri (u,v) x nz-pt Gauss (w).  Exact for the degree<=4
    prism-L2-shape x monomial products (the prism L2 order-1 shapes and _MONS_WEDGE are both tri-degree<=1
    (x) z-degree<=1)."""
    tp, tw = np.asarray(_SYM5_TRI[0]), np.asarray(_SYM5_TRI[1])
    gz, gwz = _g01(nz)
    P, W = [], []
    for (tu, tv), twt in zip(tp, tw):
        for z, wz in zip(gz, gwz):
            P.append((tu, tv, z)); W.append(twt * wz)
    return np.array(P), np.array(W)


def _charge_basis_wedge(fes):
    """WEDGE (prism) analogue of `_charge_basis_hex`: charge map B + 18-node prism cell nodes + MIXED
    tri(6-node)/quad(9-node) face nodes (packed in 27-double 9-node slots, a tri fills the first 6) + a
    per-face type array, all via GetTrafo (flat + curved ONE path).  fes = HDiv(prismmesh, order=1).
    CALLER wraps TaskManager.  Piola-exact charge model, exactly as the hex path (the J's cancel)."""
    mesh = fes.mesh
    assert fes.globalorder == 1, "RT1-wedge is HDiv order-1 only"
    nn = ng.specialcf.normal(mesh.dim)
    L2v = ng.L2(mesh, order=1)               # div_ref(HDiv-prism order1) = tri-P1 (x) z-P1 = 6 monomials
    L2b = ng.SurfaceL2(mesh, order=1)        # tri face P1 (3) / quad face Q1 (4)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv = _csr(bv); Bb = _csr(bb); M_mass = _csr(mh)

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = [ng.ElementId(ng.BND, i) for i in range(mesh.GetNE(ng.BND))]
    php, phw = _prism_cob_quad()                                  # prism ref change-of-basis quadrature
    qp2, qw2 = _ref_prod_gauss(3, 2)                             # quad-face ref quadrature ([0,1]^2)
    tp2, tw2 = np.asarray(_SYM5_TRI[0]), np.asarray(_SYM5_TRI[1])  # tri-face ref quadrature (ref tri)
    ir_wedge = ng.IntegrationRule(_WEDGE_Q2_LATTICE, [1.0] * 18)
    ir_quad = ng.IntegrationRule(_Q2_LATTICE_2D, [1.0] * 9)
    ir_tri = ng.IntegrationRule(_TRI6_LAT, [1.0] * 6)

    Brows, host, kind, expo = [], [], [], []
    cell_nodes, face_nodes, face_type = [], [], []
    for c, e in enumerate(vels):
        cell_nodes.append(_trafo_lattice_nodes(mesh, e, ir_wedge))    # (18, 3)
        fe = L2v.GetFE(e)
        Phi = np.array([fe.CalcShape(*pt) for pt in php])
        Mref = (Phi * phw[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bv[list(L2v.GetDofNrs(e)), :].toarray())
        Sv = _change_of_basis_ref(fe, _MONS_WEDGE, php, phw, dim=3)
        blk = sp.csr_matrix(Sv @ rows)
        for a, (i, j, k) in enumerate(_MONS_WEDGE):
            Brows.append(blk[a]); host.append(c); kind.append(0); expo += [i, j, k]
    n_el = len(vels)
    for f, e in enumerate(bels):
        if len(list(mesh[e].vertices)) == 3:                          # TRI face -> 6-node, P1, 1 sub-tri
            nd = _trafo_lattice_nodes(mesh, e, ir_tri)                # (6, 3)
            slot = np.zeros((9, 3)); slot[:6] = nd
            face_nodes.append(slot); face_type.append(0)
            fe = L2b.GetFE(e)
            Phi = np.array([fe.CalcShape(pt[0], pt[1]) for pt in tp2])
            Mref = (Phi * tw2[:, None]).T @ Phi
            rows = np.linalg.solve(Mref, Bb[list(L2b.GetDofNrs(e)), :].toarray())
            Ss = _change_of_basis_ref(fe, _MONS_TRI3D, tp2, tw2, dim=2)
            blk = sp.csr_matrix(Ss @ rows)
            for a, (i, j) in enumerate(_MONS_TRI3D):
                Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
        else:                                                         # QUAD face -> 9-node, Q1, 2 sub-tris
            face_nodes.append(_trafo_lattice_nodes(mesh, e, ir_quad)) # (9, 3)
            face_type.append(1)
            fe = L2b.GetFE(e)
            Phi = np.array([fe.CalcShape(pt[0], pt[1]) for pt in qp2])
            Mref = (Phi * qw2[:, None]).T @ Phi
            rows = np.linalg.solve(Mref, Bb[list(L2b.GetDofNrs(e)), :].toarray())
            Ss = _change_of_basis_ref(fe, _MONS_QUAD, qp2, qw2, dim=2)
            blk = sp.csr_matrix(Ss @ rows)
            for a, (i, j) in enumerate(_MONS_QUAD):
                Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
    B = sp.vstack(Brows).tocsr()
    return dict(B=B, M_mass=M_mass, host=host, kind=kind, expo=expo, n_el=n_el, n_bf=len(bels),
                cell_nodes=np.concatenate([n.ravel() for n in cell_nodes]).tolist(),
                face_nodes=np.concatenate([n.ravel() for n in face_nodes]).tolist(),
                face_type=face_type)


def _build_charge_gram_wedge(fes, glout_n=6, glin_n=5, near_grade=0.6, far_inner=1.5,
                             eps=1e-12, leafsize=64, eta=2.0, image_masks=None, image_signs=None):
    """Pure-prism RT1 charge Gram via the wedge-mode C++ _ChargeGramHMatrix (mirror of _build_charge_gram_hex;
    FLAT + Curve(2) share ONE path).  numpy de-risk eig(M_mass^-1 N) in [0,1]: 0.989 @ n=2, 0.997 @ n=3;
    demag_z ~ 1/3.  The wedge mode shares the hex block memo / symmetric-fill build, so the golden hex path
    is byte-for-byte untouched."""
    cb = _charge_basis_wedge(fes)
    glo, gwo = _g01(glout_n); gli, gwi = _g01(glin_n)
    G = _rp._ChargeGramHMatrix(
        wedge_cell_nodes=cb["cell_nodes"], face_nodes=cb["face_nodes"], face_type=list(cb["face_type"]),
        n_el=int(cb["n_el"]), n_bf=int(cb["n_bf"]),
        charge_host=list(cb["host"]), charge_kind=list(cb["kind"]), charge_expo=list(cb["expo"]),
        sym_tet_pts=np.asarray(_SYM5_TET[0]).ravel().tolist(), sym_tet_w=np.asarray(_SYM5_TET[1]).tolist(),
        sym_tri_pts=np.asarray(_SYM5_TRI[0]).ravel().tolist(), sym_tri_w=np.asarray(_SYM5_TRI[1]).tolist(),
        gl_out=glo.tolist(), gw_out=gwo.tolist(), gl_in=gli.tolist(), gw_in=gwi.tolist(),
        far_tet_pts=np.asarray(_SYM5_TET[0]).ravel().tolist(), far_tet_w=np.asarray(_SYM5_TET[1]).tolist(),
        far_tri_pts=np.asarray(_SYM5_TRI[0]).ravel().tolist(), far_tri_w=np.asarray(_SYM5_TRI[1]).tolist(),
        near_grade=near_grade, far_inner_factor=far_inner,
        image_masks=([] if image_masks is None else list(image_masks)),
        image_signs=([] if image_signs is None else list(image_signs)),
        eps=eps, leaf=leafsize, eta=eta)
    chk = G.hex_state_check()
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "wedge charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption (0xc0000374 class) -- "
            "do NOT trust this Gram; rerun, and report the incident.")
    return cb["B"], G, cb["M_mass"]


def build_charge_gram(fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0, far_quad=3, ho_far_factor=2.0,
                      inner_quad=None, curve_order=None, curve_gauss=8, nonlinear=False,
                      image_masks=None, image_signs=None):
    """From an HDiv FESpace (order p, the order from the fes), build the monomial charge-density map
    B (scipy CSR, n_charge x ndof), the C++ charge-Gram H-matrix G, and the HDiv mass M_mass (CSR).
    order=0 is the degenerate constant-monomial case (== RT0).  The CALLER wraps in TaskManager.

    curve_order (None=flat, or 2=isoparametric P2): when set, build the CURVED charge Gram on the
    mesh.Curve(curve_order) geometry -- curved charge map B (reference-frame change-of-basis) + the C++ curved
    Duffy Gram (curve_gauss = the inner Gauss-Legendre pts/dim, ~8 -> Duffy ~1e-4).  curve_order helps
    near-surface FIELD / FLUX accuracy (sigma=M.n on the true curved surface), NOT the demag FACTOR (which is
    curving-insensitive on a sphere, ~3e-5 in the de-risk sweep).  Only P2
    (curve_order=2) is wired; the mesh MUST already be mesh.Curve(2)'d by the caller.

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
    if p != 1:
        raise ValueError(
            "vim.ChargeGram: HDiv-VIM is RT1 (HDiv order=1) only -- RT0 (order=0) is retired (per-element "
            "inaccurate) and RT2+ is retired (no "
            "per-element gain over RT1, slower).  Build the FESpace as HDiv(mesh, order=1).  (The geometry "
            "curve_order is a SEPARATE knob: curve_order=2 isoparametric P2 is still allowed.)")
    image_masks = [] if image_masks is None else list(image_masks)   # robust for NumPy arrays (truth-value)
    image_signs = [] if image_signs is None else list(image_signs)
    if len(image_masks) != len(image_signs):
        raise ValueError("vim.ChargeGram: image_masks and image_signs must have the same length")
    if image_masks:
        # IMA (mirror-image charge folding): wired for the FLAT pure-TET (C++ m_highorder QuadDotRefl->PhiInner)
        # AND pure-HEX / pure-WEDGE (the QuadBlockHex/Wedge(mask) reflected block) RT1 Grams.  Fail loud on the
        # still-unwired 2D-planar + CURVED reduced models (No-Fallbacks) rather than silently drop the image.
        _ivt = {len(el.vertices) for el in mesh.Elements(ng.VOL)} if mesh.dim == 3 else None
        _curved = (curve_order is not None) or (mesh.dim == 3 and mesh.GetCurveOrder() >= 2)
        if mesh.dim != 3 or _ivt not in ({4}, {8}, {6}) or _curved:
            raise ValueError(
                "vim.ChargeGram: image_masks (IMA) is wired for the FLAT pure-TET / pure-HEX / pure-WEDGE "
                "RT1 Gram only; 2D-planar and CURVED reduced models are not supported.  "
                "(got dim=%s, vtypes=%s, curve_order=%r)."
                % (mesh.dim, sorted(_ivt) if _ivt else None, curve_order))
    if mesh.dim == 2:
        # 2D PLANAR (motor cross-section) layer: tri/quad cells + boundary-edge charges, log kernel.
        return _build_charge_gram_2d(fes, eps=eps, leafsize=leafsize, eta=eta)
    _vtypes = set(len(el.vertices) for el in mesh.Elements(ng.VOL))
    if _vtypes == {8}:
        # PURE-HEX RT1: the hex-mode charge Gram (Q1 volume charge + Q2 geometry; FLAT or Curve(2) one path).
        # curve_order is IGNORED for hex -- curved is automatic (GetTrafo picks up mesh.Curve(2)); the caller
        # Curve(2)'s the mesh before this call, exactly like the tet curved path.  Uses the hex-gated params.
        return _build_charge_gram_hex(
            fes, eps=eps, leafsize=leafsize, eta=eta,
            image_masks=image_masks, image_signs=image_signs)
    if _vtypes == {6}:
        # PURE-WEDGE (PRISM) RT1: the wedge-mode charge Gram (6-monomial volume charge + mixed tri/quad-face
        # surface charge; 18-node Q2 geometry; FLAT or Curve(2) one path).  curve_order is IGNORED (curved is
        # automatic via GetTrafo picking up mesh.Curve(2)), same as the hex path.
        return _build_charge_gram_wedge(
            fes, eps=eps, leafsize=leafsize, eta=eta,
            image_masks=image_masks, image_signs=image_signs)
    if _vtypes != {4}:
        raise ValueError(
            "vim.ChargeGram: HDiv-VIM is TET (tri-face), pure-HEX (quad-face), or pure-WEDGE/prism "
            "(6-vertex) -- a MIXED-element mesh (e.g. tet+hex) needs HDiv-pyramid transition elements "
            "(NGSolve 6.2.2604 does NOT implement them yet).  Got vertex counts %s." % sorted(_vtypes))
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
    # DEFAULT depends on the USE (2026-06-30).  The LINEAR demag only needs the 3*p PSD FLOOR (RT1 -> quad=3):
    # that makes the NEAR build ~1.8-2.1x cheaper (the near U-list is 98% of the build = the dominant lever) and
    # is VALIDATED to preserve demag (7e-6), per-element leak + magnetic moment, and PSD (min eig ~0).  The
    # NONLINEAR energy-Newton KEEPS the +1 margin (max(3*p,4) -> quad=4 for RT1): the energy HESSIAN wants more
    # than the linear near accuracy at DEEP saturation -- with the OLD product _tet_ref(3) (effective degree 3)
    # deep saturation still converged to M->Msat but took ~2x more Newton iters (195 vs <100).  quad=4 buys the
    # effective degree 5 the energy solve needs.  NOTE (2026-07-02): both quad==3 and quad==4 now route through
    # the SYMMETRIC degree-5 Keast-15/Dunavant-7 outer rule (_outer_tet/_outer_tri) -- the product _tet_ref(4)
    # (64 pts) is itself only effective degree 5 (the Duffy Jacobian lowers 2o-1 to 2o-3), so the 15-pt
    # symmetric rule matches its degree at 4.3x fewer points and the nonlinear solve converges in the SAME or
    # FEWER iters (deep saturation 15 vs 64-pt 17; golden test_hdiv_vim_energy_newton).  So quad=4 is still the
    # nonlinear choice (degree 5 > the old degree-3 product-27), just delivered by the symmetric rule.  intorder
    # still overrides for an explicit choice.  See tests/feec/test_hdiv_vim_psd.py + _symmetric_outer_quad.py.
    quad = (max(p + 2, (intorder + 1) // 2) if intorder is not None
            else (max(3 * p, 4) if nonlinear else 3 * p))
    if curve_order is not None:
        # CURVED (isoparametric P2) Gram: curved charge map B (reference-frame change-of-basis) + the C++
        # curved-Duffy charge Gram.  Only P2 is wired (the C++ CurvedTet/TriPotential are P2); the mesh must
        # already be mesh.Curve(2)'d (the caller does it, per the NGSolve convention).
        if int(curve_order) != 2:
            raise NotImplementedError("vim.ChargeGram: only curve_order=2 (isoparametric P2) is wired "
                                      "(the C++ CurvedTet/TriPotential are P2); got %r." % (curve_order,))
        cbk = _charge_basis_curved(fes, quad)
        rtp, rtw = _tet_ref(quad); rsp, rsw = _tri_ref(quad)
        gx, gw = _g01(int(curve_gauss))
        G = _rp._ChargeGramHMatrix(
            cell_nodes=cbk["cell_nodes"], face_nodes=cbk["face_nodes"], n_el=cbk["n_el"], curve_order=2,
            charge_host=cbk["host"], charge_kind=cbk["kind"], charge_expo=cbk["expo"],
            ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
            ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(),
            curve_gl=gx.tolist(), curve_gw=gw.tolist(), eps=eps, leaf=leafsize, eta=eta)
        return cbk["B"], G, cbk["M_mass"]
    # INNER subtraction quad (B2 speedup): the subtraction remainder (m_src(y)-m_src(p)) is SMOOTH (the
    # singular part is carried EXACTLY by the analytic PhiTet/TriPotential base), so the inner sum uses a
    # COARSER rule than the outer -> another ~1.5-2x on the O(quad_out^3 * quad_in^3) near entries.  Floor at
    # max(quad-2, p+1); only passed to C++ when iq < quad (else inner = outer).  Validated to hold the same
    # demag accuracy as inner=outer by the nearfar/operator goldens + the uniform-1/3 metric.
    iq = inner_quad if inner_quad is not None else max(quad - 2, p + 1)
    cb = _charge_basis(fes, quad)
    B, M_mass, host, kind, expo = cb["B"], cb["M_mass"], cb["host"], cb["kind"], cb["expo"]
    cell_verts, face_verts, n_el = cb["cell_verts"], cb["face_verts"], cb["n_el"]
    # OUTER Gram quadrature: symmetric degree-5 (Keast-15/Dunavant-7) at quad in {3,4}; else product.
    rtp, rtw = _outer_tet(quad)
    rsp, rsw = _outer_tri(quad)
    if p == 0:
        # order-0 = CONSTANT charges -> use the FAST ANALYTIC Gram (Wilton/PhiTet, exact), NOT quadrature.
        # The high-order QUADRATURE constructor (Sauter-Schwab over charge pairs) is ~100x slower per
        # H-matrix entry and is pure waste for constant charges -- it was THE DemagOperator build bottleneck
        # (>195s @ 5310 tets; the change-of-basis is only ~4s).  Charge order [cells..., faces...] matches B's
        # row order, so this is the same B^T G B as the validated solve_nonlinear_newton_scalable order-0 path.
        G = _rp._ChargeGramHMatrix(cell_verts=cell_verts, face_verts=face_verts, n_el=n_el,
                                   eps=eps, leaf=leafsize, eta=eta, near_factor=1e30)
    else:
        # order p>0: monomial-charge quadrature Gram.  By DEFAULT (ho_far_factor=2.0) the near/far split is
        # ON -- well-separated pairs use the cheap LOW-quad QuadDotFar, near/self keep the full high-quad
        # subtraction (validated accuracy-preserving by the nearfar golden).  Pass ho_far_factor=inf to
        # DISABLE the split (every pair high-quad = the exact reference build; the _lo rules are then not
        # built or passed, binding the same way the pre-far-split overload did).
        kw = dict(cell_verts=cell_verts, face_verts=face_verts,
                  n_el=n_el, charge_host=host, charge_kind=kind, charge_expo=expo,
                  ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
                  ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(),
                  eps=eps, leaf=leafsize, eta=eta)
        if np.isfinite(ho_far_factor):
            rtp_lo, rtw_lo = _outer_tet(far_quad)      # symmetric degree-5 at far_quad in {3,4}, else product
            rsp_lo, rsw_lo = _outer_tri(far_quad)
            kw.update(ref_tet_pts_lo=rtp_lo.ravel().tolist(), ref_tet_w_lo=rtw_lo.tolist(),
                      ref_tri_pts_lo=rsp_lo.ravel().tolist(), ref_tri_w_lo=rsw_lo.tolist(),
                      ho_far_factor=ho_far_factor)
        if iq < quad:
            rtp_in, rtw_in = _tet_ref(iq)
            rsp_in, rsw_in = _tri_ref(iq)
            kw.update(ref_tet_pts_in=rtp_in.ravel().tolist(), ref_tet_w_in=rtw_in.tolist(),
                      ref_tri_pts_in=rsp_in.ravel().tolist(), ref_tri_w_in=rsw_in.tolist())
        if image_masks:                                    # IMA (flat pure-tet): fold mirror-image charges
            kw.update(image_masks=list(image_masks), image_signs=list(image_signs))
        G = _rp._ChargeGramHMatrix(**kw)
    return B, G, M_mass


def _gauss_point_cloud(cb, qpts):
    """High-order Gauss point cloud from the shared charge basis cb (=_charge_basis output).  Places
    qpts/dim Gauss-Duffy points per HOST element; the scatter coefficient of charge a (host h, expo) onto
    its host point p is coef[a][p] = W_p * lam_p^expo (W_p the physical quad weight, lam_p the ref point).

    Returns (point_coords [n_point*3], P_pt, P_chg, P_coef  -- the COO scatter for _ChargeGaussHMatrix),
    plus per-charge (centroid, size) for the near criterion and per-charge (point idx, coef) arrays for the
    point-quadrature near reference."""
    host, kind, expo = cb["host"], cb["kind"], cb["expo"]
    vV, bV = cb["vV"], cb["bV"]
    n_charge = len(host)
    rtp, rtw = _tet_ref(qpts)
    rsp, rsw = _tri_ref(qpts)
    # per (kind,host) point block, built once and shared by all that host's monomial charges
    pt_coords, host_block = [], {}      # (kind,h) -> (base, lam, W)
    for a in range(n_charge):
        key = (kind[a], host[a])
        if key in host_block:
            continue
        V = (vV if kind[a] == 0 else bV)[host[a]]
        if kind[a] == 0:
            lam, rw = rtp, rtw
            X = V[0] + lam @ (V[1:] - V[0])
            W = rw * 6.0 * abs(np.linalg.det(V[1:] - V[0])) / 6.0
        else:
            lam, rw = rsp, rsw
            X = V[0] + lam @ (V[1:] - V[0])
            W = rw * 2.0 * (0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0])))
        host_block[key] = (len(pt_coords) // 3, lam, W)
        pt_coords.extend(X.ravel().tolist())
    point_coords = np.asarray(pt_coords, float).reshape(-1, 3)
    P_pt, P_chg, P_coef = [], [], []
    chg_pts, chg_coef, cent, size = [], [], np.zeros((n_charge, 3)), np.zeros(n_charge)
    for a in range(n_charge):
        base, lam, W = host_block[(kind[a], host[a])]
        e = expo[3 * a:3 * a + 3]
        if kind[a] == 0:
            mon = lam[:, 0] ** e[0] * lam[:, 1] ** e[1] * lam[:, 2] ** e[2]
        else:
            mon = lam[:, 0] ** e[0] * lam[:, 1] ** e[1]
        coef = W * mon
        idx = np.arange(base, base + len(lam))
        P_pt.extend(idx.tolist()); P_chg.extend([a] * len(lam)); P_coef.extend(coef.tolist())
        chg_pts.append(idx); chg_coef.append(coef)
        Xh = point_coords[idx]
        cent[a] = Xh.mean(0); size[a] = float(np.max(np.linalg.norm(Xh - cent[a], axis=1)))
    return dict(point_coords=point_coords, P_pt=P_pt, P_chg=P_chg, P_coef=P_coef,
                chg_pts=chg_pts, chg_coef=chg_coef, cent=cent, size=size)


def build_charge_gauss(fes, qpts=3, near_factor=1.0, eps=1e-5, leafsize=64, eta=2.0):
    """High-order HDiv-VIM charge Gram as a GAUSS POINT operator (G ~= P^T K_point P + sparse near
    correction) -- the scalable alternative to vim.ChargeGram's analytic Gram for high order, where the
    analytic per-pair entry cost O((3p)^6) explodes.  Returns (B, G_gauss, M_mass) with the SAME B / M_mass
    as vim.ChargeGram (so N = B^T G_gauss B is the same demag operator, validated to ~1e-4 at p<=2).

    The point H-matrix carries the FAR field (cheap 1/r); the sparse near correction (pairs within
    near_factor*(size_a+size_b)) restores the exact analytic entry via a build=False oracle.  qpts = the
    Gauss-Duffy points/dim of the point cloud (qpts=3 validated ~1e-4 at near_factor=1.0); scope: tet meshes.
    CALLER wraps in TaskManager."""
    from scipy.spatial import cKDTree
    raise NotImplementedError(
        "vim.ChargeGramGauss is RETIRED: the Gauss point-operator charge Gram was an RT0 build-speed "
        "experiment.  HDiv-VIM (RT1, tet) uses the analytic charge Gram (vim.ChargeGram).")
    p = fes.globalorder
    quad = max(3 * p, 4)
    cb = _charge_basis(fes, quad)
    B, M_mass = cb["B"], cb["M_mass"]
    n_charge = len(cb["host"])
    pc = _gauss_point_cloud(cb, qpts)
    cent, size = pc["cent"], pc["size"]
    chg_pts, chg_coef = pc["chg_pts"], pc["chg_coef"]
    X = pc["point_coords"]
    inv4pi = 1.0 / (4.0 * np.pi)

    # exact analytic ENTRY ORACLE (build=False -> no full H-matrix build): p==0 analytic, p>0 high-order.
    if p == 0:
        oracle = _rp._ChargeGramHMatrix(cell_verts=cb["cell_verts"], face_verts=cb["face_verts"],
                                        n_el=cb["n_el"], near_factor=1e30, build=False)
    else:
        rtp, rtw = _tet_ref(quad); rsp, rsw = _tri_ref(quad)
        oracle = _rp._ChargeGramHMatrix(
            cell_verts=cb["cell_verts"], face_verts=cb["face_verts"], n_el=cb["n_el"],
            charge_host=cb["host"], charge_kind=cb["kind"], charge_expo=cb["expo"],
            ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
            ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(), build=False)

    def point_entry(a, b):
        ia, ib = chg_pts[a], chg_pts[b]
        Xa, Xb = X[ia], X[ib]
        D = np.linalg.norm(Xa[:, None, :] - Xb[None, :, :], axis=2)
        with np.errstate(divide="ignore"):
            K = np.where(D > 1e-300, inv4pi / D, 0.0)
        if a == b:
            np.fill_diagonal(K, 0.0)               # self pair excludes coincident points (carried by oracle)
        return float(chg_coef[a] @ (K @ chg_coef[b]))

    # near pairs (O(N)): exact analytic - point quadrature.  Same KDTree neighbourhood as _solve.
    tree = cKDTree(cent)
    max_size = float(np.max(size)) if n_charge else 0.0
    corr_i, corr_j, corr_v = [], [], []
    for a in range(n_charge):
        for b in tree.query_ball_point(cent[a], near_factor * (size[a] + max_size)):
            if b < a:
                continue
            r = float(np.linalg.norm(cent[a] - cent[b]))
            if a == b or r <= near_factor * (size[a] + size[b]):
                delta = float(oracle.entry(int(a), int(b))) - point_entry(a, b)
                corr_i.append(int(a)); corr_j.append(int(b)); corr_v.append(delta)
                if a != b:
                    corr_i.append(int(b)); corr_j.append(int(a)); corr_v.append(delta)
    G = _rp._ChargeGaussHMatrix(
        point_coords=X.ravel().tolist(), P_pt=list(map(int, pc["P_pt"])),
        P_chg=list(map(int, pc["P_chg"])), P_coef=list(map(float, pc["P_coef"])),
        n_charge=n_charge, corr_i=corr_i, corr_j=corr_j, corr_v=corr_v,
        eps=eps, leaf=leafsize, eta=eta)
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
    wraps construction + DemagFactor in `with TaskManager():`.

    gram_backend selects the charge-Gram H-matrix:
      "analytic" (default) -- the exact analytic charge Gram (vim.ChargeGram); the per-pair entry cost is
        O((3p)^6) at order p (the singular outer x inner subtraction quadrature), so the BUILD explodes with p.
      "gauss" -- the Gauss POINT operator (vim.ChargeGramGauss): G ~= P^T K_point P + sparse near correction,
        the cheap-1/r point H-matrix carrying the far field, the analytic entry used ONLY on the O(N) near
        pairs.  Validated to match "analytic" demag to ~1e-4 at p<=2 (qpts/gauss_near_factor control the
        accuracy).  Aimed at the high-order / curved regime where the analytic build is the bottleneck.
    Both backends produce the SAME B / M_mass and a backend-agnostic `.mat` (N = B^T G B); DemagFactor and
    the NGSolve-composability are identical."""

    def __init__(self, fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0,
                 far_quad=3, ho_far_factor=2.0, inner_quad=None,
                 gram_backend="analytic", qpts=3, gauss_near_factor=1.0,
                 curve_order=None, curve_gauss=8):
        if gram_backend != "analytic":
            raise ValueError(
                "DemagOperator: gram_backend must be 'analytic' -- the 'gauss' point-operator backend is "
                "retired (RT0 build-speed experiment).  (got %r)" % (gram_backend,))
        if fes.globalorder != 1:
            raise ValueError(
                "DemagOperator: HDiv-VIM is RT1 (HDiv order=1) only -- RT0 (order=0) is retired (per-element "
                "inaccurate) and RT2+ is retired (no gain over RT1).  Build the FESpace as HDiv(mesh, order=1).")
        self.space = fes
        self.gram_backend = gram_backend
        # AUTO-MATCH the Gram curve order to the MESH geometry order (mesh.GetCurveOrder()).  A STRAIGHT Gram on
        # a CURVED mesh (where B/M_mass are NGSolve curved integrals) is geometry-inconsistent and the demag
        # factor DRIFTS with geometry order (sphere: straight-Gram 0.336/0.308/0.279 at curve 1/2/3; the matched
        # curved Gram restores ~1/3 -- 0.338 at curve 2).  curve_order=None => auto from GetCurveOrder(); pass an
        # explicit int to override (curve_order=0 forces the STRAIGHT Gram, e.g. a deliberate flat-Gram probe).
        if curve_order is None and gram_backend == "analytic":
            _k = fes.mesh.GetCurveOrder()
            curve_order = _k if _k >= 2 else None
        elif curve_order == 0:
            curve_order = None
        self.curve_order = curve_order
        if gram_backend == "gauss":
            self._B, self._G, self._Mmass = build_charge_gauss(
                fes, qpts=qpts, near_factor=gauss_near_factor, eps=eps, leafsize=leafsize, eta=eta)
        else:
            self._B, self._G, self._Mmass = build_charge_gram(
                fes, intorder=intorder, eps=eps, leafsize=leafsize, eta=eta,
                far_quad=far_quad, ho_far_factor=ho_far_factor, inner_quad=inner_quad,
                curve_order=curve_order, curve_gauss=curve_gauss)
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
