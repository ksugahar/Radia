"""radia.vim._vim -- an ngsolve.bem-STYLE API for the HDiv-type VIM demag operator.

Mirrors the ngsolve.bem design (SingleLayerPotentialOperator etc.): construct the operator from an NGSolve
FESpace -- the production space is `HDiv(mesh, order=1)` -- and expose
`.mat`, an H-matrix-backed NGSolve `BaseMatrix` that composes with NGSolve's solvers / BlockMatrix just
like `SingleLayerPotentialOperator(fes, ...).mat`:

    from ngsolve import *
    from ngsolve.krylovspace import GMRes
    from radia.vim import DemagOperator

    mesh = Mesh(...)
    fes  = HDiv(mesh, order=1)
    with TaskManager():
        N = DemagOperator(fes, eps=1e-7)             # like SingleLayerPotentialOperator(fes, ...)
        # N.mat : BaseMatrix == the demag operator B^T G B (G = the C++ charge-Gram H-matrix)
        u, v = fes.TnT()
        M = BilinearForm(u*v*dx).Assemble()          # the HDiv mass
        A = (1.0/chi)*M.mat + N.mat                  # physical SPD +N system
        gfm = GridFunction(fes)
        gfm.vec.data = GMRes(A=A, b=rhs.vec, tol=1e-8, maxsteps=400)

Convenience: ``N.DemagFactor(CF((0,0,1)))`` -> the Rayleigh quotient (the demag factor, ~1/3 for a sphere).

Backend: the C++ charge-Gram H-matrix (radia._radia_pybind._ChargeGramHMatrix).  The charge basis is
element-local monomials (host reference coordinates); N = B^T G B is basis-invariant, so the demag matches
the NGSolve-L2-basis dense reference.  Python declares NGSolve spaces/forms and prepares the one-time sparse
charge topology.  The assembled charge map, mass matrices, persistent NGSolve ``BaseMatrix``, Krylov solve,
and field source/evaluation live in C++; no dense O(N^2) operator is formed and no Python list is used on the
solve/field hot path.

TaskManager: per the caller-wraps policy, this module does NOT open a TaskManager; the CALLER wraps the
DemagOperator construction + DemagFactor / solve in `with TaskManager():` (the ngsolve.bem idiom).
"""
import numpy as np
import scipy.sparse as sp

from ._capabilities import validate_hdiv_configuration
import ngsolve as ng
import time

import radia._radia_pybind as _rp


def _f64_buffer(value):
    """Return the canonical contiguous float64 buffer for the C++ boundary."""
    return np.ascontiguousarray(value, dtype=np.float64).reshape(-1)


def _i32_buffer(value):
    """Return the canonical contiguous int32 buffer for the C++ boundary."""
    return np.ascontiguousarray(value, dtype=np.int32).reshape(-1)


_EMPTY_F64 = np.empty(0, dtype=np.float64)
_EMPTY_I32 = np.empty(0, dtype=np.int32)


def _volume_vertex_counts(mesh):
    """Unique volume-element vertex counts from the native MeshAccess."""
    return frozenset(int(value) for value in _rp._volume_element_vertex_counts(mesh))


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


# Fully symmetric degree-10 simplex rules (Witherden--Vincent 2015).  Curved
# BDM2 uses these instead of a vertex-order-dependent collapsed product rule:
# 81 tet / 25 tri points versus 216 / 36 at product order 6.  Orbit expansion
# keeps the point set exactly closed under every reference-vertex permutation.
def _tet_ref_sym10():
    orbits = [
        ((0.25, 0.25, 0.25, 0.25), 0.04739977355602074),
    ]
    for weight, a in zip(
            (0.0269370599922687, 0.009869159716793382),
            (0.3122500686951887, 0.1143096538573461)):
        orbits.append(((a, a, a, 1.0 - 3.0*a), weight))
    weights = (0.01139388122019523, 0.0003619443443392536,
               0.02573973198045607, 0.01013587167975579,
               0.006576147277035904, 0.01290703579886199)
    aa = (0.4104307392189654, 0.006138008824790653,
          0.1210501811455894, 0.03277946821644262,
          0.03248528156482305, 0.174979342183939)
    bb = (0.1654860256196111, 0.9429887673452049,
          0.4771903799042804, 0.594256269480007,
          0.8011772846583444, 0.628071845475366)
    for weight, a, b in zip(weights, aa, bb):
        orbits.append(((a, a, b, 1.0 - 2.0*a - b), weight))
    points, rule_weights = [], []
    for bary, weight in orbits:
        for point in _sym_orbit(bary, 3):
            points.append(point)
            rule_weights.append(weight / 6.0)
    return np.array(points), np.array(rule_weights)


def _tri_ref_sym10():
    orbits = [
        ((1.0/3.0, 1.0/3.0, 1.0/3.0), 0.081743329146285973),
    ]
    for weight, a in zip(
            (0.013352968813149567, 0.045957963604744731),
            (0.032055373216943517, 0.14216110105656438)):
        orbits.append(((a, a, 1.0 - 2.0*a), weight))
    weights = (0.025297757707288385, 0.034184648162959429,
               0.063904906396424044)
    aa = (0.028367665339938453, 0.029619889488729734,
          0.14813288578382056)
    bb = (0.1637017337371825, 0.36914678182781102,
          0.32181299528883545)
    for weight, a, b in zip(weights, aa, bb):
        orbits.append(((a, b, 1.0 - a - b), weight))
    points, rule_weights = [], []
    for bary, weight in orbits:
        for point in _sym_orbit(bary, 2):
            points.append(point)
            rule_weights.append(weight * 0.5)
    return np.array(points), np.array(rule_weights)


_SYM10_TET = _tet_ref_sym10()
_SYM10_TRI = _tri_ref_sym10()


def _curved_outer_rules(order):
    """Permutation-invariant curved Gram rules for production BDM1/BDM2."""
    if int(order) == 1:
        return _SYM5_TET, _SYM5_TRI
    if int(order) == 2:
        return _SYM10_TET, _SYM10_TRI
    raise ValueError("curved HDiv production order must be 1 or 2")


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
    a few IntegrationRule evaluations.  pt is in NGSolve's reference frame.

    NGSolve's scalar ``ElementTransformation`` path can return a transient,
    uninitialised point on its first touch after a different element topology
    was used in the same TaskManager.  Accept two consecutive finite,
    bit-identical fits only; the retry is cheap compared with the subsequent
    per-element change-of-basis and prevents a bad affine map from poisoning
    the HDiv charge basis.
    """
    ir = ng.IntegrationRule(eltype, 2)
    rp = np.array([list(p.point)[:ndim] for p in ir])
    A = np.hstack([rp, np.ones((len(rp), 1))])
    previous = None
    for attempt in range(16):
        Pp = np.array([list(trafo(p).point) for p in ir])
        X, *_ = np.linalg.lstsq(A, Pp, rcond=None)           # [Jng^T (ndim rows) ; P0]
        fitted = np.hstack([X[ndim, :], X[:ndim, :].T.ravel()])
        if (previous is not None and np.all(np.isfinite(fitted))
                and np.array_equal(previous, fitted)):
            return X[ndim, :], X[:ndim, :].T                  # (P0, Jng 3xndim)
        previous = fitted
        if attempt:
            time.sleep(0.002)
    raise RuntimeError(
        "NGSolve affine ElementTransformation remained unstable after 16 "
        "finite/consecutive retries; refusing to build a non-deterministic "
        "HDiv charge basis.")


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


def _exterior_bnd_elements(mesh):
    """Surface elements on the TRUE exterior boundary: their facet is owned by exactly ONE volume element.

    The charge layer treats every listed surface element as a charged boundary face (sigma = M.n,
    single-sided).  That is correct only for the exterior skin.  An INTERNAL surface element (a conforming
    multi-region material interface from glued CAD, fd with domin != 0 and domout != 0) sits between two
    solved volume elements whose normal trace is HDiv-continuous, so its true single-layer charge is the
    JUMP = 0; charging it single-sided injects a spurious charge sheet.  Measured 2026-07-28 (two-region
    glued ball, conforming, uniform chi=100): <Mz> reads -12 % low with the interface flux clamped ~23x;
    dropping the internal faces restores the exact single-region result (2.9190/2.9134/2.9125 at maxh
    .35/.22/.14 vs exact 2.91262) and makes conforming multi-region meshes valid inputs.

    Matching is by sorted corner-vertex tuples: BND-element facet NUMBERS are not comparable with the
    volume facet numbering, and netgen fd/index conventions are not relied on.  Works for tri/quad faces
    (tet/hex/wedge) and for boundary edges of planar 2D meshes.  Cost is one O(n_el) dict pass, negligible
    vs the Gram build.  A surface element matching NO volume facet means a non-conforming or corrupt mesh:
    raise (fail loud), never charge it."""
    owners = {}
    for el in mesh.Elements(ng.VOL):
        for fa in el.facets:
            key = tuple(sorted(v.nr for v in mesh[fa].vertices))
            owners[key] = owners.get(key, 0) + 1
    kept = []
    for i in range(mesh.GetNE(ng.BND)):
        e = ng.ElementId(ng.BND, i)
        n = owners.get(tuple(sorted(v.nr for v in mesh[e].vertices)), 0)
        if n == 1:
            kept.append(e)
        elif n != 2:
            raise ValueError(
                "HDiv-VIM charge basis: surface element %d matches %d volume facets "
                "(expected 1 = exterior, or 2 = internal interface which carries no "
                "single-sided charge); the mesh is non-conforming or corrupt." % (i, n))
    return kept


def _assert_broken_hdiv(fes):
    """Require element-local HDiv unknowns for explicit interface charges.

    An internal charge sheet is the jump of ``M.n`` between the two element
    sides.  A conforming HDiv space identifies those trace DoFs, so asking for
    explicit internal charges on it would silently produce zero.  Detect that
    invalid combination from the actual element DoF ownership instead of an
    undocumented NGSolve flag.
    """
    owners = {}
    for el in fes.mesh.Elements(ng.VOL):
        for facet in el.facets:
            owners.setdefault(int(facet.nr), []).append(el)
    for adjacent in owners.values():
        if len(adjacent) == 2:
            left = set(int(d) for d in fes.GetDofNrs(adjacent[0]) if int(d) >= 0)
            right = set(int(d) for d in fes.GetDofNrs(adjacent[1]) if int(d) >= 0)
            if left & right:
                raise ValueError(
                    "vim.ChargeGram: internal_interfaces=True requires "
                    "ng.HDiv(..., discontinuous=True); the supplied space "
                    "shares normal-trace DoFs across an internal facet")
            return
    if fes.mesh.ne > 1:
        raise ValueError(
            "vim.ChargeGram: could not find a two-sided internal facet while "
            "checking internal_interfaces=True")


def _broken_tet_face_charge_basis(fes, p):
    """NGSolve-assembled jump ``[M.n]`` on every straight triangular facet.

    The trial space is broken HDiv and the test space is NGSolve's shared
    ``FacetFESpace``.  ``dx(element_boundary=True)`` visits both sides of an
    internal facet with their outward normals; assembly into the shared facet
    test DoFs therefore forms the signed jump without reconstructing an HDiv
    basis in Python.  The small per-facet moment conversion only changes the
    NGSolve facet polynomial into the C++ triangle's monomial frame.
    """
    _assert_broken_hdiv(fes)
    mesh = fes.mesh
    facet_space = ng.FacetFESpace(mesh, order=int(p))
    u = fes.TrialFunction()
    q = facet_space.TestFunction()
    normal = ng.specialcf.normal(mesh.dim)
    jump = ng.BilinearForm(trialspace=fes, testspace=facet_space)
    jump += (u * normal) * q * ng.dx(element_boundary=True)
    jump.Assemble()
    jump_moments = _csr(jump)

    # Physical monomial moments are assembled by NGSolve.  Restricting these
    # global polynomials to one affine triangle and changing to the C++ local
    # (xi,eta) monomials is geometry algebra, not an FE-basis reimplementation.
    xyz_mons = _monos_vol(int(p))
    moment_vectors = []
    for i, j, k in xyz_mons:
        physical_monomial = ng.x ** i * ng.y ** j * ng.z ** k
        linear = ng.LinearForm(facet_space)
        linear += physical_monomial * q * ng.dx(element_boundary=True)
        linear.Assemble()
        moment_vectors.append(
            np.asarray(linear.vec.FV().NumPy(), dtype=float).copy())

    owner_count = {int(facet.nr): 0 for facet in mesh.facets}
    for el in mesh.Elements(ng.VOL):
        for facet in el.facets:
            owner_count[int(facet.nr)] += 1

    facets = list(mesh.facets)
    mons = _monos_surf(int(p))
    ref_points, _ = _tri_ref(max(int(p) + 1, 2))
    Q = np.asarray([[x ** i * y ** j for i, j in mons]
                    for x, y in ref_points], dtype=float)
    rows = []
    face_vertices = []
    for facet in facets:
        nr = int(facet.nr)
        dofs = [int(d) for d in facet_space.GetDofNrs(facet)]
        if len(dofs) != len(mons):
            raise RuntimeError(
                "vim.ChargeGram: FacetFESpace order %d has %d DoFs on "
                "triangle %d, expected %d" % (p, len(dofs), nr, len(mons)))
        owners = owner_count[nr]
        if owners not in (1, 2):
            raise ValueError(
                "vim.ChargeGram: facet %d has %d volume owners (expected "
                "1 exterior or 2 internal)" % (nr, owners))
        vertices = np.asarray([mesh[v].point for v in facet.vertices], dtype=float)
        if vertices.shape != (3, 3):
            raise ValueError(
                "vim.ChargeGram: internal TET interface path requires "
                "triangular facets")
        face_vertices.append(vertices)

        # Linear-form assembly visits an internal facet twice.  Divide those
        # test moments by the owner count to recover one physical-face moment;
        # do NOT divide jump_moments, whose two signed contributions are the
        # desired jump.
        D = np.column_stack([values[dofs] / owners
                             for values in moment_vectors])
        physical_points = (vertices[0][None, :]
                           + ref_points[:, 0, None]
                           * (vertices[1] - vertices[0])[None, :]
                           + ref_points[:, 1, None]
                           * (vertices[2] - vertices[0])[None, :])
        P = np.asarray([[x ** i * y ** j * z ** k for i, j, k in xyz_mons]
                        for x, y, z in physical_points], dtype=float)
        # P = Q @ local_to_physical.  Since the global monomial family is
        # redundant on a plane, recover the unique NGSolve-test/local-monomial
        # moment matrix with the Moore-Penrose right inverse.
        local_to_physical = np.linalg.lstsq(Q, P, rcond=None)[0]
        C = D @ np.linalg.pinv(local_to_physical)
        if np.linalg.cond(C) > 1e12:
            raise RuntimeError(
                "vim.ChargeGram: ill-conditioned facet moment conversion "
                "on facet %d (cond=%g)" % (nr, np.linalg.cond(C)))
        rows.append(sp.csr_matrix(np.linalg.solve(C, jump_moments[dofs, :].toarray())))

    return dict(
        B=sp.vstack(rows).tocsr(),
        face_vertices=face_vertices,
        mons=mons,
    )


def _charge_basis(fes, quad, *, materialize_mass=True,
                  internal_interfaces=False):
    """Shared geometry + monomial charge-density map for the BDM1 HDiv-VIM operator.

    Returns B (CSR
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
    M_mass = _csr(mh) if materialize_mass else None

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = (_exterior_bnd_elements(mesh) if not internal_interfaces
            else [])
    # A linear topological TET mesh can carry a live NGSolve deformation.
    # ``mesh[v].point`` remains the undeformed topological coordinate in that
    # state, whereas every assembled FE form above uses ``GetTrafo``.  Mixing
    # those geometries makes the rebuilt ChargeGram inconsistent with its
    # mass/RHS and invalidates shape-derivative regressions.  Read the four/
    # three reference vertices through the mapped-rule path whenever a live
    # deformation is installed.  The raw vertex path remains the cheap source
    # of truth for an undeformed affine mesh.
    live_deformation = mesh.deformation is not None
    vV = ([_trafo_lattice_nodes(mesh, e, _IR_TET_VERTICES)[
               _TET_REF_TO_ELEMENT_VERTEX_ORDER]
           for e in vels]
          if live_deformation else
          [np.array([mesh[v].point for v in mesh[e].vertices]) for e in vels])
    bV = ([_trafo_lattice_nodes(mesh, e, _IR_TRI_VERTICES)[
               _TRI_REF_TO_ELEMENT_VERTEX_ORDER]
           for e in bels]
          if live_deformation else
          [np.array([mesh[v].point for v in mesh[e].vertices]) for e in bels])
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
        # BDM1 fast path: the VOLUME change-of-basis is the IDENTITY (Sv == [[1.0]], geometry-INDEPENDENT -- the
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
    if internal_interfaces:
        broken_faces = _broken_tet_face_charge_basis(fes, p)
        bV = broken_faces["face_vertices"]
        Bface = broken_faces["B"]
        for f in range(len(bV)):
            for i, j in broken_faces["mons"]:
                host.append(f); kind.append(1); expo += [i, j, 0]
        B = sp.vstack([sp.vstack(Brows).tocsr(), Bface]).tocsr()
    else:
        for f in range(len(bels)):
            Ss = sp.csr_matrix(_change_of_basis(L2b.GetFE(bels[f]), mons_s, *rsq, dim=2,
                                                trafo=mesh.GetTrafo(bels[f]), Vmesh=bV[f]))
            blk = Ss @ Bb_d[bdof[f], :]                         # sparse (nmons_s x ndof)
            for a, (i, j) in enumerate(mons_s):
                Brows.append(blk[a]); host.append(f); kind.append(1); expo += [i, j, 0]
        B = sp.vstack(Brows).tocsr()                            # (n_charge, ndof)
    return dict(B=B, M_mass=M_mass, M_mass_ngsolve=mh.mat,
                host=host, kind=kind, expo=expo, vV=vV, bV=bV,
                cell_verts=_f64_buffer(np.concatenate([V.ravel() for V in vV])),
                face_verts=_f64_buffer(np.concatenate([V.ravel() for V in bV])),
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
_IR_TET_VERTICES = ng.IntegrationRule(
    [tuple(r) for r in _TET_REFNODES[:4]], [1.0] * 4
)
_IR_TRI_VERTICES = ng.IntegrationRule(
    [(r[0], r[1]) for r in _TRI_REFNODES[:3]], [1.0] * 3
)
# NGSolve's simplex reference convention starts at the final topological
# element vertex.  Reorder mapped reference corners back to
# ``mesh[e].vertices`` order, which is the established C++ flat-simplex
# geometry contract used by ``_change_of_basis`` and ChargeGram.
_TET_REF_TO_ELEMENT_VERTEX_ORDER = np.array([1, 2, 3, 0], dtype=int)
_TRI_REF_TO_ELEMENT_VERTEX_ORDER = np.array([1, 2, 0], dtype=int)


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


def _shape_signature_ref(fe, pts, dim):
    """Small orientation signature for scalar DG FE shapes.

    Computing the full Q1 projection matrix for every hex/quad is pure waste on structured meshes, but blindly
    reusing one transform would be unsafe if NGSolve changes the local scalar-shape ordering.  A few asymmetric
    reference points distinguish the orientation/ordering cheaply; a cache miss then computes the exact
    full-quadrature transform for that signature.
    """
    if dim == 3:
        vals = np.array([fe.CalcShape(float(p[0]), float(p[1]), float(p[2])) for p in pts], dtype=float)
    else:
        vals = np.array([fe.CalcShape(float(p[0]), float(p[1])) for p in pts], dtype=float)
    return fe.ndof, np.round(vals, 14).tobytes()


def _ref_monomial_moment_transform(fe, mons, refP, refW, dim):
    """Direct map from ref-shape moments to ref monomial coefficients.

    Existing code did two dense operations per element:
      shape moments -> shape coefficients via M_shape^{-1}, then
      shape coefficients -> monomial coefficients via S.
    The product S M_shape^{-1} depends only on the scalar DG FE orientation, so the hex path caches this tiny
    matrix and applies it to each element's sparse moment rows.
    """
    if dim == 3:
        Phi = np.array([fe.CalcShape(float(pt[0]), float(pt[1]), float(pt[2])) for pt in refP], dtype=float)
        Mono = np.array([[pt[0] ** i * pt[1] ** j * pt[2] ** k for (i, j, k) in mons]
                         for pt in refP], dtype=float)
    else:
        Phi = np.array([fe.CalcShape(float(pt[0]), float(pt[1])) for pt in refP], dtype=float)
        Mono = np.array([[pt[0] ** i * pt[1] ** j for (i, j) in mons] for pt in refP], dtype=float)
    M_shape = (Phi * refW[:, None]).T @ Phi
    M_mono = (Mono * refW[:, None]).T @ Mono
    cross = (Mono * refW[:, None]).T @ Phi
    shape_to_mono = np.linalg.solve(M_mono, cross)
    return np.linalg.solve(M_shape.T, shape_to_mono.T).T


def _charge_basis_curved(fes, quad, *, materialize_mass=True):
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
    M_mass = _csr(mh) if materialize_mass else None

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = _exterior_bnd_elements(mesh)          # internal material-interface faces carry no charge
    vdof = [list(L2v.GetDofNrs(e)) for e in vels]
    bdof = [list(L2b.GetDofNrs(e)) for e in bels]
    mons_v, mons_s = _monos_vol(pv), _monos_surf(p)
    (rtp, rtw), (rsp, rsw) = _curved_outer_rules(p)

    Brows, host, kind, expo = [], [], [], []
    cell_nodes, face_nodes = [], []
    cell_vertices = [v.nr for e in vels for v in mesh[e].vertices]
    face_vertices = [v.nr for e in bels for v in mesh[e].vertices]
    for c, e in enumerate(vels):
        cell_nodes.append(_trafo_lattice_nodes(mesh, e, _IR_TET_NODES))   # P2 nodes (curved geom, kept)
        if pv == 0:                                            # BDM1: volume Sv == [[1]] (identity) -> Bv_d row direct
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
    # The reference Piola extraction has exact structural zeros.  The dense
    # local solves leave O(1e-14..1e-11) fill in those slots at BDM2; retaining
    # it breaks reflection parity after a high-mu solve.  Snap only values at
    # the roundoff floor relative to this dimensionless reference map.
    if B.nnz:
        zero_floor = 4096*np.finfo(float).eps*np.max(np.abs(B.data))
        B.data[np.abs(B.data) <= zero_floor] = 0.0
        B.eliminate_zeros()
    return dict(B=B, M_mass=M_mass, M_mass_ngsolve=mh.mat,
                host=host, kind=kind, expo=expo,
                cell_nodes=(_f64_buffer(np.concatenate([V.ravel() for V in cell_nodes]))
                            if cell_nodes else _EMPTY_F64),
                face_nodes=(_f64_buffer(np.concatenate([V.ravel() for V in face_nodes]))
                            if face_nodes else _EMPTY_F64),
                cell_vertices=_i32_buffer(cell_vertices), face_vertices=_i32_buffer(face_vertices),
                n_el=len(vels))


# ----------------------------------------------------------- HEX BDM1/BDM2 charge basis ---
# The hex volume charge = -div(HDiv-hex order-p) lives in tensor Qp: 8 modes at BDM1 and 27 at BDM2.
# Testing only the tet-like total-degree space pollutes ker(N), so keep the complete tensor product.
# Geometry = 27-node triquadratic (Q2) lattice via GetTrafo -> ONE code path for FLAT (trilinear subset of Q2,
# exact) AND CURVED (mesh.Curve(2)) hexes.  Face charge is SurfaceL2(order=p), with 4/9 modes per
# quad face at BDM1/BDM2.  The C++
# RadHACApKChargeGram hex mode decomposes each hex into 6 sub-tets (quad face -> 2 sub-tris) as an INTEGRATION
# device and does the both-domains-graded Duffy singular quadrature that keeps eig(M_mass^-1 N) <= 1.  See
# the de-risk + block-memo perf fix in memory hdiv-tet-hex-coupling-pyramid-gated.
def _mons_hex(order):
    return [(i, j, k) for k in range(order+1)
            for j in range(order+1) for i in range(order+1)]


def _mons_quad(order):
    return [(i, j) for j in range(order+1) for i in range(order+1)]
_Q2_LATTICE_3D = [(ix / 2.0, iy / 2.0, iz / 2.0)
                  for iz in range(3) for iy in range(3) for ix in range(3)]  # n = ix + 3*iy + 9*iz
_Q2_LATTICE_2D = [(iu / 2.0, iv / 2.0) for iv in range(3) for iu in range(3)]  # n = iu + 3*iv
_HEX_SHAPE_SIG_PTS = np.array([(0.137, 0.239, 0.361), (0.713, 0.211, 0.557),
                               (0.319, 0.823, 0.173)], dtype=float)
_QUAD_SHAPE_SIG_PTS = np.array([(0.137, 0.239), (0.713, 0.421), (0.319, 0.823)], dtype=float)
_HEX_NGSOLVE_LINEAR_ORDER = np.array([0, 1, 3, 2, 4, 5, 7, 6], dtype=int)
_QUAD_NGSOLVE_LINEAR_ORDER = np.array([0, 1, 3, 2], dtype=int)
_HEX_Q2_LINEAR_WEIGHTS = np.array([
    [(1.0-u)*(1.0-v)*(1.0-w), u*(1.0-v)*(1.0-w), (1.0-u)*v*(1.0-w), u*v*(1.0-w),
     (1.0-u)*(1.0-v)*w,       u*(1.0-v)*w,       (1.0-u)*v*w,       u*v*w]
    for u, v, w in _Q2_LATTICE_3D
], dtype=float)
_QUAD_Q2_LINEAR_WEIGHTS = np.array([
    [(1.0-u)*(1.0-v), u*(1.0-v), (1.0-u)*v, u*v]
    for u, v in _Q2_LATTICE_2D
], dtype=float)


def _mesh_vertices_array(mesh, e):
    return np.array([mesh[v].point for v in mesh[e].vertices], dtype=float)


def _hex_q2_lattice_nodes_ngsolve_linear(mesh, e):
    """27 lattice nodes for a linear NGSolve .vol hex, in the C++ Q2 lattice order.

    This is deliberately NGSolve-reference ordering, not Cubit/GMSH ordering.  Probed against
    `mesh.GetTrafo(e)` on structured .vol hexes:
      local vertices 0,1,3,2,4,5,7,6 map to ref corners
      (0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1).
    For curved .vol meshes this helper is not used; GetTrafo remains the source of truth.
    """
    V = _mesh_vertices_array(mesh, e)[_HEX_NGSOLVE_LINEAR_ORDER]
    return _HEX_Q2_LINEAR_WEIGHTS @ V


def _quad_q2_lattice_nodes_ngsolve_linear(mesh, e):
    """9 lattice nodes for a linear NGSolve .vol boundary quad, in C++ Q2 face order."""
    V = _mesh_vertices_array(mesh, e)[_QUAD_NGSOLVE_LINEAR_ORDER]
    return _QUAD_Q2_LINEAR_WEIGHTS @ V


def _ref_prod_gauss(n, dim):
    """Product Gauss rule on the reference tensor cell for BDM1/BDM2 basis transforms."""
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


def _trafo_lattice_nodes(mesh, e, ir, max_tries=4):
    """Extract physical lattice nodes through NGSolve's vectorized mapped-rule path.

    The scalar ``tr(ip).point`` API can intermittently expose uninitialized coordinates on first touch.
    Evaluating the global coordinate CoefficientFunction on the complete mapped integration rule follows
    the same SIMD geometry path as NGSolve assembly and avoids that failure mode.  Two consecutive,
    bit-identical finite evaluations remain the determinism contract; unstable geometry fails loudly.
    """
    coords = (ng.x, ng.y) if int(mesh.dim) == 2 else (ng.x, ng.y, ng.z)
    coordinate_cf = ng.CF(coords)
    prev = None
    for _ in range(max_tries):
        tr = mesh.GetTrafo(e)
        cur = np.asarray(coordinate_cf(tr(ir)), dtype=float).reshape(-1, len(coords))
        if prev is not None and np.array_equal(prev, cur) and np.all(np.isfinite(cur)):
            return cur
        prev = cur
    raise RuntimeError(
        f"GetTrafo lattice evaluation unstable for element {e} after {max_tries} tries "
        "(NGSolve mapped-rule path returned differing node coordinates -- do not trust this mesh "
        "extraction; abort and report the incident).")


def _broken_hex_face_charge_basis(fes, p):
    """NGSolve-assembled normal jump on every quadrilateral facet.

    ``FacetFESpace`` owns element-side orientation.  The conversion below only
    changes its scalar facet basis into the ``u^i v^j`` tensor monomials used
    by the existing C++ HEX ChargeGram.  No HDiv shape or Piola transform is
    reconstructed in Python.

    For RT0, the reference normal flux is one constant per facet and its
    coefficient is exactly the assembled total flux: physical ``dS`` and the
    Piola surface factor cancel.  Thus mapped bilinear facets need no geometric
    projection.  Higher-order facet polynomials retain the affine-only change
    of basis below; silently projecting them on a warped facet would be wrong.
    """
    _assert_broken_hdiv(fes)
    mesh = fes.mesh
    facet_space = ng.FacetFESpace(mesh, order=int(p))
    u = fes.TrialFunction()
    q = facet_space.TestFunction()
    normal = ng.specialcf.normal(mesh.dim)
    jump = ng.BilinearForm(trialspace=fes, testspace=facet_space)
    jump += (u * normal) * q * ng.dx(element_boundary=True)
    jump.Assemble()
    jump_moments = _csr(jump)

    # On an affine quad, every tensor monomial u^i v^j (i,j <= p) is
    # representable by global monomials of total degree <= 2p.  Assemble the
    # physical moments with NGSolve and solve only the small geometry change of
    # basis per facet.
    xyz_mons = _monos_vol(2 * int(p))
    moment_vectors = []
    for i, j, k in xyz_mons:
        physical_monomial = ng.x ** i * ng.y ** j * ng.z ** k
        linear = ng.LinearForm(facet_space)
        linear += physical_monomial * q * ng.dx(element_boundary=True)
        linear.Assemble()
        moment_vectors.append(
            np.asarray(linear.vec.FV().NumPy(), dtype=float).copy())

    owner_count = {int(facet.nr): 0 for facet in mesh.facets}
    for element in mesh.Elements(ng.VOL):
        for facet in element.facets:
            owner_count[int(facet.nr)] += 1

    mons = _mons_quad(int(p))
    gauss, _ = _g01(max(int(p) + 2, 3))
    ref_points = np.asarray([(uu, vv) for vv in gauss for uu in gauss])
    Q = np.asarray([[uu ** i * vv ** j for i, j in mons]
                    for uu, vv in ref_points], dtype=float)
    rows = []
    face_nodes = []
    for facet in mesh.facets:
        nr = int(facet.nr)
        dofs = [int(d) for d in facet_space.GetDofNrs(facet)]
        if len(dofs) != len(mons):
            raise RuntimeError(
                "vim.ChargeGram: FacetFESpace order %d has %d DoFs on "
                "quadrilateral %d, expected %d"
                % (p, len(dofs), nr, len(mons)))
        owners = owner_count[nr]
        if owners not in (1, 2):
            raise ValueError(
                "vim.ChargeGram: facet %d has %d volume owners (expected "
                "1 exterior or 2 internal)" % (nr, owners))
        corners_cyclic = np.asarray(
            [mesh[vertex].point for vertex in facet.vertices], dtype=float)
        if corners_cyclic.shape != (4, 3):
            raise ValueError(
                "vim.ChargeGram: internal HEX interface path requires "
                "quadrilateral facets")
        # NGSolve facets provide cyclic vertices.  C++ Q2 lattice order uses
        # [00,10,01,11], hence the 0,1,3,2 permutation.
        corners = corners_cyclic[[0, 1, 3, 2]]
        face_nodes.append(_QUAD_Q2_LINEAR_WEIGHTS @ corners)
        if int(p)==0:
            # FacetFESpace's constant test is one.  Assembly therefore gives
            # int_face M.n dS = int_[0,1]^2 q_ref du dv = q_ref directly,
            # including the signed two-sided jump on an internal interface.
            rows.append(jump_moments[dofs, :].tocsr())
            continue
        scale = max(np.linalg.norm(corners[1] - corners[0]),
                    np.linalg.norm(corners[2] - corners[0]), 1.0)
        affine_residual = np.linalg.norm(
            corners[3] - corners[1] - corners[2] + corners[0]) / scale
        if affine_residual > 1e-10:
            raise NotImplementedError(
                "vim.ChargeGram: internal HEX interfaces currently require "
                "affine straight quadrilateral facets; facet %d residual=%g"
                % (nr, affine_residual))
        weights = np.asarray([
            [(1.0-uu)*(1.0-vv), uu*(1.0-vv),
             (1.0-uu)*vv, uu*vv]
            for uu, vv in ref_points], dtype=float)
        physical_points = weights @ corners
        P = np.asarray([
            [x ** i * y ** j * z ** k for i, j, k in xyz_mons]
            for x, y, z in physical_points], dtype=float)
        local_to_physical = np.linalg.lstsq(Q, P, rcond=None)[0]
        D = np.column_stack(
            [values[dofs] / owners for values in moment_vectors])
        C = D @ np.linalg.pinv(local_to_physical)
        condition = np.linalg.cond(C)
        if not np.isfinite(condition) or condition > 1e12:
            raise RuntimeError(
                "vim.ChargeGram: ill-conditioned HEX facet moment conversion "
                "on facet %d (cond=%g)" % (nr, condition))
        # C contains physical-surface moments (dS = Js du dv), whereas the
        # HEX C++ kernel consumes the Piola-exact reference flux
        # q_ref = sigma * Js.  Js is constant on the affine facet.
        surface_jacobian = np.linalg.norm(
            np.cross(corners[1] - corners[0],
                     corners[2] - corners[0]))
        if not np.isfinite(surface_jacobian) or surface_jacobian <= 0.0:
            raise ValueError(
                "vim.ChargeGram: degenerate HEX facet %d" % nr)
        rows.append(sp.csr_matrix(
            surface_jacobian
            * np.linalg.solve(C, jump_moments[dofs, :].toarray())))

    return dict(
        B=sp.vstack(rows).tocsr(),
        face_nodes=face_nodes,
        mons=mons,
    )


def _charge_basis_hex(fes, cob_quad=3, *, materialize_mass=True,
                      internal_interfaces=False):
    """HEX analogue of `_charge_basis_curved`: charge map B + 27/9-node Q2 geometry nodes (via GetTrafo ->
    flat + curved ONE path).  fes = HDiv(hexmesh, order=0|1|2).  CALLER wraps TaskManager.  Flat NGSolve `.vol`
    hexes use direct NGSolve-reference lattice interpolation; curved meshes keep the GetTrafo source of truth.

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
    p = int(fes.globalorder)
    if p not in (0, 1, 2):
        raise ValueError("HEX HDiv-VIM supports HDiv order in {0,1,2}")
    mons_hex = _mons_hex(p)
    mons_quad = _mons_quad(p)
    t0 = time.perf_counter()
    nn = ng.specialcf.normal(mesh.dim)
    L2v = ng.L2(mesh, order=p)
    L2b = ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv = _csr(bv)                            # REF-measure moments of q_ref (the physical J cancels)
    Bb = _csr(bb)
    M_mass = _csr(mh) if materialize_mass else None
    t_assembly = time.perf_counter()

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = _exterior_bnd_elements(mesh)          # internal material-interface faces carry no charge
    rhp, rhw = _ref_prod_gauss(cob_quad, 3)
    rqp, rqw = _ref_prod_gauss(cob_quad, 2)
    ir_hex = ng.IntegrationRule(_Q2_LATTICE_3D, [1.0] * 27)
    ir_quad = ng.IntegrationRule(_Q2_LATTICE_2D, [1.0] * 9)
    # A linear topological mesh can still carry a live NGSolve deformation.
    # In that case the cheap vertex interpolation describes the undeformed
    # geometry, so the charge Gram lattice must come from GetTrafo.
    linear_lattice = (mesh.GetCurveOrder() < 2 and mesh.deformation is None)
    t_topology = time.perf_counter()

    host, kind, expo = [], [], []
    cell_nodes, face_nodes = [], []
    vol_transform_cache = {}
    face_transform_cache = {}
    vol_rows, vol_cols, vol_data = [], [], []
    face_rows, face_cols, face_data = [], [], []
    vol_lattice_s = 0.0
    vol_project_s = 0.0
    face_lattice_s = 0.0
    face_project_s = 0.0
    for c, e in enumerate(vels):
        _ts = time.perf_counter()
        cell_nodes.append(_hex_q2_lattice_nodes_ngsolve_linear(mesh, e)
                          if linear_lattice else _trafo_lattice_nodes(mesh, e, ir_hex))
        vol_lattice_s += time.perf_counter() - _ts
        _ts = time.perf_counter()
        fe = L2v.GetFE(e)
        key = _shape_signature_ref(fe, _HEX_SHAPE_SIG_PTS, dim=3)
        T = vol_transform_cache.get(key)
        if T is None:
            T = _ref_monomial_moment_transform(fe, mons_hex, rhp, rhw, dim=3)
            vol_transform_cache[key] = T
        dofs = list(L2v.GetDofNrs(e))
        base = c * len(mons_hex)
        for a, (i, j, k) in enumerate(mons_hex):
            host.append(c); kind.append(0); expo += [i, j, k]
            for b, col in enumerate(dofs):
                val = float(T[a, b])
                if val != 0.0:
                    vol_rows.append(base + a); vol_cols.append(int(col)); vol_data.append(val)
        vol_project_s += time.perf_counter() - _ts
    n_el = len(vels)
    _ts = time.perf_counter()
    Tv = sp.csr_matrix((vol_data, (vol_rows, vol_cols)), shape=(n_el * len(mons_hex), Bv.shape[0]))
    Bvol = Tv @ Bv
    vol_project_s += time.perf_counter() - _ts
    for f, e in enumerate(bels):
        _ts = time.perf_counter()
        face_nodes.append(_quad_q2_lattice_nodes_ngsolve_linear(mesh, e)
                          if linear_lattice else _trafo_lattice_nodes(mesh, e, ir_quad))
        face_lattice_s += time.perf_counter() - _ts
        _ts = time.perf_counter()
        fe = L2b.GetFE(e)
        key = _shape_signature_ref(fe, _QUAD_SHAPE_SIG_PTS, dim=2)
        T = face_transform_cache.get(key)
        if T is None:
            T = _ref_monomial_moment_transform(fe, mons_quad, rqp, rqw, dim=2)
            face_transform_cache[key] = T
        dofs = list(L2b.GetDofNrs(e))
        base = f * len(mons_quad)
        for a, (i, j) in enumerate(mons_quad):
            host.append(f); kind.append(1); expo += [i, j, 0]
            for b, col in enumerate(dofs):
                val = float(T[a, b])
                if val != 0.0:
                    face_rows.append(base + a); face_cols.append(int(col)); face_data.append(val)
        face_project_s += time.perf_counter() - _ts
    _ts = time.perf_counter()
    if bels:
        Tf = sp.csr_matrix((face_data, (face_rows, face_cols)), shape=(len(bels) * len(mons_quad), Bb.shape[0]))
        Bface = Tf @ Bb
    else:
        Bface = sp.csr_matrix((0, Bb.shape[1]))
    face_project_s += time.perf_counter() - _ts
    if internal_interfaces:
        internal = _broken_hex_face_charge_basis(fes, p)
        Bface = internal["B"]
        face_nodes = internal["face_nodes"]
        n_volume_charges = n_el * len(mons_hex)
        del host[n_volume_charges:]
        del kind[n_volume_charges:]
        del expo[3*n_volume_charges:]
        for f in range(len(face_nodes)):
            for i, j in mons_quad:
                host.append(f); kind.append(1); expo += [i, j, 0]
    t_before_vstack = time.perf_counter()
    B = sp.vstack([Bvol, Bface]).tocsr()
    t_after_vstack = time.perf_counter()
    return dict(B=B, M_mass=M_mass, M_mass_ngsolve=mh.mat,
                host=host, kind=kind, expo=expo, n_el=n_el,
                n_bf=len(face_nodes),
                cell_nodes=_f64_buffer(np.concatenate([n.ravel() for n in cell_nodes])),
                face_nodes=(_f64_buffer(np.concatenate([n.ravel() for n in face_nodes]))
                            if face_nodes else _EMPTY_F64),
                _timings={
                    "charge_basis_assembly_wall_s": t_assembly - t0,
                    "charge_basis_topology_wall_s": t_topology - t_assembly,
                    "charge_basis_vol_lattice_wall_s": vol_lattice_s,
                    "charge_basis_vol_project_wall_s": vol_project_s,
                    "charge_basis_face_lattice_wall_s": face_lattice_s,
                    "charge_basis_face_project_wall_s": face_project_s,
                    "charge_basis_vstack_wall_s": t_after_vstack - t_before_vstack,
                    "charge_basis_pack_wall_s": time.perf_counter() - t_after_vstack,
                    "charge_basis_lattice_mode": "ngsolve-linear-vol" if linear_lattice else "gettrafo",
                    "charge_basis_vol_transform_cache_size": len(vol_transform_cache),
                    "charge_basis_face_transform_cache_size": len(face_transform_cache),
                })


def _build_charge_gram_hex(fes, glout_n=None, glin_n=None, near_grade=0.5, far_inner=1.0,
                           eps=1e-12, leafsize=64, eta=2.0, image_masks=None, image_signs=None,
                           materialize_mass=True, build_hmatrix=True,
                           internal_interfaces=False):
    """Pure-hex BDM1/BDM2 charge Gram via the hex-mode C++ _ChargeGramHMatrix.  FLAT and CURVED (mesh.Curve(2))
    share ONE path (the 27-node Q2 lattice is extracted via GetTrafo either way -- the caller Curve(2)'s the
    mesh for curved).  glout_n = the 1D outer rule.  BDM1 keeps the validated default 4.  Flat BDM2 uses
    the TET-style analytic polynomial source moments plus a whole-host tensor outer for self/near hosts.
    Smooth far hosts use a reflection-invariant tensor-product rule on both complete reference domains,
    avoiding degree-six moment recurrences without changing the accepted spectrum or IMA contracts.
    Order 5 keeps the one-cell generalized Gram spectrum at 1 + 2e-6.
    Curved BDM2 retains order 6 and the reference-frame graded Duffy path.  The BDM1 default 4 was
    selected by the BDM1 hex spectrum/demag gates: 3 breaks the mass-normalized Gram spectrum, while 5/6
    were slower without improving the accepted affine/distorted regression cases.  glin_n = the 1D rule
    of the REF-frame RADIAL near/self inner (the PhiAtHO_Duffy port
    -- robust on distorted/curved hexes, where graded clouds left eig > 1 on the real cylinder mesh);
    eig(M_mass^-1 N) <= 1 gated on box AND cylinder meshes; block-memo build (~59x vs naive per-entry).
    far_inner = the PER-OUTER-POINT radial reach: an outer point farther than far_inner*size from a source
    sub-simplex integrates it with the cheap CACHED far cloud instead of the per-point radial rule.  1.0
    keeps the radial on self / genuinely touching near geometry while moving smooth shell points onto the
    cached far cloud; the affine/distorted BDM1 spectrum gates and the cube N=8/10 demag regression keep this
    from silently becoming a coarse approximation.  The far TET cloud is the SAME Keast-15 degree-5 rule as
    the outer (a degree-3 WV rule at reach 1.5 ate the flat-cylinder eig margin, 1.0005 -> 1.0044; Keast-15
    restores 1.0006 for +4%% build).  The build also skips the strictly-lower H-matrix leaves (symmetric fill
    -- every apply of the Gram routes through the exactly-symmetric matvec, so they are never read)."""
    t0 = time.perf_counter()
    p = int(fes.globalorder)
    default_outer = 4 if p == 1 else (5 if fes.mesh.GetCurveOrder() < 2 else 6)
    glout_n = default_outer if glout_n is None else int(glout_n)
    glin_n = (5 if p == 1 else 7) if glin_n is None else int(glin_n)
    cb = _charge_basis_hex(
        fes, cob_quad=max(3, p+1), materialize_mass=materialize_mass,
        internal_interfaces=bool(internal_interfaces))
    t1 = time.perf_counter()
    glo, gwo = _g01(glout_n)
    gli, gwi = _g01(glin_n)
    ftp = np.asarray(_SYM5_TET[0]); ftw = np.asarray(_SYM5_TET[1])
    G = _rp._ChargeGramHMatrix(
        hex_cell_nodes=cb["cell_nodes"], quad_face_nodes=cb["face_nodes"],
        n_el=int(cb["n_el"]), n_bf=int(cb["n_bf"]),
        charge_host=_i32_buffer(cb["host"]), charge_kind=_i32_buffer(cb["kind"]),
        charge_expo=_i32_buffer(cb["expo"]),
        sym_tet_pts=_f64_buffer(_SYM5_TET[0]), sym_tet_w=_f64_buffer(_SYM5_TET[1]),
        sym_tri_pts=_f64_buffer(_SYM5_TRI[0]), sym_tri_w=_f64_buffer(_SYM5_TRI[1]),
        gl_out=_f64_buffer(glo), gw_out=_f64_buffer(gwo),
        gl_in=_f64_buffer(gli), gw_in=_f64_buffer(gwi),
        far_tet_pts=_f64_buffer(ftp), far_tet_w=_f64_buffer(ftw),
        far_tri_pts=_f64_buffer(_SYM5_TRI[0]), far_tri_w=_f64_buffer(_SYM5_TRI[1]),
        near_grade=near_grade, far_inner_factor=far_inner,
        image_masks=(_EMPTY_I32 if image_masks is None else _i32_buffer(image_masks)),
        image_signs=(_EMPTY_F64 if image_signs is None else _f64_buffer(image_signs)),
        eps=eps, leaf=leafsize, eta=eta, build=bool(build_hmatrix))
    t2 = time.perf_counter()
    chk = G.hex_state_check()
    t3 = time.perf_counter()
    build_charge_gram.last_timings = {
        "charge_basis_wall_s": t1 - t0,
        "charge_gram_cpp_wall_s": t2 - t1,
        "hex_state_check_wall_s": t3 - t2,
    }
    build_charge_gram.last_timings.update(cb.get("_timings", {}))
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "hex charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption "
            "(0xc0000374 class) -- do NOT trust this Gram; abort and report the incident.")
    return cb["B"], G, cb["M_mass"], cb["M_mass_ngsolve"]


_TRI6_LAT = [(1, 0), (0, 1), (0, 0), (0.5, 0.5), (0, 0.5), (0.5, 0)]


def _tri_mons_2d(degree):
    return [(total-j, j) for total in range(degree+1) for j in range(total+1)]


def _quad_mons_2d(degree):
    return [(i, j) for j in range(degree+1) for i in range(degree+1)]


def _edge_mons_2d(degree):
    return [(i,) for i in range(degree+1)]


def _fit_geometry_map_2d(mesh, element, cell_type, degree):
    """Fit the NGSolve isoparametric map in the monomial basis consumed by C++."""
    if cell_type == 0:
        mons = _tri_mons_2d(degree)
        points = [(i/degree, j/degree)
                  for j in range(degree+1) for i in range(degree+1-j)]
    elif cell_type == 1:
        mons = _quad_mons_2d(degree)
        points = [(i/degree, j/degree)
                  for j in range(degree+1) for i in range(degree+1)]
    else:
        mons = _edge_mons_2d(degree)
        points = [(i/degree,) for i in range(degree+1)]
    ir = ng.IntegrationRule(
        [(p[0], p[1] if len(p) > 1 else 0.0, 0.0) for p in points],
        [1.0]*len(points))
    physical = _trafo_lattice_nodes(mesh, element, ir)[:, :2]
    if cell_type < 2:
        vandermonde = np.array([
            [(p[0]**i)*(p[1]**j) for i, j in mons] for p in points])
    else:
        vandermonde = np.array([[p[0]**i for (i,) in mons] for p in points])
    return np.linalg.solve(vandermonde, physical)


def _charge_basis_2d(fes, cob_quad=3, *, materialize_mass=True):
    """2D analogue of `_charge_basis_hex` (motor cross-sections; memory hdiv-vim-tri-quad-motor): charge
    map B + Q1..Q3 polynomial geometry for tri/quad cells and boundary edges, all in the NGSolve REF frame with
    the Piola-exact extraction (the dimension-independent J-cancellation identity).  Kernel side is the
    2D log Gram (C++ dim2 mode).  CALLER wraps TaskManager."""
    mesh = fes.mesh
    p = int(fes.globalorder)
    if p not in (1, 2):
        raise ValueError("2D HDiv-VIM supports HDiv order in {1,2}")
    geometry_order = max(1, int(mesh.GetCurveOrder()))
    tri_mons = _tri_mons_2d(p-1)
    quad_mons = _quad_mons_2d(p)
    edge_mons = _edge_mons_2d(p)
    bonus_intorder = 2*geometry_order
    nn2 = ng.specialcf.normal(2)
    L2v = ng.L2(mesh, order=p)
    Sb2 = ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v)
    bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx(bonus_intorder=bonus_intorder)
    bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=Sb2)
    bb += (u.Trace() * nn2) * Sb2.TestFunction() * ng.ds(bonus_intorder=bonus_intorder)
    bb.Assemble()
    mh = ng.BilinearForm(fes)
    mh += u * fes.TestFunction() * ng.dx(bonus_intorder=bonus_intorder)
    mh.Assemble()
    Bv = _csr(bv); Bb = _csr(bb)
    M_mass = _csr(mh) if materialize_mass else None
    if M_mass is not None and M_mass.nnz:
        mass_zero_floor = 4096*np.finfo(float).eps*np.max(np.abs(M_mass.data))
        M_mass.data[np.abs(M_mass.data) <= mass_zero_floor] = 0.0
        M_mass.eliminate_zeros()

    g, gw = _g01(cob_quad)
    tp = np.array([[uu, vv*(1 - uu)] for uu in g for vv in g])            # Duffy on the NGSolve tri ref
    tw = np.array([wu*wv*(1 - uu) for uu, wu in zip(g, gw) for vv, wv in zip(g, gw)])
    qp = np.array([[uu, vv] for uu in g for vv in g])
    qw = np.array([wu*wv for _, wu in zip(g, gw) for _, wv in zip(g, gw)])
    ep = g.reshape(-1, 1); ew = gw

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = _exterior_bnd_elements(mesh)          # internal material-interface faces carry no charge
    Brows, host, kind, expo = [], [], [], []
    cell_map, cell_type, edge_map = [], [], []
    for c, e in enumerate(vels):
        quad = len(mesh[e].vertices) == 4
        ct = 1 if quad else 0
        coeff = _fit_geometry_map_2d(mesh, e, ct, geometry_order)
        slot = np.zeros(((geometry_order+1)**2, 2))
        slot[:coeff.shape[0]] = coeff
        cell_map.append(slot); cell_type.append(ct)
        fe = L2v.GetFE(e)
        rp, rw = (qp, qw) if quad else (tp, tw)
        Phi = np.array([fe.CalcShape(pt[0], pt[1], 0.0) for pt in rp])
        Mref = (Phi * rw[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bv[list(L2v.GetDofNrs(e)), :].toarray())
        mons = quad_mons if quad else tri_mons
        Sv = _change_of_basis_ref(fe, mons, rp, rw, dim=2)
        blk = Sv @ rows
        for a, (mi, mj) in enumerate(mons):
            Brows.append(sp.csr_matrix(blk[a])); host.append(c); kind.append(0); expo += [mi, mj, 0]
    for f, e in enumerate(bels):
        edge_map.append(_fit_geometry_map_2d(mesh, e, 2, geometry_order))
        fe = Sb2.GetFE(e)
        Phi = np.array([fe.CalcShape(t[0], 0.0, 0.0) for t in ep])
        Mref = (Phi * ew[:, None]).T @ Phi
        rows = np.linalg.solve(Mref, Bb[list(Sb2.GetDofNrs(e)), :].toarray())
        # dim=1 change of basis (S = M_mono^-1 C), mirroring _change_of_basis_ref
        Mm = np.array([[np.sum(ew * ep[:, 0]**(mi + mj)) for (mj,) in edge_mons] for (mi,) in edge_mons])
        Cm = np.array([[np.sum(ew * (ep[:, 0]**mi) * Phi[:, k]) for k in range(Phi.shape[1])]
                       for (mi,) in edge_mons])
        Se = np.linalg.solve(Mm, Cm)
        blk = Se @ rows
        for a, (mi,) in enumerate(edge_mons):
            Brows.append(sp.csr_matrix(blk[a])); host.append(f); kind.append(1); expo += [mi, 0, 0]
    B = sp.vstack(Brows).tocsr()
    return dict(B=B, M_mass=M_mass, M_mass_ngsolve=mh.mat,
                host=host, kind=kind, expo=expo,
                n_el=len(vels), n_be=len(bels),
                geometry_order=geometry_order,
                cell_map=_f64_buffer(np.concatenate([n.ravel() for n in cell_map])),
                cell_type=_i32_buffer(cell_type),
                edge_map=_f64_buffer(np.concatenate([n.ravel() for n in edge_map])))


def _prod_tri01(n):
    """product-Gauss bary rule on the unit simplex (lam1 = u, lam2 = v(1-u); W sums 1/2)."""
    g, gw = _g01(n)
    P, W = [], []
    for u, wu in zip(g, gw):
        for v, wv in zip(g, gw):
            P.append((u, v*(1 - u))); W.append(wu*wv*(1 - u))
    return np.array(P), np.array(W)


def _build_charge_gram_2d(fes, outer_n=None, glin_n=None, gledge_n=None, near_grade=0.6, far_inner=1.5,
                          eps=1e-12, leafsize=64, eta=2.0,
                          image_masks=None, image_signs=None, materialize_mass=True,
                          build_hmatrix=True):
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
    p = int(fes.globalorder)
    outer_n = (4 if p == 1 else 6) if outer_n is None else int(outer_n)
    glin_n = (8 if p == 1 else 10) if glin_n is None else int(glin_n)
    gledge_n = (12 if p == 1 else 16) if gledge_n is None else int(gledge_n)
    cb = _charge_basis_2d(fes, cob_quad=max(3, p+1), materialize_mass=materialize_mass)
    otp, otw = _prod_tri01(outer_n)
    glq, gwq = _g01(outer_n)
    gli, gwi = _g01(glin_n)
    gle, gwe = _g01(gledge_n)
    G = _rp._ChargeGramHMatrix(
        dim2=2,
        geometry_order=int(cb["geometry_order"]), cell_map=cb["cell_map"],
        cell_type=cb["cell_type"], edge_map=cb["edge_map"],
        n_el=int(cb["n_el"]), n_be=int(cb["n_be"]),
        charge_host=_i32_buffer(cb["host"]), charge_kind=_i32_buffer(cb["kind"]),
        charge_expo=_i32_buffer(cb["expo"]),
        sym_tri_pts=_f64_buffer(otp), sym_tri_w=_f64_buffer(otw),
        gl_quad=_f64_buffer(glq), gw_quad=_f64_buffer(gwq),
        gl_edge=_f64_buffer(gle), gw_edge=_f64_buffer(gwe),
        gl_in=_f64_buffer(gli), gw_in=_f64_buffer(gwi),
        far_tri_pts=_f64_buffer(_SYM5_TRI[0]), far_tri_w=_f64_buffer(_SYM5_TRI[1]),
        near_grade=near_grade, far_inner_factor=far_inner,
        image_masks=(_EMPTY_I32 if image_masks is None else _i32_buffer(image_masks)),
        image_signs=(_EMPTY_F64 if image_signs is None else _f64_buffer(image_signs)),
        # Small BDM2 reduced models must remain dense: an ACA leaf can amplify a
        # nominal 1e-12 block error through high-mu solves and destroy the IMA
        # roundoff contract.  Larger models still split/compress normally.
        eps=eps, leaf=max(64, int(leafsize)), eta=eta, build=bool(build_hmatrix))
    chk = G.hex_state_check()
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "2D charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption "
            "(0xc0000374 class) -- do NOT trust this Gram; abort and report the incident.")
    return cb["B"], G, cb["M_mass"], cb["M_mass_ngsolve"]


# WEDGE (PRISM) BDM1/BDM2: the prism div-image is tri-Pp (x) z-Pp: 6 modes at BDM1, 18 at BDM2;
# boundary faces are MIXED tri (SurfaceL2 Pp, 3/6) + quad (SurfaceL2 Qp, 4/9).  Geometry = the 18-node tri-P2
# (x) z-P2 lattice (n = t + 6*iz, t = _TRI6_LAT node) via GetTrafo -> flat + curved ONE path.
def _mons_tri(order):
    return [(total-j, j) for total in range(order+1) for j in range(total+1)]


def _mons_wedge(order):
    return [(i, j, k) for k in range(order+1) for i, j in _mons_tri(order)]
_WEDGE_Q2_LATTICE = [(_TRI6_LAT[t][0], _TRI6_LAT[t][1], iz / 2.0) for iz in range(3) for t in range(6)]
_TRI_SHAPE_SIG_PTS = np.array([(0.137, 0.239), (0.521, 0.211), (0.319, 0.173)], dtype=float)
_WEDGE_SHAPE_SIG_PTS = np.array([(0.137, 0.239, 0.361), (0.521, 0.211, 0.557),
                                 (0.319, 0.173, 0.823)], dtype=float)
_WEDGE_Q2_LINEAR_WEIGHTS = np.array([
    [(1.0-w)*u, (1.0-w)*v, (1.0-w)*(1.0-u-v), w*u, w*v, w*(1.0-u-v)]
    for u, v, w in _WEDGE_Q2_LATTICE
], dtype=float)
_TRI_Q2_LINEAR_WEIGHTS = np.array([[u, v, 1.0-u-v] for u, v in _TRI6_LAT], dtype=float)


def _wedge_q2_lattice_nodes_ngsolve_linear(mesh, e):
    """18 lattice nodes for a linear NGSolve prism, in the C++ wedge Q2 lattice order."""
    return _WEDGE_Q2_LINEAR_WEIGHTS @ _mesh_vertices_array(mesh, e)


def _tri_q2_lattice_nodes_ngsolve_linear(mesh, e):
    """6 lattice nodes for a linear NGSolve boundary triangle, in the C++ tri6 lattice order."""
    return _TRI_Q2_LINEAR_WEIGHTS @ _mesh_vertices_array(mesh, e)


def _prism_cob_quad(nz=3):
    """Ref-prism change-of-basis quadrature: SYM5 tri (u,v) x nz-pt Gauss (w).

    It integrates the degree<=4 BDM2 L2-shape x monomial products exactly in
    the triangular and axial reference coordinates.
    """
    tp, tw = np.asarray(_SYM5_TRI[0]), np.asarray(_SYM5_TRI[1])
    gz, gwz = _g01(nz)
    P, W = [], []
    for (tu, tv), twt in zip(tp, tw):
        for z, wz in zip(gz, gwz):
            P.append((tu, tv, z)); W.append(twt * wz)
    return np.array(P), np.array(W)


def _charge_basis_wedge(fes, *, materialize_mass=True):
    """WEDGE (prism) analogue of `_charge_basis_hex`: charge map B + 18-node prism cell nodes + MIXED
    tri(6-node)/quad(9-node) face nodes (packed in 27-double 9-node slots, a tri fills the first 6) + a
    per-face type array, all via GetTrafo (flat + curved ONE path).  fes = HDiv(prismmesh, order=1|2).
    CALLER wraps TaskManager.  Piola-exact charge model, exactly as the hex path (the J's cancel)."""
    t0 = time.perf_counter()
    mesh = fes.mesh
    p = int(fes.globalorder)
    if p not in (1, 2):
        raise ValueError("WEDGE HDiv-VIM supports HDiv order in {1,2}")
    mons_wedge = _mons_wedge(p)
    mons_tri = _mons_tri(p)
    mons_quad = _mons_quad(p)
    nn = ng.specialcf.normal(mesh.dim)
    L2v = ng.L2(mesh, order=p)
    L2b = ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds; bb.Assemble()
    mh = ng.BilinearForm(fes); mh += u * fes.TestFunction() * ng.dx; mh.Assemble()
    Bv = _csr(bv); Bb = _csr(bb)
    M_mass = _csr(mh) if materialize_mass else None
    t_forms = time.perf_counter()

    vels = [ng.ElementId(ng.VOL, i) for i in range(mesh.GetNE(ng.VOL))]
    bels = _exterior_bnd_elements(mesh)          # internal material-interface faces carry no charge
    php, phw = _prism_cob_quad()                                  # prism ref change-of-basis quadrature
    qp2, qw2 = _ref_prod_gauss(3, 2)                             # quad-face ref quadrature ([0,1]^2)
    tp2, tw2 = np.asarray(_SYM5_TRI[0]), np.asarray(_SYM5_TRI[1])  # tri-face ref quadrature (ref tri)
    ir_wedge = ng.IntegrationRule(_WEDGE_Q2_LATTICE, [1.0] * 18)
    ir_quad = ng.IntegrationRule(_Q2_LATTICE_2D, [1.0] * 9)
    ir_tri = ng.IntegrationRule(_TRI6_LAT, [1.0] * 6)

    Brows, host, kind, expo = [], [], [], []
    cell_nodes, face_nodes, face_type = [], [], []
    t_setup = time.perf_counter()
    # See the HEX path above: SetDeformation changes GetTrafo without changing
    # GetCurveOrder, hence a deformed linear WEDGE needs transformed nodes.
    linear_lattice = (mesh.GetCurveOrder() < 2 and mesh.deformation is None)
    vol_transform_cache = {}
    face_transform_cache = {}
    vol_rows, vol_cols, vol_data = [], [], []
    face_rows, face_cols, face_data = [], [], []
    vol_lattice_s = 0.0
    vol_project_s = 0.0
    face_lattice_s = 0.0
    face_project_s = 0.0
    for c, e in enumerate(vels):
        _ts = time.perf_counter()
        cell_nodes.append(_wedge_q2_lattice_nodes_ngsolve_linear(mesh, e)
                          if linear_lattice else _trafo_lattice_nodes(mesh, e, ir_wedge))
        vol_lattice_s += time.perf_counter() - _ts
        _ts = time.perf_counter()
        fe = L2v.GetFE(e)
        key = _shape_signature_ref(fe, _WEDGE_SHAPE_SIG_PTS, dim=3)
        T = vol_transform_cache.get(key)
        if T is None:
            T = _ref_monomial_moment_transform(fe, mons_wedge, php, phw, dim=3)
            vol_transform_cache[key] = T
        dofs = list(L2v.GetDofNrs(e))
        base = c * len(mons_wedge)
        for a, (i, j, k) in enumerate(mons_wedge):
            host.append(c); kind.append(0); expo += [i, j, k]
            for b, col in enumerate(dofs):
                val = float(T[a, b])
                if val != 0.0:
                    vol_rows.append(base + a); vol_cols.append(int(col)); vol_data.append(val)
        vol_project_s += time.perf_counter() - _ts
    n_el = len(vels)
    _ts = time.perf_counter()
    Tv = sp.csr_matrix((vol_data, (vol_rows, vol_cols)), shape=(n_el * len(mons_wedge), Bv.shape[0]))
    Bvol = Tv @ Bv
    vol_project_s += time.perf_counter() - _ts
    t_vol = time.perf_counter()
    face_charge_count = 0
    for f, e in enumerate(bels):
        _ts = time.perf_counter()
        if len(list(mesh[e].vertices)) == 3:                          # TRI face -> 6-node, P1, 1 sub-tri
            nd = (_tri_q2_lattice_nodes_ngsolve_linear(mesh, e)
                  if linear_lattice else _trafo_lattice_nodes(mesh, e, ir_tri))
            slot = np.zeros((9, 3)); slot[:6] = nd
            face_nodes.append(slot); face_type.append(0)
            face_lattice_s += time.perf_counter() - _ts
            _ts = time.perf_counter()
            fe = L2b.GetFE(e)
            key = ("tri", _shape_signature_ref(fe, _TRI_SHAPE_SIG_PTS, dim=2))
            T = face_transform_cache.get(key)
            if T is None:
                T = _ref_monomial_moment_transform(fe, mons_tri, tp2, tw2, dim=2)
                face_transform_cache[key] = T
            dofs = list(L2b.GetDofNrs(e))
            base = face_charge_count
            for a, (i, j) in enumerate(mons_tri):
                host.append(f); kind.append(1); expo += [i, j, 0]
                for b, col in enumerate(dofs):
                    val = float(T[a, b])
                    if val != 0.0:
                        face_rows.append(base + a); face_cols.append(int(col)); face_data.append(val)
            face_charge_count += len(mons_tri)
            face_project_s += time.perf_counter() - _ts
        else:                                                         # QUAD face -> 9-node, Q1, 2 sub-tris
            face_nodes.append(_quad_q2_lattice_nodes_ngsolve_linear(mesh, e)
                              if linear_lattice else _trafo_lattice_nodes(mesh, e, ir_quad))
            face_type.append(1)
            face_lattice_s += time.perf_counter() - _ts
            _ts = time.perf_counter()
            fe = L2b.GetFE(e)
            key = ("quad", _shape_signature_ref(fe, _QUAD_SHAPE_SIG_PTS, dim=2))
            T = face_transform_cache.get(key)
            if T is None:
                T = _ref_monomial_moment_transform(fe, mons_quad, qp2, qw2, dim=2)
                face_transform_cache[key] = T
            dofs = list(L2b.GetDofNrs(e))
            base = face_charge_count
            for a, (i, j) in enumerate(mons_quad):
                host.append(f); kind.append(1); expo += [i, j, 0]
                for b, col in enumerate(dofs):
                    val = float(T[a, b])
                    if val != 0.0:
                        face_rows.append(base + a); face_cols.append(int(col)); face_data.append(val)
            face_charge_count += len(mons_quad)
            face_project_s += time.perf_counter() - _ts
    _ts = time.perf_counter()
    if face_charge_count:
        Tf = sp.csr_matrix((face_data, (face_rows, face_cols)), shape=(face_charge_count, Bb.shape[0]))
        Bface = Tf @ Bb
    else:
        Bface = sp.csr_matrix((0, Bb.shape[1]))
    face_project_s += time.perf_counter() - _ts
    t_face = time.perf_counter()
    B = sp.vstack([Bvol, Bface]).tocsr()
    t_vstack = time.perf_counter()
    return dict(B=B, M_mass=M_mass, M_mass_ngsolve=mh.mat,
                host=host, kind=kind, expo=expo, n_el=n_el, n_bf=len(bels),
                cell_nodes=_f64_buffer(np.concatenate([n.ravel() for n in cell_nodes])),
                face_nodes=_f64_buffer(np.concatenate([n.ravel() for n in face_nodes])),
                face_type=_i32_buffer(face_type),
                _timings={
                    "charge_basis_forms_wall_s": t_forms - t0,
                    "charge_basis_setup_wall_s": t_setup - t_forms,
                    "charge_basis_vol_wall_s": t_vol - t_setup,
                    "charge_basis_face_wall_s": t_face - t_vol,
                    "charge_basis_vol_lattice_wall_s": vol_lattice_s,
                    "charge_basis_vol_project_wall_s": vol_project_s,
                    "charge_basis_face_lattice_wall_s": face_lattice_s,
                    "charge_basis_face_project_wall_s": face_project_s,
                    "charge_basis_vstack_wall_s": t_vstack - t_face,
                    "charge_basis_pack_wall_s": time.perf_counter() - t_vstack,
                    "charge_basis_lattice_mode": "ngsolve-linear-prism" if linear_lattice else "gettrafo",
                    "charge_basis_vol_transform_cache_size": len(vol_transform_cache),
                    "charge_basis_face_transform_cache_size": len(face_transform_cache),
                    "charge_basis_wall_s": time.perf_counter() - t0,
                })


def _build_charge_gram_wedge(fes, glout_n=None, glin_n=None, near_grade=0.6, far_inner=1.5,
                             eps=1e-12, leafsize=64, eta=2.0, image_masks=None, image_signs=None,
                             materialize_mass=True, build_hmatrix=True):
    """Pure-prism BDM1/BDM2 charge Gram via the wedge-mode C++ _ChargeGramHMatrix (mirror of _build_charge_gram_hex;
    FLAT + Curve(2) share ONE path).  numpy de-risk eig(M_mass^-1 N) in [0,1]: 0.989 @ n=2, 0.997 @ n=3;
    demag_z ~ 1/3.  The wedge mode shares the hex block memo / symmetric-fill build, so the golden hex path
    is byte-for-byte untouched."""
    t0 = time.perf_counter()
    p = int(fes.globalorder)
    glout_n = (6 if p == 1 else 7) if glout_n is None else int(glout_n)
    glin_n = (5 if p == 1 else 8) if glin_n is None else int(glin_n)
    cb = _charge_basis_wedge(fes, materialize_mass=materialize_mass)
    t1 = time.perf_counter()
    glo, gwo = _g01(glout_n); gli, gwi = _g01(glin_n)
    field_tri_rule = ng.IntegrationRule(ng.ET.TRIG, 5)
    field_tri_pts = np.asarray(
        [(ip.point[0], ip.point[1]) for ip in field_tri_rule], dtype=float)
    field_tri_w = np.asarray([ip.weight for ip in field_tri_rule], dtype=float)
    G = _rp._ChargeGramHMatrix(
        wedge_cell_nodes=cb["cell_nodes"], face_nodes=cb["face_nodes"], face_type=cb["face_type"],
        n_el=int(cb["n_el"]), n_bf=int(cb["n_bf"]),
        charge_host=_i32_buffer(cb["host"]), charge_kind=_i32_buffer(cb["kind"]),
        charge_expo=_i32_buffer(cb["expo"]),
        sym_tet_pts=_f64_buffer(_SYM5_TET[0]), sym_tet_w=_f64_buffer(_SYM5_TET[1]),
        sym_tri_pts=_f64_buffer(_SYM5_TRI[0]), sym_tri_w=_f64_buffer(_SYM5_TRI[1]),
        field_tri_pts=_f64_buffer(field_tri_pts), field_tri_w=_f64_buffer(field_tri_w),
        gl_out=_f64_buffer(glo), gw_out=_f64_buffer(gwo),
        gl_in=_f64_buffer(gli), gw_in=_f64_buffer(gwi),
        far_tet_pts=_f64_buffer(_SYM5_TET[0]), far_tet_w=_f64_buffer(_SYM5_TET[1]),
        far_tri_pts=_f64_buffer(_SYM5_TRI[0]), far_tri_w=_f64_buffer(_SYM5_TRI[1]),
        near_grade=near_grade, far_inner_factor=far_inner,
        image_masks=(_EMPTY_I32 if image_masks is None else _i32_buffer(image_masks)),
        image_signs=(_EMPTY_F64 if image_signs is None else _f64_buffer(image_signs)),
        eps=eps, leaf=leafsize, eta=eta, build=bool(build_hmatrix))
    t2 = time.perf_counter()
    chk = G.hex_state_check()
    t3 = time.perf_counter()
    build_charge_gram.last_timings = {
        "charge_basis_wall_s": t1 - t0,
        "charge_gram_cpp_wall_s": t2 - t1,
        "hex_state_check_wall_s": t3 - t2,
    }
    build_charge_gram.last_timings.update(cb.get("_timings", {}))
    if chk["ctor"] != chk["now"]:
        raise RuntimeError(
            "wedge charge Gram instance state was corrupted between construction and use "
            f"(canary ctor={chk['ctor']!r} != now={chk['now']!r}): heap corruption (0xc0000374 class) -- "
            "do NOT trust this Gram; abort and report the incident.")
    return cb["B"], G, cb["M_mass"], cb["M_mass_ngsolve"]


def _configure_cpp_operator(B, G, M_mass, M_mass_ngsolve):
    """Pin geometry/FESpace sparse topology in the persistent C++ operator."""
    B = sp.csr_matrix(B)
    if M_mass is not None:
        M_mass = sp.csr_matrix(M_mass)
    G.configure_charge_map(
        _i32_buffer(B.indptr), _i32_buffer(B.indices), _f64_buffer(B.data), int(B.shape[1]))
    G.configure_mass_matrix_ngsolve(M_mass_ngsolve)
    G.configure_geometry_mass_matrix_ngsolve(M_mass_ngsolve)
    return B, G, M_mass


def build_charge_gram(fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0, far_quad=3, ho_far_factor=2.0,
                      inner_quad=None, curve_order=None, curve_gauss=8, nonlinear=False,
                      image_masks=None, image_signs=None, _materialize_mass=True,
                      _build_hmatrix=True, internal_interfaces=False):
    """From an HDiv FESpace (order p, the order from the fes), build the monomial charge-density map
    B (scipy CSR, n_charge x ndof), the C++ charge-Gram H-matrix G, and the HDiv mass M_mass (CSR).
    The CALLER wraps in TaskManager.

    For TET, curve_order=None matches the mesh, 0 forces the flat diagnostic path, and 2 selects the
    isoparametric-P2 curved charge Gram.  HEX/WEDGE always read the active geometry through GetTrafo, so their
    flat/curved route is automatic.  The TET curved path uses a reference-frame charge map B and the C++ Duffy
    Gram (curve_gauss = the inner Gauss-Legendre pts/dim, 8 for the production rule).  Matching by default
    prevents callers such as the history solver from combining a curved NGSolve mass matrix with a straight
    TET charge Gram.  curve_order helps
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
    build_charge_gram.last_timings = {}
    mesh = fes.mesh
    p = int(fes.globalorder)
    _vtypes = _volume_vertex_counts(mesh)
    validate_hdiv_configuration(mesh.dim, _vtypes, p, mesh.GetCurveOrder())
    if internal_interfaces and not (
            mesh.dim == 3 and _vtypes in ({4}, {8})
            and mesh.GetCurveOrder() < 2):
        raise NotImplementedError(
            "vim.ChargeGram: internal_interfaces=True is currently wired "
            "for straight pure-TET and affine-facet pure-HEX meshes; "
            "WEDGE/curved topology interfaces require their matching facet "
            "geometry kernels")
    if mesh.dim == 3 and _vtypes == {4}:
        if curve_order is None:
            mesh_curve_order = int(mesh.GetCurveOrder())
            curve_order = mesh_curve_order if mesh_curve_order >= 2 else None
        elif int(curve_order) == 0:
            curve_order = None
    if curve_order is not None and int(curve_order) > 0 and mesh.GetCurveOrder() < int(curve_order):
        raise ValueError(
            "vim.ChargeGram: curve_order=%d requires mesh geometry order >= %d; got %d."
            % (int(curve_order), int(curve_order), int(mesh.GetCurveOrder())))
    image_masks = [] if image_masks is None else list(image_masks)   # robust for NumPy arrays (truth-value)
    image_signs = [] if image_signs is None else list(image_signs)
    if len(image_masks) != len(image_signs):
        raise ValueError("vim.ChargeGram: image_masks and image_signs must have the same length")
    if image_masks:
        # IMA (mirror-image charge folding): wired for flat/Curve(2) pure-TET (C++ m_highorder QuadDotRefl->PhiInner)
        # AND pure-HEX / pure-WEDGE (the QuadBlockHex/Wedge(mask) reflected block), plus the planar log kernel.
        # Fail loud on unsupported topology rather than silently dropping an image term.
        _ivt = _vtypes if mesh.dim == 3 else None
        if mesh.dim == 2:
            if any(int(mask) < 1 or int(mask) > 3 for mask in image_masks):
                raise ValueError("vim.ChargeGram: planar IMA masks use x/y bits only (1..3).")
        elif _ivt not in ({4}, {8}, {6}):
            raise ValueError(
                "vim.ChargeGram: image_masks (IMA) is wired for flat/Curve(2) pure-TET / pure-HEX / pure-WEDGE "
                "BDM1/BDM2 Gram; 2D-planar reduced models use the planar image path.  "
                "(got dim=%s, vtypes=%s, curve_order=%r)."
                % (mesh.dim, sorted(_ivt) if _ivt else None, curve_order))
    if mesh.dim == 2:
        # 2D PLANAR (motor cross-section) layer: BDM1/BDM2 tri/quad cells + boundary-edge charges, log kernel.
        return _configure_cpp_operator(
            *_build_charge_gram_2d(
                fes, eps=eps, leafsize=leafsize, eta=eta,
                image_masks=image_masks, image_signs=image_signs,
                materialize_mass=_materialize_mass, build_hmatrix=_build_hmatrix))
    if _vtypes == {8}:
        # PURE-HEX BDM1/BDM2: tensor Qp charge basis + Q2 geometry; FLAT or Curve(2) one path.
        # curve_order is IGNORED for hex -- curved is automatic (GetTrafo picks up mesh.Curve(2)); the caller
        # Curve(2)'s the mesh before this call, exactly like the tet curved path.  Uses the hex-gated params.
        return _configure_cpp_operator(*_build_charge_gram_hex(
            fes, eps=eps, leafsize=leafsize, eta=eta,
            image_masks=image_masks, image_signs=image_signs,
            materialize_mass=_materialize_mass, build_hmatrix=_build_hmatrix,
            internal_interfaces=bool(internal_interfaces)))
    if _vtypes == {6}:
        # PURE-WEDGE (PRISM) BDM1/BDM2: tri-Pp x z-Pp volume charge + mixed tri/quad-face
        # surface charge; 18-node Q2 geometry; FLAT or Curve(2) one path).  curve_order is IGNORED (curved is
        # automatic via GetTrafo picking up mesh.Curve(2)), same as the hex path.
        return _configure_cpp_operator(*_build_charge_gram_wedge(
            fes, eps=eps, leafsize=leafsize, eta=eta,
            image_masks=image_masks, image_signs=image_signs,
            materialize_mass=_materialize_mass, build_hmatrix=_build_hmatrix))
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
    # DEFAULT depends on the USE (2026-06-30).  The LINEAR demag only needs the 3*p PSD FLOOR (BDM1 -> quad=3):
    # that makes the NEAR build ~1.8-2.1x cheaper (the near U-list is 98% of the build = the dominant lever) and
    # is VALIDATED to preserve demag (7e-6), per-element leak + magnetic moment, and PSD (min eig ~0).  The
    # NONLINEAR energy-Newton KEEPS the +1 margin (max(3*p,4) -> quad=4 for BDM1): the energy HESSIAN wants more
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
            else (max(3 * p, 4) if nonlinear else max(2, 3 * p)))
    if curve_order is not None:
        # CURVED (isoparametric P2) Gram: curved charge map B (reference-frame change-of-basis) + the C++
        # curved-Duffy charge Gram.  Only P2 is wired (the C++ CurvedTet/TriPotential are P2); the mesh must
        # already be mesh.Curve(2)'d (the caller does it, per the NGSolve convention).
        if int(curve_order) != 2:
            raise NotImplementedError("vim.ChargeGram: only curve_order=2 (isoparametric P2) is wired "
                                      "(the C++ CurvedTet/TriPotential are P2); got %r." % (curve_order,))
        cbk = _charge_basis_curved(fes, quad, materialize_mass=_materialize_mass)
        (rtp, rtw), (rsp, rsw) = _curved_outer_rules(p)
        gx, gw = _g01(int(curve_gauss))
        kw = dict(
            cell_nodes=cbk["cell_nodes"], face_nodes=cbk["face_nodes"],
            cell_vertices=cbk["cell_vertices"], face_vertices=cbk["face_vertices"],
            n_el=cbk["n_el"], curve_order=2,
            charge_host=_i32_buffer(cbk["host"]), charge_kind=_i32_buffer(cbk["kind"]),
            charge_expo=_i32_buffer(cbk["expo"]),
            ref_tet_pts=_f64_buffer(rtp), ref_tet_w=_f64_buffer(rtw),
            ref_tri_pts=_f64_buffer(rsp), ref_tri_w=_f64_buffer(rsw),
            curve_gl=_f64_buffer(gx), curve_gw=_f64_buffer(gw),
            image_masks=_i32_buffer(image_masks), image_signs=_f64_buffer(image_signs),
            eps=eps, leaf=leafsize, eta=eta, build=bool(_build_hmatrix))
        # Keep a permutation-invariant low rule for smooth, well-separated Gram blocks.  The persistent
        # field evaluator retains curved P2 elements directly and does not reuse this approximation.
        rtp_lo, rtw_lo = _outer_tet(far_quad)
        rsp_lo, rsw_lo = _outer_tri(far_quad)
        kw.update(ref_tet_pts_lo=_f64_buffer(rtp_lo), ref_tet_w_lo=_f64_buffer(rtw_lo),
                  ref_tri_pts_lo=_f64_buffer(rsp_lo), ref_tri_w_lo=_f64_buffer(rsw_lo),
                  ho_far_factor=ho_far_factor)
        G = _rp._ChargeGramHMatrix(**kw)
        return _configure_cpp_operator(
            cbk["B"], G, cbk["M_mass"], cbk["M_mass_ngsolve"])
    # INNER subtraction quad (B2 speedup): the subtraction remainder (m_src(y)-m_src(p)) is SMOOTH (the
    # singular part is carried EXACTLY by the analytic PhiTet/TriPotential base), so the inner sum uses a
    # COARSER rule than the outer -> another ~1.5-2x on the O(quad_out^3 * quad_in^3) near entries.  Floor at
    # max(quad-2, p+1); only passed to C++ when iq < quad (else inner = outer).  Validated to hold the same
    # demag accuracy as inner=outer by the nearfar/operator goldens + the uniform-1/3 metric.
    iq = inner_quad if inner_quad is not None else max(quad - 2, p + 1)
    cb = _charge_basis(
        fes, quad, materialize_mass=_materialize_mass,
        internal_interfaces=bool(internal_interfaces))
    B, M_mass, host, kind, expo = cb["B"], cb["M_mass"], cb["host"], cb["kind"], cb["expo"]
    cell_verts, face_verts, n_el = cb["cell_verts"], cb["face_verts"], cb["n_el"]
    # OUTER Gram quadrature: symmetric degree-5 (Keast-15/Dunavant-7) at quad in {3,4}; else product.
    rtp, rtw = _outer_tet(quad)
    rsp, rsw = _outer_tri(quad)
    # BDM1/BDM2 monomial-charge quadrature Gram.  By default, well-separated pairs
    # use the low-order rule while near/self pairs keep subtraction quadrature.
    kw = dict(cell_verts=cell_verts, face_verts=face_verts,
              n_el=n_el, charge_host=_i32_buffer(host), charge_kind=_i32_buffer(kind),
              charge_expo=_i32_buffer(expo),
              ref_tet_pts=_f64_buffer(rtp), ref_tet_w=_f64_buffer(rtw),
              ref_tri_pts=_f64_buffer(rsp), ref_tri_w=_f64_buffer(rsw),
              eps=eps, leaf=leafsize, eta=eta)
    if np.isfinite(ho_far_factor):
        rtp_lo, rtw_lo = _outer_tet(far_quad)
        rsp_lo, rsw_lo = _outer_tri(far_quad)
        kw.update(ref_tet_pts_lo=_f64_buffer(rtp_lo), ref_tet_w_lo=_f64_buffer(rtw_lo),
                  ref_tri_pts_lo=_f64_buffer(rsp_lo), ref_tri_w_lo=_f64_buffer(rsw_lo),
                  ho_far_factor=ho_far_factor)
    if iq < quad:
        rtp_in, rtw_in = _tet_ref(iq)
        rsp_in, rsw_in = _tri_ref(iq)
        kw.update(ref_tet_pts_in=_f64_buffer(rtp_in), ref_tet_w_in=_f64_buffer(rtw_in),
                  ref_tri_pts_in=_f64_buffer(rsp_in), ref_tri_w_in=_f64_buffer(rsw_in))
    if image_masks:
        kw.update(image_masks=_i32_buffer(image_masks), image_signs=_f64_buffer(image_signs))
    kw["build"] = bool(_build_hmatrix)
    G = _rp._ChargeGramHMatrix(**kw)
    return _configure_cpp_operator(B, G, M_mass, cb["M_mass_ngsolve"])


class DemagOperator:
    """ngsolve.bem-style production HDiv-VIM demag operator.

    Construct from an order-1/order-2 HDiv FESpace, or straight TET/HEX
    order-0 RT0 for the explicit broken-interface material-topology path;
    ``.mat`` is the analytic C++
    charge-Gram operator ``N = B^T G B``.  The same operator is used by
    :func:`radia.vim.Solve`, so the diagnostic and material paths cannot drift
    into separate research backends.  The caller wraps construction and
    ``DemagFactor`` in ``with TaskManager():``.
    """

    def __init__(self, fes, intorder=None, eps=1e-7, leafsize=16, eta=2.0,
                 far_quad=3, ho_far_factor=2.0, inner_quad=None,
                 curve_order=None, curve_gauss=8,
                 internal_interfaces=False, image_masks=None,
                 image_signs=None):
        p = int(fes.globalorder)
        vtypes = _volume_vertex_counts(fes.mesh)
        validate_hdiv_configuration(fes.mesh.dim, vtypes, p, fes.mesh.GetCurveOrder())
        self.space = fes
        # AUTO-MATCH the Gram curve order to the MESH geometry order (mesh.GetCurveOrder()).  A STRAIGHT Gram on
        # a CURVED mesh (where B/M_mass are NGSolve curved integrals) is geometry-inconsistent and the demag
        # factor DRIFTS with geometry order (sphere: straight-Gram 0.336/0.308/0.279 at curve 1/2/3; the matched
        # curved Gram restores ~1/3 -- 0.338 at curve 2).  curve_order=None => auto from GetCurveOrder(); pass an
        # explicit int to override (curve_order=0 forces the STRAIGHT Gram, e.g. a deliberate flat-Gram probe).
        if curve_order is None:
            _k = fes.mesh.GetCurveOrder()
            curve_order = _k if _k >= 2 else None
        elif curve_order == 0:
            curve_order = None
        self.curve_order = curve_order
        self.internal_interfaces = bool(internal_interfaces)
        self.image_masks = tuple(() if image_masks is None else image_masks)
        self.image_signs = tuple(() if image_signs is None else image_signs)
        self._B, self._G, self._Mmass = build_charge_gram(
            fes, intorder=intorder, eps=eps, leafsize=leafsize, eta=eta,
            far_quad=far_quad, ho_far_factor=ho_far_factor, inner_quad=inner_quad,
            curve_order=curve_order, curve_gauss=curve_gauss,
            internal_interfaces=self.internal_interfaces,
            image_masks=self.image_masks, image_signs=self.image_signs)
        self.mat = self._G.demag_matrix()

    @property
    def ndof(self):
        return self.space.ndof

    def DemagFactor(self, M_cf):
        """Rayleigh quotient (the demag factor) for a magnetization CoefficientFunction M_cf:
        <c, G c> / <m, M_mass m>, c = B m, m = the HDiv projection of M_cf.  ~1/3 for a sphere/cube."""
        gfu = ng.GridFunction(self.space)
        gfu.Set(M_cf)
        demag = gfu.vec.CreateVector()
        self.mat.Mult(gfu.vec, demag)
        numerator = float(ng.InnerProduct(gfu.vec, demag))
        denominator = float(ng.Integrate(ng.InnerProduct(gfu, gfu), self.space.mesh))
        return numerator / denominator
