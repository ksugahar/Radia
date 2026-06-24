"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r), and
`hdiv_demag_solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the FEEC/HDiv counterpart to the moment-yano MSC hex/wedge/pyramid soft-iron demag in `rad.Solve`.
Both modes take an ARBITRARY applied field `H_ext` (any NGSolve CoefficientFunction -- e.g. a coil's
Biot-Savart field `rad.RadiaField(coil,'h')`, the C-type electromagnet driver) and return per-element M.

## Formulation (verified-first, 2026-06-15)
ONE projected weak form everywhere -- the magnetization M is the RT0 primary, the constitutive law
M = M(H) is imposed in the L2 sense (M_mass m = INT M(H).v dx), and H = h_ext - M_mass^-1 N m is the
weak total field (N = B^T G B, h_ext = H_ext L2-projected onto RT0 order 0).  LINEAR soft iron is the
CONSTANT-chi special case M = chi H, giving the form-1 system

    (M_mass + M_chi M_mass^-1 N) m = M_chi h_ext ,   M_chi = INT chi(x) u.v dx   (the chi-weighted mass),

solved +N (the -N system is mu_r-independent but NON-physical -- wrong sign).  For a SINGLE region
(uniform chi) this is identical bit-for-bit to ((1/chi) M_mass + N) m = M_mass h_ext (verified vs the
analytic sphere to 1.4e-4); for PER-REGION chi it is the consistent generalization, and it makes a
linear region agree with its nonlinear-table equivalent (the 1/chi-weighted form differs at material
interfaces by O(h)).  NONLINEAR iron solves the same projected residual by Newton (below).

NONLINEAR -- a **damped matrix-free Newton** on the constitutive residual.  F(m) = M_mass m - b_M(H),
H = h_ext - M_mass^-1 N m, b_M(H) = INT M(H).v dx (the RT0 projection of M(H)); consistent tensor
Jacobian J = M_mass + T D_op (T = INT (dM/dH) u.v dx the tensor-tangent mass, D_op = M_mass^-1 N),
applied MATRIX-FREE J v = M_mass v + T (M_mass^-1 (N v)); each Newton step is ONE M_mass^-1-GMRES on J
+ an Armijo line search, after a scalar-chi LINEAR warmstart into the basin.  WHY Newton, not the
earlier Anderson-Hantila fixed point: Hantila/Picard's contraction
rho = (chi_max - chi_min)/(chi_max + chi_min) -> 1 as the iron saturates (chi_max ~ 1e4 unsaturated,
chi_min ~ 10 at the pole), so it STALLS on real silicon steel (measured: ~1e-2 residual after 300
iters, NOT converged) -- Newton's quadratic step is immune to the chi-range.  VERIFIED: the real CEFC
Si-steel C-type at 3000 AT converges relF 0.55 -> 8e-7 in ~24 iters, gap B within ~1% of the FEM truth
(Anderson-Hantila could not solve it at all).  N = B^T G B stays SPD-PSD (G the SPD Coulomb Gram), the
de-Rham loops sit in ker(N) and are carried by M_mass (no loop-star), and the scalable C++ charge-Gram
H-matvec is the only O(N log N) cost.

## Preconditioner -- the FULL HDiv mass inverse (mesh + mu_r robust)
Every +N / constant-LHS solve is GMRES preconditioned with `M_mass^{-1}` (a one-time sparse LU of the
local RT0 HDiv mass), which DEFLATES the de-Rham loops (ker B) -> the iteration count stays bounded as
the mesh refines AND as mu_r grows (the de-Rham structural advantage: no loop-star, distortion-robust).
A Jacobi diagonal is NOT mesh-robust (the +N count blew past 5000 iters at 7224 faces -- measured).
GMRES (not CG/MINRES) absorbs the residual ACA asymmetry of the charge-Gram H-matvec.

KELVIN-LESS: the 1/r charge Gram IS the open boundary (a volume integral method like MMM/MSC); only
the iron is meshed -- no air box / Kelvin needed.  The NONLINEAR path uses the analytic charge Gram
(scalable `_ChargeGramHMatrix` at tight gram_eps), REQUIRED for div M != 0 (non-uniform M) bodies.

## Scope (M1)
Per-region soft iron, LINEAR (`mu_r` scalar or `{material: mu_r}` dict) AND NONLINEAR (`bh_table` one
[[H,B]] table or `{material: [[H,B]]}` dict).  N = B^T G B is geometry-only, so multi-grade iron enters
ONLY through the (1/chi)-weighted HDiv mass (linear) / the per-element constitutive law (nonlinear).
Mixed PM+iron (fixed-M source regions) + the 165k-DOF-scale preconditioner + the M0 parity gate are the
remaining productionization steps (docs/hdiv_vim/PRODUCTIONIZATION.md).  Until they land, moment-yano MSC
stays the `rad.Solve` default demag backend (`radia.set_demag_backend`); this entry does not touch it.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- this library helper does NOT
open a TaskManager; the caller wraps the call in `with ng.TaskManager():`.
"""
import inspect
from math import pi as _PI

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import ngsolve as ng

import radia._radia_pybind as _rp
from . import _core as _tet
from ._nonlinear import _bh_table_funcs, _table_tensor_tangent, _table_tensor_tangent_multi

_MU0 = 4e-7 * _PI
# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect once.
_GMRES_TOL = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, bh_table=None, pm_M=None,
                     image=None, gram_eps=1e-10, leaf=32, eta=2.0, near_factor=1e30, tol=1e-8,
                     maxit=4000, gmres_restart=400, nl_maxit=300, nl_tol=1e-6, anderson_window=6):
    """HDiv-type VIM soft-iron demag solve (the +N physical material system).

    Soft-iron spec (EXACTLY ONE, unless every region is a permanent magnet -> both may be omitted):
      mu_r     : float > 1 -> LINEAR isotropic soft iron (ONE region), OR a dict {material_name: mu_r}
                 for PER-REGION linear soft iron (multiple iron grades; each mu_r > 1).  N = B^T G B is
                 geometry-only, so per-region enters ONLY through the chi-weighted HDiv mass (form 1).
      bh_table : [[H,B], ..] (A/m, T; the MatSatIsoTab data) -> NONLINEAR isotropic soft iron (ONE
                 region), OR a dict {material_name: [[H,B], ..]} for PER-REGION nonlinear soft iron.
                 N = B^T G B is geometry-only, so per-region nonlinear enters ONLY through the per-element
                 constitutive law + warmstart.
      pm_M     : dict {material_name: [Mx,My,Mz]} (A/m) -> PERMANENT-MAGNET (fixed-M source) regions
                 mixed alongside soft iron (or alone).  A PM region is the M = M_pm (tangent 0) special
                 case of the same projected statement: it contributes a source b_pm = INT M_pm.v dx and
                 its field reaches the iron through the full-m demag, while its own M stays pinned to
                 M_pm.  With pm_M, the soft iron is given as EITHER mu_r (linear) OR bh_table (nonlinear)
                 and covers the NON-PM regions (a scalar applies to every non-PM region; a dict names
                 them).  PM directly TOUCHING soft iron is rejected (conforming RT0 cannot represent the
                 PM-iron magnetization discontinuity) -- separate them with an air gap.
    H_ext      : NGSolve CoefficientFunction, the applied field (A/m) -- uniform, analytic, or a coil's
                 Biot-Savart field rad.RadiaField(coil,'h').  REQUIRED.
    Returns dict: M (n_el,3) per-element magnetization, M_avg (3,), iters, demag (Rayleigh factor),
    ndof, n_el, n_charge, nonlinear(bool).  The caller must open `with ng.TaskManager():`.
    """
    if H_ext is None:
        raise ValueError("hdiv_demag_solve: H_ext (applied-field CoefficientFunction) is required")
    if pm_M is None:
        if (mu_r is None) == (bh_table is None):
            raise ValueError("hdiv_demag_solve: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")
    else:
        if (mu_r is not None) and (bh_table is not None):
            raise ValueError("hdiv_demag_solve: with pm_M, give the iron as EITHER mu_r (linear) OR "
                             "bh_table (nonlinear), not both")

    # ---- demag operator N: ALWAYS the scalable C++ charge-Gram H-matrix (_ChargeGramHMatrix) ----
    # The C++ kernel handles tet (cell_verts/face_verts) AND hex/wedge (the polytope triangle-soup ctor),
    # and folds IMA mirror symmetry as image charges (image_masks/image_signs) so the reduced (1/2,1/4,1/8)
    # model reproduces the full one -- validated == the dense IMA Gram to ~1e-10.  The `scalable` knob + the
    # dense O(N^2) Python Gram path were REMOVED (the dense path was ~70x slower at ~1000 elements;
    # No-Fallbacks: one supported C++ path).
    all_tet = all(len(el.vertices) == 4 for el in mesh.Elements(ng.VOL))
    image_masks, image_signs = [], []
    if image is not None:
        for axes, sign in _tet.image_group(_tet.parse_image_string(image)):
            image_masks.append(int(sum(1 << a for a in axes)))
            image_signs.append(float(sign))
    d = _tet.build_demag(mesh)
    Mm, B = d["M_mass"], d["B_csr"]
    if all_tet:
        H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                   n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta,
                                   near_factor=near_factor, image_masks=image_masks, image_signs=image_signs)
    else:
        # HEX/WEDGE: the polytope triangle-soup charge Gram (build_demag emits d["poly"] for non-tet).
        p = d["poly"]
        H = _rp._ChargeGramHMatrix(
            cell_tris=list(p["cell_tris"]), cell_troff=list(p["cell_troff"]),
            cell_cent=list(p["cell_cent"]), cell_meas=list(p["cell_meas"]),
            face_tris=list(p["face_tris"]), face_troff=list(p["face_troff"]),
            face_cent=list(p["face_cent"]), face_meas=list(p["face_meas"]),
            n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta, near_factor=near_factor,
            image_masks=image_masks, image_signs=image_signs)

    def N_apply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    n_face, n_el = int(d["ndof"]), int(d["n_el"])
    mu = d["m_unit"]
    Mfac = spla.splu(sp.csc_matrix(Mm))
    Mprec = spla.LinearOperator((n_face, n_face), matvec=lambda v: Mfac.solve(np.asarray(v, float)))

    # ---- applied field projected onto RT0 ----
    fes = ng.HDiv(mesh, order=0)
    gfHext = ng.GridFunction(fes); gfHext.Set(H_ext)
    h_ext = gfHext.vec.FV().NumPy().copy()
    denom = float(mu @ np.asarray(Mm @ mu).ravel())
    D = float((mu @ N_apply(mu)) / denom)

    if pm_M is not None and bh_table is not None:
        m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec,
                                    n_face, fes, tol, gmres_restart, nl_maxit, nl_tol, anderson_window,
                                    pm_M=pm_M)
    elif pm_M is not None:
        M_chi, b_pm = _build_mixed_pm(mesh, fes, mu_r, pm_M, Mm, n_face)
        m, iters = _solve_linear(Mm, M_chi, N_apply, Mfac, Mprec, n_face, h_ext, tol, maxit,
                                 gmres_restart, b_extra=b_pm)
    elif mu_r is not None:
        M_chi = _build_chi_mass(mesh, fes, mu_r, Mm, n_face)
        m, iters = _solve_linear(Mm, M_chi, N_apply, Mfac, Mprec, n_face, h_ext, tol, maxit, gmres_restart)
    else:
        m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec,
                                    n_face, fes, tol, gmres_restart, nl_maxit, nl_tol, anderson_window)

    # ---- per-element M (VectorL2(0) projection; component-major DOFs -> (n_el,3)) + averages ----
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol

    return dict(M=M_el, M_avg=M_avg, iters=int(iters), demag=D, ndof=n_face, n_el=n_el,
                n_charge=int(d["n_charge"]), nonlinear=(bh_table is not None))


def _build_chi_mass(mesh, fes, mu_r, Mm, n_face):
    """The chi-weighted HDiv mass M_chi = INT chi(x) u.v dx for the form-1 linear system
    A = M_mass + M_chi M_mass^-1 N (the projected M = chi H statement; see the module docstring).
    `mu_r` is either a scalar (one isotropic soft-iron region -> M_chi = chi M_mass) or a dict
    {material: mu_r} (per-region soft iron -> a region-wise coefficient).  N = B^T G B is geometry-only,
    so per-region soft iron enters ONLY through this mass.  Fail-loud (No-Fallbacks): every mesh
    material must be specified, each mu_r > 1."""
    if isinstance(mu_r, dict):
        mats = list(mesh.GetMaterials())
        missing = sorted(set(mats) - set(mu_r))
        if missing:
            raise ValueError("hdiv_demag_solve: mu_r dict missing region(s) %s; mesh materials are %s"
                             % (missing, mats))
        bad = {r: mu_r[r] for r in mats if float(mu_r[r]) <= 1.0}
        if bad:
            raise ValueError("hdiv_demag_solve: every region mu_r must be > 1 (got %s)" % bad)
        chi_cf = mesh.MaterialCF({r: float(mu_r[r]) - 1.0 for r in mats})
        u, v = fes.TnT()
        a = ng.BilinearForm(fes); a += chi_cf * u * v * ng.dx; a.Assemble()
        r, c, val = a.mat.COO()
        return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(n_face, n_face))
    chi = float(mu_r) - 1.0
    if chi <= 0.0:
        raise ValueError("hdiv_demag_solve: mu_r must be > 1 (got %r)" % (mu_r,))
    return chi * Mm


def _build_mixed_pm(mesh, fes, mu_r, pm_M, Mm, n_face):
    """Permanent-magnet (fixed-M) + LINEAR soft-iron material build for the form-1 mixed system
    (M_mass + M_chi M_mass^-1 N) m = M_chi h_ext + b_pm.  A PM region is the M = M_pm (no H-response,
    chi = 0) special case of the projected statement: it contributes the source load b_pm = INT M_pm.v dx
    and its field reaches the iron through the full-m demag, while its own M stays pinned to M_pm by the
    projection (chi = 0 there -> M_mass m = b_pm).  `pm_M` is {material: [Mx,My,Mz]} (A/m); the soft-iron
    spec covers the NON-PM regions (scalar mu_r -> every non-PM region; dict mu_r -> the named iron
    regions; mu_r = None -> pure PM, no iron).  Returns (M_chi sparse mass, b_pm source vector).
    Fail-loud (No-Fallbacks): every material is classified exactly once (iron XOR PM), PM regions exist,
    each iron mu_r > 1."""
    mats = list(mesh.GetMaterials())
    pm_regions = set(pm_M)
    missing_pm = sorted(pm_regions - set(mats))
    if missing_pm:
        raise ValueError("hdiv_demag_solve: pm_M region(s) %s not in mesh materials %s" % (missing_pm, mats))
    if mu_r is None:
        iron_mu = {}
    elif isinstance(mu_r, dict):
        iron_mu = {r: float(mu_r[r]) for r in mu_r}
    else:
        iron_mu = {r: float(mu_r) for r in mats if r not in pm_regions}
    iron_regions = set(iron_mu)
    overlap = sorted(iron_regions & pm_regions)
    if overlap:
        raise ValueError("hdiv_demag_solve: region(s) %s are both iron (mu_r) and PM (pm_M)" % overlap)
    unspec = sorted(set(mats) - iron_regions - pm_regions)
    if unspec:
        raise ValueError("hdiv_demag_solve: region(s) %s are neither iron (mu_r) nor PM (pm_M); "
                         "mesh materials are %s" % (unspec, mats))
    bad = {r: iron_mu[r] for r in iron_regions if iron_mu[r] <= 1.0}
    if bad:
        raise ValueError("hdiv_demag_solve: every iron mu_r must be > 1 (got %s)" % bad)
    _check_pm_iron_not_touching(mesh, fes, pm_regions, iron_regions)
    u, v = fes.TnT()
    chi_cf = mesh.MaterialCF({r: (iron_mu[r] - 1.0 if r in iron_regions else 0.0) for r in mats})
    a = ng.BilinearForm(fes); a += chi_cf * u * v * ng.dx; a.Assemble()
    rr, cc, val = a.mat.COO()
    M_chi = sp.csr_matrix((np.array(val), (np.array(rr), np.array(cc))), shape=(n_face, n_face))
    pm_cf = mesh.MaterialCF({r: ng.CoefficientFunction(tuple(pm_M[r])) if r in pm_regions
                             else ng.CoefficientFunction((0.0, 0.0, 0.0)) for r in mats})
    lf = ng.LinearForm(fes); lf += pm_cf * v * ng.dx; lf.Assemble()
    b_pm = lf.vec.FV().NumPy().copy()
    return M_chi, b_pm


def _check_pm_iron_not_touching(mesh, fes, pm_regions, iron_regions):
    """Fail-loud if a PM region shares an RT0 facet with a soft-iron region.  The magnetization is
    DISCONTINUOUS across a PM<->iron boundary (M_pm != chi H there), but the conforming HDiv (RT0,
    order 0) field forces the normal component continuous on every shared facet, so a directly-touching
    PM-iron interface produces a spurious interface artifact that DIVERGES under refinement (measured:
    ~20% external-field error vs rad.Solve, growing with mesh density).  A non-touching PM+iron
    (an air gap between them -- the usual magnetic-circuit / PM-motor case) shares no facets and matches
    rad.Solve to ~1e-3, converging.  Detecting a shared PM-iron facet -> raise (the broken-RT0 /
    interface-DG treatment of a touching PM-iron boundary is the next productionization step)."""
    cls = {}                                   # HDiv order-0 DOF (== facet) -> set of {'pm','iron'}
    for el in mesh.Elements(ng.VOL):
        c = "pm" if el.mat in pm_regions else ("iron" if el.mat in iron_regions else None)
        if c is None:
            continue
        for dof in fes.GetDofNrs(el):
            if dof >= 0:
                cls.setdefault(dof, set()).add(c)
    n_shared = sum(1 for s in cls.values() if "pm" in s and "iron" in s)
    if n_shared:
        raise NotImplementedError(
            "hdiv_demag_solve: a permanent-magnet region directly TOUCHES a soft-iron region "
            "(%d shared RT0 facets).  The conforming HDiv field cannot represent the magnetization "
            "discontinuity across a PM-iron boundary, so the result would DIVERGE under refinement "
            "(~20%% external-field error vs rad.Solve).  Separate the PM and iron with an air gap "
            "(verified to match rad.Solve to ~1e-3); the touching-interface formulation is the next "
            "productionization step." % n_shared)


def _build_nl_constit(mesh, fes, bh_table, Mm, n_face, pm_M=None):
    """Build the NONLINEAR constitutive callback + the zero-field-chi (form-1) warmstart mass M_chi0 +
    the permanent-magnet source b_pm.

    `bh_table` is either ONE [[H,B]] table (a single nonlinear iron grade for every NON-PM region) or a
    dict {material_name: [[H,B]]} (PER-REGION nonlinear iron grades).  N = B^T G B is geometry-only, so
    per-region nonlinear enters ONLY through (a) the per-element constitutive law (each element uses its
    region's PCHIP BH table via `_table_tensor_tangent_multi`) and (b) the chi0-weighted warmstart mass
    M_chi0 (region-wise zero-field susceptibility chi0 = B'(0)/mu0 - 1).

    `pm_M` (optional) = {material: [Mx,My,Mz]} fixed-magnet regions mixed with the nonlinear iron: each
    PM region is the M = M_pm (tangent 0) special case, overriding the constitutive law via a MaterialCF
    and contributing the source b_pm = INT M_pm.v dx (the same projected statement as the linear PM path,
    now inside Newton).  PM directly touching iron is rejected (the conforming RT0 limitation).

    Returns (constit_fn(gfH, Id) -> (M target CF, tensor-tangent CF), M_chi0 sparse/scaled mass, b_pm
    source vector or None).  Fail-loud (No-Fallbacks): every material classified, each table 2-column."""
    mats = list(mesh.GetMaterials())
    pm_regions = set(pm_M) if pm_M else set()
    if pm_M:
        missing_pm = sorted(pm_regions - set(mats))
        if missing_pm:
            raise ValueError("hdiv_demag_solve: pm_M region(s) %s not in mesh materials %s" % (missing_pm, mats))

    # ---- iron regions + per-region/single BH constitutive ----
    if isinstance(bh_table, dict):
        iron_regions = set(bh_table)
        if not pm_M:
            missing = sorted(set(mats) - iron_regions)
            if missing:
                raise ValueError("hdiv_demag_solve: bh_table dict missing region(s) %s; mesh materials are %s"
                                 % (missing, mats))
        region_names = list(bh_table)
        name_to_ridx = {nm: k for k, nm in enumerate(region_names)}
        region_funcs = []
        chi0_by_name = {}
        for nm in region_names:
            arr = np.asarray(bh_table[nm], float)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError("hdiv_demag_solve: bh_table[%r] must be [[H,B], ...] (A/m, T)" % (nm,))
            _Mof, Bpch, Bder, Hmax, Mmax = _bh_table_funcs(arr[:, 0], arr[:, 1])
            region_funcs.append((Bpch, Bder, Hmax, Mmax))
            chi0_by_name[nm] = max(float(Bder(0.0)) / _MU0 - 1.0, 1.0)
        # PM elements map to iron index 0 (their iron values are computed but OVERRIDDEN by the PM
        # MaterialCF below); iron elements map to their grade.
        elem_region = np.array([name_to_ridx.get(mesh[ng.ElementId(ng.VOL, i)].mat, 0)
                                for i in range(mesh.ne)], dtype=int)

        def _iron_tangent(gfH, Id):
            return _table_tensor_tangent_multi(gfH, mesh, region_funcs, elem_region, Id)

        def _chi0_of(r):
            return chi0_by_name[r]
    else:
        iron_regions = set(mats) - pm_regions
        arr = np.asarray(bh_table, float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("hdiv_demag_solve: bh_table must be [[H,B], ...] (A/m, T)")
        _Mof, Bpch, Bder, Hmax, Mmax = _bh_table_funcs(arr[:, 0], arr[:, 1])
        chi0 = max(float(Bder(0.0)) / _MU0 - 1.0, 1.0)

        def _iron_tangent(gfH, Id):
            return _table_tensor_tangent(gfH, mesh, Bpch, Bder, Hmax, Mmax, Id)

        def _chi0_of(r):
            return chi0

    if not pm_M:
        # iron-only: M_chi0 reuses the linear chi-weighted mass (every material is iron)
        if isinstance(bh_table, dict):
            M_chi0 = _build_chi_mass(mesh, fes, {nm: chi0_by_name[nm] + 1.0 for nm in mats}, Mm, n_face)
        else:
            M_chi0 = chi0 * Mm
        return _iron_tangent, M_chi0, None

    # ---- PM + nonlinear iron: classify XOR, reject touching, build mixed constitutive + b_pm ----
    overlap = sorted(iron_regions & pm_regions)
    if overlap:
        raise ValueError("hdiv_demag_solve: region(s) %s are both iron (bh_table) and PM (pm_M)" % overlap)
    unspec = sorted(set(mats) - iron_regions - pm_regions)
    if unspec:
        raise ValueError("hdiv_demag_solve: region(s) %s are neither iron (bh_table) nor PM (pm_M); "
                         "mesh materials are %s" % (unspec, mats))
    _check_pm_iron_not_touching(mesh, fes, pm_regions, iron_regions)
    _ZERO3 = ng.CoefficientFunction((0.0,) * 9, dims=(3, 3))

    def constit_fn(gfH, Id):
        M_iron, tang_iron = _iron_tangent(gfH, Id)
        M_target = mesh.MaterialCF({**{r: M_iron for r in iron_regions},
                                    **{r: ng.CoefficientFunction(tuple(pm_M[r])) for r in pm_regions}})
        tang_target = mesh.MaterialCF({**{r: tang_iron for r in iron_regions},
                                       **{r: _ZERO3 for r in pm_regions}})
        return M_target, tang_target

    # warmstart mass: chi0 on iron, 0 on PM
    u, v = fes.TnT()
    chi0_cf = mesh.MaterialCF({**{r: _chi0_of(r) for r in iron_regions}, **{r: 0.0 for r in pm_regions}})
    am = ng.BilinearForm(fes); am += chi0_cf * u * v * ng.dx; am.Assemble()
    rr, cc, vv = am.mat.COO()
    M_chi0 = sp.csr_matrix((np.array(vv), (np.array(rr), np.array(cc))), shape=(n_face, n_face))
    # PM source b_pm = INT M_pm . v dx
    pm_cf = mesh.MaterialCF({**{r: ng.CoefficientFunction((0.0, 0.0, 0.0)) for r in iron_regions},
                             **{r: ng.CoefficientFunction(tuple(pm_M[r])) for r in pm_regions}})
    lf = ng.LinearForm(fes); lf += pm_cf * v * ng.dx; lf.Assemble()
    b_pm = lf.vec.FV().NumPy().copy()
    return constit_fn, M_chi0, b_pm


def _solve_linear(Mm, M_chi, N_apply, Mfac, Mprec, n_face, h_ext, tol, maxit, gmres_restart, b_extra=None):
    """Form-1 linear solve: (M_mass + M_chi M_mass^-1 N) m = M_chi h_ext (+ b_extra), GMRES + M_mass^-1
    precond.  This is the constant-chi special case of the projected M = chi H statement (so a linear
    region and its nonlinear-table equivalent agree).  For uniform chi it is identical to the
    (1/chi)M_mass + N system (verified bit-for-bit); for per-region chi it is the consistent
    generalization.  `b_extra` (the permanent-magnet source load b_pm = INT M_pm.v dx) is added to the
    RHS when permanent magnets are mixed in."""
    rhs = np.asarray(M_chi @ h_ext).ravel()
    if b_extra is not None:
        rhs = rhs + np.asarray(b_extra).ravel()

    def A_apply(v):
        v = np.asarray(v, float)
        return np.asarray(Mm @ v).ravel() + np.asarray(M_chi @ Mfac.solve(N_apply(v))).ravel()

    A = spla.LinearOperator((n_face, n_face), matvec=A_apply)
    it = {"n": 0}
    cycles = max(2, int(np.ceil(maxit / gmres_restart)))
    m, info = spla.gmres(A, rhs, M=Mprec, restart=int(gmres_restart), maxiter=cycles,
                         callback=lambda _x: it.__setitem__("n", it["n"] + 1),
                         callback_type="pr_norm", **{_GMRES_TOL: tol})
    if info != 0:
        raise RuntimeError("hdiv_demag_solve (linear): form-1 GMRES did not converge (info=%d, "
                           "n_face=%d, iters=%d)" % (info, n_face, it["n"]))
    return m, it["n"]


def _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec,
                     n_face, fes, gmres_tol, gmres_restart, nl_maxit, nl_tol, anderson_window, pm_M=None):
    """Damped matrix-free NEWTON on the constitutive residual (replaces the Anderson-Hantila fixed
    point, which STALLS on real silicon steel -- the Hantila/Picard contraction
    rho = (chi_max-chi_min)/(chi_max+chi_min) -> 1 as the iron saturates, so even Anderson cannot
    recover it within a bounded iter count).  Newton converges where Hantila cannot: VERIFIED on the
    real CEFC Si-steel C-type at 3000 AT (relF 0.55 -> 8e-7 in ~24 iters, gap B within ~1% of the FEM
    truth) where Anderson-Hantila stalled at ~1e-2 after 300.  `anderson_window` is now unused (kept in
    the signature for call-site stability).

    Residual    F(m) = M_mass m - b_M(H),  H = h_ext - M_mass^-1 N m,  b_M(H) = INT M(H).v dx (the RT0
                L2 projection of the constitutive M).
    Consistent TENSOR Jacobian  J = M_mass + T D_op,  T = INT (dM/dH) u.v dx (the tensor-tangent mass),
                D_op = M_mass^-1 N -- applied MATRIX-FREE:  J v = M_mass v + T (M_mass^-1 (N v)).
    Each Newton step is ONE M_mass^-1-preconditioned GMRES on J + an Armijo backtracking line search.
    A scalar-chi LINEAR warmstart  ((1/chi0) M_mass + N) m = M_mass h_ext  lands inside the Newton basin
    (the unsaturated chi0 response).  The tensor tangent (`_table_tensor_tangent`) + matrix-free N reuse
    the EXACT pieces the Hantila path used -- no new operator, fully scalable (the C++ charge-Gram
    H-matvec is the only O(N log N) cost).  Fail-loud on non-convergence (CLAUDE.md No-Fallbacks)."""
    del anderson_window                      # unused by Newton (signature kept for the caller)
    constit_fn, M_chi0, b_pm = _build_nl_constit(mesh, fes, bh_table, Mm, n_face, pm_M=pm_M)
    Id = ng.Id(3); vf = fes.TestFunction(); uf = fes.TrialFunction(); gfH = ng.GridFunction(fes)

    def _constit(m):                         # (M(H) CF, tensor-tangent CF) at H = h_ext - M_mass^-1 N m
        gfH.vec.FV().NumPy()[:] = h_ext - Mfac.solve(N_apply(m))
        return constit_fn(gfH, Id)

    def _bM(M_cf):                           # RT0 L2 projection load INT M(H).v dx
        lf = ng.LinearForm(fes); lf += M_cf * vf * ng.dx; lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def _Fnorm(m):
        M_cf, _ = _constit(m)
        return float(np.linalg.norm(np.asarray(Mm @ m).ravel() - _bM(M_cf)))

    # zero-field-chi LINEAR warmstart -> the Newton basin (the unsaturated chi0 response, per region for a
    # bh_table dict).  Form-1 warmstart (M_mass + M_chi0 M_mass^-1 N) m = M_chi0 h_ext, where M_chi0 is the
    # chi0-weighted HDiv mass -- the SAME projected system the LINEAR path solves (so linear == the
    # constant-chi nonlinear limit), just at the zero-field susceptibility chi0.
    rhs0 = np.asarray(M_chi0 @ h_ext).ravel()
    if b_pm is not None:                     # PM source enters the warmstart RHS too (mixed PM+iron)
        rhs0 = rhs0 + np.asarray(b_pm).ravel()
    A0 = spla.LinearOperator((n_face, n_face), matvec=lambda v:
                             np.asarray(Mm @ np.asarray(v, float)).ravel()
                             + np.asarray(M_chi0 @ Mfac.solve(N_apply(v))).ravel())
    m, _ = spla.gmres(A0, rhs0, M=Mprec, restart=int(gmres_restart), maxiter=4, **{_GMRES_TOL: 1e-6})

    converged = False; nit = 0; rel = float("inf")
    for it in range(nl_maxit):
        nit = it + 1
        M_cf, tang = _constit(m)
        F = np.asarray(Mm @ m).ravel() - _bM(M_cf)
        nF = float(np.linalg.norm(F))
        scale = float(np.linalg.norm(np.asarray(Mm @ m).ravel())) + 1e-30
        rel = nF / scale
        if rel < nl_tol:
            converged = True; break
        T = ng.BilinearForm(fes); T += ng.InnerProduct(tang * uf, vf) * ng.dx; T.Assemble()
        r, c, v = T.mat.COO()
        Tcsr = sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(n_face, n_face))

        def _Jv(x, _Tcsr=Tcsr):              # J v = M_mass v + T (M_mass^-1 (N v))   (matrix-free)
            x = np.asarray(x, float)
            return np.asarray(Mm @ x).ravel() + np.asarray(_Tcsr @ Mfac.solve(N_apply(x))).ravel()

        Jop = spla.LinearOperator((n_face, n_face), matvec=_Jv)
        # Eisenstat-Walker inexact-Newton forcing: solve the Newton step only as tightly as the current
        # outer residual warrants (loose when far, tight when near) -- the inner J-GMRES is the dominant
        # cost (~130 iters/step, M_mass^-1 is a weak preconditioner for J = M_mass + T D_op), so a loose
        # early tol cuts wasted inner work without changing the (tangent-limited, ~0.65-linear) outer rate.
        inner_tol = float(min(1e-2, max(gmres_tol, 0.1 * rel)))
        dm, _info = spla.gmres(Jop, -F, M=Mprec, restart=int(gmres_restart), maxiter=4, **{_GMRES_TOL: inner_tol})
        lam = 1.0                            # Armijo backtracking line search (globalises Newton)
        while lam > 1e-6 and _Fnorm(m + lam * dm) >= nF:
            lam *= 0.5
        m = m + lam * dm
    if not converged:
        raise RuntimeError("hdiv_demag_solve (nonlinear Newton): did NOT converge -- rel=%.2e > "
                           "nl_tol=%.1e after %d iters (returning M would be a silent wrong result)"
                           % (rel, nl_tol, nit))
    return m, nit
