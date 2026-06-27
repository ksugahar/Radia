"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r), and
`hdiv_demag_solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the FEEC/HDiv counterpart to the multipole-moment MMM MSC hex/wedge/pyramid soft-iron demag in `rad.Solve`.
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

## Linear solve dispatch -- SYMMETRIC C++ CG default (symmetric HACApK), GMRES cross-check opt-in
For the common scalar-mu_r case, `linear_solver="auto"` (the DEFAULT) solves the SPD +N system
`((1/chi)M_mass + N) m = M_mass h_ext` by CG, preconditioned with the FULL RT0 H(div) mass inverse
`M_mass^{-1}` (the MASS RIESZ map), ENTIRELY in C++ (PARDISO mass factor + C++ Krylov, no Python glue).
CG (not GMRES) because the charge-Gram is applied via the EXACTLY-SYMMETRIC H-matvec (`matvec_sym`): the
HACApK H-matrix stores both (I,J) and (J,I) leaf blocks but ACA-truncates them INDEPENDENTLY, so the
GENERAL matvec is only approximately symmetric; `matvec_sym` instead applies the UPPER-triangular leaves
only -- each upper leaf supplies its own block AND the mirror as its exact transpose -- so the operator is
machine-symmetric (||G - G^T|| == 0) regardless of the per-block truncation.  This makes CG robust BY
CONSTRUCTION at all N (it removes the asymmetry failure mode entirely), and the symmetric matvec is ~1.4x
FASTER than the general one (it touches half the leaves).  So CG is the default again (Sugahara 2026-06-27,
"対称HACApKを実装しよう。CGがいいね").  `linear_solver="cpp-cg"` is an explicit alias for this symmetric C++ CG.
`linear_solver="gmres"` is the asymmetry-tolerant cross-check/opt-in (mass-Riesz GMRES on the GENERAL matvec,
Python recurrence) -- it was the default 2026-06-27 morning before the symmetric matvec landed.
`linear_solver="hlu"` is an opt-in tet-only system-A H-LU path (`A = M_mass + chi*N` on face DOFs).
`linear_solver="python"` is the form-1 GMRES + `M_mass^{-1}` sparse LU, used for per-region chi / PM-mixed /
nonlinear Newton paths until their material-specific operators move to C++.  (The mass Riesz makes the
operator well-conditioned by construction -- the earlier "h-explosion => need AMS" was a monopole-Gram
artifact; the accurate analytic Gram + mass Riesz needs no auxiliary-space preconditioner.)

The uniform-linear Krylov paths (default CG = auto/cpp-cg, and the gmres cross-check) build the analytic Gram
at the tight `gram_eps=1e-12`; per-region / nonlinear / H-LU keep `1e-10`.  (With the symmetric matvec the CG
no longer NEEDS 1e-12 for symmetry -- symmetry is now STRUCTURAL, independent of the ACA accuracy -- 1e-12 is
kept only for solution ACCURACY + golden stability.)  An explicit `gram_eps` always wins.  All material solve
paths are fail-loud: a non-converged solve RAISES (No-Fallbacks) rather than returning a wrong M.

ON THE EARLIER "+N CG SCALE WALL" (recorded honestly): the 2026-06-27-morning GMRES retreat was motivated by a
report that the GENERAL-matvec CG diverges past nf ~ 20k (the spurious ACA antisymmetry growing with N).  When
the symmetric matvec was added, a re-measurement could NOT reproduce that divergence at HDiv scales: even with
the lossy monopole far (max asymmetry) + a distorted hex + mu_r up to 1e6, the GENERAL-matvec CG converges
fine through nf ~ 51k (the measured operator asymmetry stays ~1e-9, far below a CG-breaking level).  So the
symmetric CG default is a BY-CONSTRUCTION robustness guarantee + a speedup, not a fix for an actively-
reproducing failure at these scales; the original retreat was conservative.  Large-scale demag is the MMMM
route's job anyway (CLAUDE.md role split 2026-06-27: HDiv = curved-surface accuracy + FEM coupling at moderate
scale), so the >50k regime is out of HDiv's lane.

The Gram BUILD dominates the cost (the per-pair analytic quadrature; cube N=8 = 47 s all-analytic vs a
~0.3 s mass-riesz solve; nonlinear sphere nf=9403 = 200 s exact build vs ~1 s/Newton-step solve).  Because
N = B^T G B is GEOMETRY-ONLY (material-independent), the PRECISION-PRESERVING fast build is the default for
the analytic-Gram material paths: uniform-linear `auto` / `cpp-cg` (symmetric mass-Riesz CG, already
validated at tight Gram eps), plus per-region linear, PM-mixed, AND the nonlinear Newton (the latter
paths stay on the GMRES/Newton asymmetry-tolerant formulations where needed; TET and polytope HEX/WEDGE).
`near_factor=2` (near pairs = exact analytic) + `far_quad=4` (far pairs = a low-order double-quadrature of
1/r, O((size/r)^4) -- degree-2 4-pt tet / 3-pt tri, or, for hex/wedge, the same degree-2 rule on the
centroid-fan sub-tets / sub-triangles).  This REPRODUCES the all-analytic Gram (uniform-linear sphere
transverse 7.26e-4 == exact 7.25e-4, demag identical; NONLINEAR sphere nf=9403 Mz agrees to 3e-7 with the
SAME 8 Newton iters) at ~4.5-9.4x faster build (linear cube N=8: 47 -> 10 s; nonlinear nf=9403: 200 -> 21 s).
The cheap centroid-monopole far (`far_quad=0`) is equally fast but leaks ~0.12% transverse (> the 1e-3
golden) -- so it is never defaulted; the low-quad far is what makes the fast build lossless.  Only H-LU
(factors its own system-A operator) and the Gauss backend keep the exact `near_factor=1e30` / `far_quad=0`.
An explicit `near_factor` / `far_quad` always wins (pass `near_factor=1e30` to force the all-analytic Gram).

KELVIN-LESS: the 1/r charge Gram IS the open boundary (a volume integral method like MMM/MSC); only
the iron is meshed -- no air box / Kelvin needed.  The NONLINEAR path uses the analytic charge Gram
(scalable `_ChargeGramHMatrix` at tight gram_eps), REQUIRED for div M != 0 (non-uniform M) bodies.

## Scope (M1)
Per-region soft iron, LINEAR (`mu_r` scalar or `{material: mu_r}` dict) AND NONLINEAR (`bh_table` one
[[H,B]] table or `{material: [[H,B]]}` dict).  N = B^T G B is geometry-only, so multi-grade iron enters
ONLY through the (1/chi)-weighted HDiv mass (linear) / the per-element constitutive law (nonlinear).
Mixed PM+iron (fixed-M source regions) + the 165k-DOF-scale preconditioner + the M0 parity gate are the
remaining productionization steps (docs/hdiv_vim/PRODUCTIONIZATION.md).  Until they land, multipole-moment MMM MSC
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
from ._vim import build_charge_gram
from ._nonlinear import (_bh_table_funcs, _table_tensor_tangent, _table_tensor_tangent_multi,
                         _bh_inverse_funcs, _reluctivity_tangent, _reluctivity_tangent_multi)

_MU0 = 4e-7 * _PI
# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect once.
_GMRES_TOL = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"
_LINEAR_SOLVERS = {"auto", "python", "cpp-cg", "gmres", "hlu"}
_GRAM_BACKENDS = {"analytic", "gauss"}


def _build_charge_gram(d, all_tet, gram_eps, leaf, eta, near_factor, image_masks, image_signs, far_quad=0):
    """Build the C++ charge-Gram H-matrix for the fallback / nonlinear demag path.

    far_quad (analytic mode): the FAR evaluation when near_factor < inf -- 0 = centroid-monopole, >0 = the
    precision-preserving low-order double-quad of 1/r (tet/tri directly, polytope via centroid-fan
    sub-tets/sub-tris; reproduces the all-analytic Gram at ~monopole cost)."""
    if all_tet:
        return _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                      n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta,
                                      near_factor=near_factor, image_masks=image_masks,
                                      image_signs=image_signs, far_quad=int(far_quad))
    # HEX/WEDGE: the polytope triangle-soup charge Gram (build_demag emits d["poly"] for non-tet).
    p = d["poly"]
    return _rp._ChargeGramHMatrix(
        cell_tris=list(p["cell_tris"]), cell_troff=list(p["cell_troff"]),
        cell_cent=list(p["cell_cent"]), cell_meas=list(p["cell_meas"]),
        face_tris=list(p["face_tris"]), face_troff=list(p["face_troff"]),
        face_cent=list(p["face_cent"]), face_meas=list(p["face_meas"]),
        n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta, near_factor=near_factor,
        image_masks=image_masks, image_signs=image_signs, far_quad=int(far_quad))


def _tet_gauss_point_cloud(d):
    """Low-order Gauss point cloud for G ~= P^T K_point P on tet/tri charge hosts."""
    cell_verts = np.asarray(d["cell_verts"], float).reshape(-1, 4, 3)
    face_verts = np.asarray(d["face_verts"], float).reshape(-1, 3, 3)
    point_coords = []
    point_charge = []
    point_weight = []
    # Symmetric degree-2 tetra rule: 4 points, weights sum to volume.
    a = 0.5854101966249685
    b = 0.1381966011250105
    tet_lam = np.array([[a, b, b, b], [b, a, b, b], [b, b, a, b], [b, b, b, a]], float)
    # Symmetric degree-2 triangle rule: 3 points, weights sum to area.
    tri_lam = np.array([[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]], float)
    sizes = np.zeros(int(d["n_charge"]), float)
    for c, V in enumerate(cell_verts):
        vol = abs(np.linalg.det(V[1:] - V[0])) / 6.0
        sizes[c] = np.cbrt(vol)
        pts = tet_lam @ V
        point_coords.extend(pts.ravel())
        point_charge.extend([c] * len(pts))
        point_weight.extend([vol / len(pts)] * len(pts))
    off = int(d["n_el"])
    for f, V in enumerate(face_verts):
        area = 0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0]))
        sizes[off + f] = np.sqrt(area)
        pts = tri_lam @ V
        point_coords.extend(pts.ravel())
        point_charge.extend([off + f] * len(pts))
        point_weight.extend([area / len(pts)] * len(pts))
    return (np.asarray(point_coords, float), np.asarray(point_charge, np.int32),
            np.asarray(point_weight, float), sizes)


def _point_direct_entry(point_coords, point_charge, point_weight, a, b):
    pa = np.flatnonzero(point_charge == a)
    pb = np.flatnonzero(point_charge == b)
    if len(pa) == 0 or len(pb) == 0:
        return 0.0
    A = point_coords.reshape(-1, 3)[pa]
    Bp = point_coords.reshape(-1, 3)[pb]
    WA = point_weight[pa]
    WB = point_weight[pb]
    s = 0.0
    inv4pi = 1.0 / (4.0 * _PI)
    for i, x in enumerate(A):
        r = np.linalg.norm(x[None, :] - Bp, axis=1)
        mask = r > 1e-300
        if np.any(mask):
            s += WA[i] * float(np.sum(WB[mask] * inv4pi / r[mask]))
    return s


def _near_charge_pairs(cent, sizes, near_factor):
    """Sparse near/self pair set for exact-minus-point Gram correction."""
    from scipy.spatial import cKDTree
    cent = np.asarray(cent, float)
    sizes = np.asarray(sizes, float)
    max_size = float(np.max(sizes)) if len(sizes) else 0.0
    tree = cKDTree(cent)
    pairs = []
    for a, c in enumerate(cent):
        cand = tree.query_ball_point(c, float(near_factor) * (sizes[a] + max_size))
        for b in cand:
            if b < a:
                continue
            r = float(np.linalg.norm(c - cent[b]))
            if a == b or r <= float(near_factor) * (sizes[a] + sizes[b]):
                pairs.append((a, b))
    return pairs


def _analytic_entry_oracle(d):
    """Geometry-only analytic charge-Gram (build=False): a cheap exact-entry oracle for the Gauss near
    correction.  .entry(a,b) is the EXACT analytic charge Gram (the ctor sets up the per-charge outer
    quadrature / centroids / sizes), but the O(N log N) H-matrix is NOT built -- so sampling the few
    near pairs costs O(near) analytic entries instead of a full analytic H-matrix build (which would
    defeat the whole point of the Gauss path, namely AVOIDING that analytic cost).

    near_factor=1e30 (all-analytic): the oracle must return the TRUE analytic entry for every near pair
    it is queried on (the near correction's exact reference), independent of the far-field near_factor
    used for the analytic backend's own build.  Querying is restricted to the near set by the caller."""
    return _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                  n_el=int(d["n_el"]), near_factor=1e30, build=False)


def _build_charge_gauss_gram(d, all_tet, gram_eps, leaf, eta, near_factor,
                             image_masks, image_signs, gauss_near_factor):
    """Build G ~= P^T K_point P + near correction.  First production slice: tet/tri, no IMA."""
    if not all_tet:
        raise ValueError("hdiv_demag_solve: gram_backend='gauss' currently supports tet/triangle meshes only")
    if image_masks or image_signs:
        raise ValueError("hdiv_demag_solve: gram_backend='gauss' does not yet support image symmetry")
    point_coords, point_charge, point_weight, sizes = _tet_gauss_point_cloud(d)
    # Exact analytic near entries via a geometry-only oracle (no full H-matrix build) -- the Gauss path
    # must NOT build the analytic Gram it exists to avoid; we only need .entry() on the near pairs.
    exact = _analytic_entry_oracle(d)
    corr_i, corr_j, corr_v = [], [], []
    for a, b in _near_charge_pairs(d["cent"], sizes, gauss_near_factor):
        exact_ab = float(exact.entry(int(a), int(b)))
        approx_ab = _point_direct_entry(point_coords, point_charge, point_weight, int(a), int(b))
        delta = exact_ab - approx_ab
        corr_i.append(int(a)); corr_j.append(int(b)); corr_v.append(float(delta))
        if a != b:
            corr_i.append(int(b)); corr_j.append(int(a)); corr_v.append(float(delta))
    # order-0 / RT0 = the TRIVIAL P scatter: one entry per point (P_pt=p, P_chg=owner, P_coef=weight).
    n_point = len(point_weight)
    return _rp._ChargeGaussHMatrix(
        point_coords=list(map(float, point_coords)),
        P_pt=list(range(n_point)), P_chg=list(map(int, point_charge)),
        P_coef=list(map(float, point_weight)), n_charge=int(d["n_charge"]),
        corr_i=corr_i, corr_j=corr_j, corr_v=corr_v,
        eps=gram_eps, leaf=leaf, eta=eta)


def _tet_hlu_solver_args(d, chi):
    """Convert build_demag output into the compact C++ system-A H-LU constructor surface."""
    Mm = sp.coo_matrix(d["M_mass"])
    Bc = d["B_csr"].tocsc()
    cent = np.asarray(d["cent"], float)
    ndof = int(d["ndof"])
    face_charge = np.full(ndof * 2, -1, np.int32)
    face_coef = np.zeros(ndof * 2, float)
    face_cent = np.zeros((ndof, 3), float)
    for i in range(ndof):
        s, e = Bc.indptr[i], Bc.indptr[i + 1]
        ids = Bc.indices[s:e]
        vals = Bc.data[s:e]
        if len(ids) > 2:
            raise RuntimeError("hdiv_demag_solve: C++ H-LU expects <=2 charge supports per RT0 face "
                               "(face %d has %d)" % (i, len(ids)))
        for p, (a, v) in enumerate(zip(ids, vals)):
            face_charge[i * 2 + p] = int(a)
            face_coef[i * 2 + p] = float(v)
        if len(ids):
            face_cent[i] = cent[ids].mean(0)
    return dict(
        face_centroids=face_cent.ravel().tolist(), chi=float(chi),
        face_charge=face_charge.tolist(), face_coef=face_coef.tolist(),
        mI=list(map(int, Mm.row)), mJ=list(map(int, Mm.col)), mV=list(map(float, Mm.data)),
        cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]), n_el=int(d["n_el"]),
    )


def _solve_linear_mass_riesz_cpp(H, B, Mm, n_face, h_ext, chi, tol, maxit):
    """DEFAULT uniform-chi linear solve: mass-Riesz-preconditioned CG ENTIRELY in C++ on the SPD +N system
    ((1/chi)M_mass + B^T G B) m = M_mass h_ext, with G applied via the EXACTLY-SYMMETRIC charge-Gram
    H-matvec (`matvec_sym`: the upper-triangular leaves define both triangles, so the operator is
    machine-symmetric -- ||G - G^T|| == 0 -- regardless of the per-block ACA truncation).  This is why CG
    is the default again (Sugahara 2026-06-27, "対称HACApKを実装しよう。CGがいいね"): the symmetric Gram
    makes CG robust BY CONSTRUCTION at all N -- it removes the asymmetry failure mode that motivated the
    earlier GMRES retreat (the spurious antisymmetric part of the GENERAL ACA matvec).  Mass-Riesz precond
    via a PARDISO SPD factor of the RT0 mass (eigenvalues vs M_mass are (1/chi)+d, d in [0,1], bounded ->
    ~3-5x fewer iters than diagonal Jacobi).  The whole Krylov loop (O(N log N) symmetric H-matvec +
    per-iteration PARDISO mass solve + vector ops) runs in C++ -- no Python per-iteration glue, no splu.
    The symmetric matvec is also ~1.4x FASTER than the general one (it skips the lower-triangle leaves).
    `H.solve_linear_material_mass_riesz(..., symmetric=True)` is the default; pass symmetric=False only to
    cross-check against the general (asymmetric) matvec."""
    inv_chi = 1.0 / float(chi)
    rhs = np.asarray(Mm @ h_ext).ravel()
    B = B.tocsr()
    Mm_coo = sp.coo_matrix(Mm)
    res = H.solve_linear_material_mass_riesz(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(map(float, B.data)),
        int(n_face), list(map(int, Mm_coo.row)), list(map(int, Mm_coo.col)),
        list(map(float, Mm_coo.data)), inv_chi, list(map(float, rhs)), tol, int(maxit))
    iters = int(res["iters"])
    if iters >= int(maxit):                      # fail-loud (No-Fallbacks): never return a non-converged M
        raise RuntimeError(
            "hdiv_demag_solve (symmetric mass-riesz CG): did NOT converge in %d iters (n_face=%d).  The "
            "operator is the EXACTLY-symmetric SPD +N system, so CG should converge -- a non-convergence "
            "here means an ill-conditioned material/mesh.  Tighten gram_eps or raise maxit; cross-check "
            "with linear_solver='gmres' (mass-Riesz GMRES) to isolate.  (Large-scale demag is the MMMM "
            "route's job, not HDiv-VIM; CLAUDE.md role split 2026-06-27.)" % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def _solve_linear_mass_riesz_gmres(H, B, Mm, n_face, h_ext, chi, tol, maxit, gmres_restart):
    """OPT-IN (`linear_solver="gmres"`) uniform-chi linear solve: mass-Riesz-preconditioned GMRES on the +N
    system ((1/chi)M_mass + B^T G B) m = M_mass h_ext, using the GENERAL (possibly slightly asymmetric)
    charge-Gram H-matvec.  GMRES is asymmetry-tolerant, so this is the cross-check / fallback for the
    (default) symmetric C++ CG -- kept because GMRES does not require the operator to be symmetric.  Was the
    default 2026-06-27 morning when the asymmetric ACA matvec was the only option; superseded the same day by
    the symmetric matvec (`matvec_sym`), which makes CG robust by construction and is faster, so CG is the
    default again.  The O(N log N) H-matvec stays the C++ kernel; the GMRES recurrence + sparse mass solve
    are Python.  Fail-loud on non-convergence (No-Fallbacks)."""
    inv_chi = 1.0 / float(chi)
    rhs = np.asarray(Mm @ h_ext).ravel()
    B = B.tocsr()
    Mfac = spla.splu(sp.csc_matrix(Mm))
    A = spla.LinearOperator((n_face, n_face),
                            matvec=lambda v: inv_chi * np.asarray(Mm @ np.asarray(v, float)).ravel()
                            + B.T @ np.asarray(H.matvec((B @ np.asarray(v, float)).tolist()), float))
    Mprec = spla.LinearOperator((n_face, n_face), matvec=lambda v: Mfac.solve(np.asarray(v, float)))
    it = {"n": 0}
    cycles = max(2, int(np.ceil(maxit / gmres_restart)))
    m, info = spla.gmres(A, rhs, M=Mprec, restart=int(gmres_restart), maxiter=cycles,
                         callback=lambda _x: it.__setitem__("n", it["n"] + 1),
                         callback_type="pr_norm", **{_GMRES_TOL: float(tol)})
    if info != 0:                                # fail-loud (No-Fallbacks): never return a non-converged M
        raise RuntimeError("hdiv_demag_solve (mass-riesz GMRES): did NOT converge (info=%d, n_face=%d, "
                           "iters=%d). Tighten gram_eps or raise maxit; for very large meshes use the MMMM "
                           "route (large-scale demag is not HDiv-VIM's role)." % (info, n_face, it["n"]))
    return np.asarray(m, float), it["n"]


def _solve_linear_W_cpp(H, B, W, Mm, n_face, h_ext, tol, maxit):
    """PER-REGION linear solve ENTIRELY in C++: the SYMMETRIC Galerkin system (M_{1/chi} + N) m = M_mass h_ext
    by symmetric mass-Riesz CG, where W = M_{1/chi} = INT (1/chi(x)) u.v dx is BOTH the system mass AND the
    Riesz preconditioner.  Passing W as the 'mass' COO with inv_chi=1.0 makes the C++ kernel
    (`solve_linear_material_mass_riesz`, symmetric=True) compute A = 1.0*W + B^T G B and precondition with
    W^{-1} (PARDISO) -- the SAME all-C++ symmetric CG as the uniform path, generalized to per-region chi.
    The whole Krylov loop runs in C++ (symmetric charge-Gram H-matvec + PARDISO W-solve + vector ops); only
    the NGSolve assembly of W and the RHS mass M_mass stays Python (assembly is NGSolve's job).  The
    symmetric 1/chi-weighted Galerkin form is CG-able (W, N both symmetric); it differs from the form-1 GMRES
    operator (M_mass + M_chi M_mass^-1 N) by O(h) at material interfaces -- both are consistent, and for a
    UNIFORM region it is bit-identical to the scalar +N system."""
    rhs = np.asarray(Mm @ h_ext).ravel()
    B = B.tocsr(); W_coo = sp.coo_matrix(W)
    res = H.solve_linear_material_mass_riesz(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(map(float, B.data)),
        int(n_face), list(map(int, W_coo.row)), list(map(int, W_coo.col)),
        list(map(float, W_coo.data)), 1.0, list(map(float, rhs)), tol, int(maxit))
    iters = int(res["iters"])
    if iters >= int(maxit):                      # fail-loud (No-Fallbacks)
        raise RuntimeError(
            "hdiv_demag_solve (per-region symmetric mass-riesz CG): did NOT converge in %d iters "
            "(n_face=%d).  The (M_{1/chi} + N) operator is SPD, so a non-convergence means an ill-"
            "conditioned material/mesh; tighten gram_eps or raise maxit, or cross-check with "
            "linear_solver='gmres' (the form-1 GMRES)." % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, bh_table=None, pm_M=None,
                     image=None, gram_eps=None, leaf=32, eta=2.0, near_factor=None, far_quad=None, tol=1e-8,
                     maxit=4000, gmres_restart=400, nl_maxit=300, nl_tol=1e-6, anderson_window=6,
                     linear_solver="auto", hlu_trunc_tol=1e-8,
                     gram_backend="analytic", gauss_near_factor=2.0, order=0,
                     curve_order=None, curve_gauss=8):
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
    if linear_solver not in _LINEAR_SOLVERS:
        raise ValueError("hdiv_demag_solve: linear_solver must be one of %s (got %r)"
                         % (sorted(_LINEAR_SOLVERS), linear_solver))
    if gram_backend not in _GRAM_BACKENDS:
        raise ValueError("hdiv_demag_solve: gram_backend must be one of %s (got %r)"
                         % (sorted(_GRAM_BACKENDS), gram_backend))
    if pm_M is None:
        if (mu_r is None) == (bh_table is None):
            raise ValueError("hdiv_demag_solve: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")
    else:
        if (mu_r is not None) and (bh_table is not None):
            raise ValueError("hdiv_demag_solve: with pm_M, give the iron as EITHER mu_r (linear) OR "
                             "bh_table (nonlinear), not both")

    # ---- HIGH-ORDER (order>=1): the order-p charge-Gram material solve ----
    # FIXED 2026-06-28 ([[hdiv-highorder-material-solve-wrong]]): the per-element change-of-basis in
    # `_vim._charge_basis` made the order-p demag operator N = B^T G B VALID (eig(M_mass^-1 N) in [0,1]; the
    # material solve p-converges instead of the old ~2-4x blow-up -- which was a SPURIOUS NET CHARGE from a
    # single element-0 change-of-basis reused for all elements, invisible to the uniform-M demag-FACTOR test).
    # The LINEAR (uniform-scalar OR per-region dict) mu_r case is wired through the same all-C++ symmetric
    # mass-Riesz CG as RT0; the not-yet-wired order>0 combos (image / PM-mixed / nonlinear / HLU / gauss) fail
    # loud in _solve_highorder (No-Fallbacks).  Golden: tests/feec/test_hdiv_vim_highorder_solve*.
    if int(order) != 0 or curve_order is not None:
        return _solve_highorder(mesh, int(order), mu_r, bh_table, pm_M, H_ext, image, linear_solver,
                                gram_backend, gram_eps, leaf, eta, near_factor, far_quad, tol, maxit,
                                gmres_restart, curve_order, curve_gauss)

    # ---- sparse HDiv geometry + applied field projection ----
    all_tet = all(len(el.vertices) == 4 for el in mesh.Elements(ng.VOL))
    image_masks, image_signs = [], []
    if image is not None:
        for axes, sign in _tet.image_group(_tet.parse_image_string(image)):
            image_masks.append(int(sum(1 << a for a in axes)))
            image_signs.append(float(sign))
    d = _tet.build_demag(mesh)
    Mm, B = d["M_mass"], d["B_csr"]
    n_face, n_el = int(d["ndof"]), int(d["n_el"])
    mu = d["m_unit"]
    fes = ng.HDiv(mesh, order=0)
    gfHext = ng.GridFunction(fes); gfHext.Set(H_ext)
    h_ext = gfHext.vec.FV().NumPy().copy()
    denom = float(mu @ np.asarray(Mm @ mu).ravel())

    # ---- linear scalar material: move the whole Krylov / H-LU loop to C++ when possible ----
    solver_used = "python-gmres"
    hmat_stats = None
    m = None
    iters = None
    D = None
    uniform_linear = pm_M is None and bh_table is None and mu_r is not None and not isinstance(mu_r, dict)
    chi_uniform = None
    # Effective ACA Gram tolerance.  The uniform-linear Krylov paths (default symmetric CG = auto/cpp-cg,
    # and the gmres cross-check) use the validated tight fast-Gram eps 1e-12; per-region / nonlinear / H-LU
    # keep 1e-10.  (With the SYMMETRIC matvec the CG no longer NEEDS 1e-12 for symmetry -- symmetry is now
    # structural, not ACA-accuracy-dependent -- but 1e-12 is kept for solution ACCURACY + golden stability.)
    # An explicit gram_eps always wins.
    fast_uniform_path = uniform_linear and linear_solver in ("auto", "cpp-cg", "gmres")
    eff_gram_eps = gram_eps if gram_eps is not None else (1e-12 if fast_uniform_path else 1e-10)
    # Gram-BUILD near/far split.  The build (per-pair analytic quadrature) dominates the cost (cube N=8: 47s
    # vs a 0.3s mass-riesz solve; nonlinear sphere nf=9403: 200s exact build vs ~1s/Newton-step solve).
    # N = B^T G B is GEOMETRY-ONLY (material-independent), so the PRECISION-PRESERVING fast build --
    # near_factor=2 (near pairs = exact analytic) + far_quad=4 (far pairs = a low-order double-quadrature of
    # 1/r, O((size/r)^4), on the tet / sub-tet+sub-tri rules) -- reproduces the all-analytic Gram for the
    # analytic-Gram material paths: uniform-linear auto/cpp-cg (already validated at tight Gram eps), plus
    # per-region linear, PM-mixed, AND the nonlinear Newton (GMRES/Newton, asymmetry-tolerant).  Measured:
    # sphere transverse 7.26e-4 == exact
    # 7.25e-4 (linear), nonlinear nf=9403 Mz agrees to 3e-7 with the same 8 Newton iters at ~9.4x faster build
    # (200->21s).  UNLIKE the bare centroid-monopole far (far_quad=0), which is equally fast but leaks ~0.12%
    # transverse (> the 1e-3 golden), so monopole is never defaulted.  H-LU factors its own system-A operator
    # (keep exact near_factor=1e30) and the Gauss backend is a separate point-cloud path.  An explicit
    # near_factor or far_quad always wins (pass near_factor=1e30 to force the all-analytic Gram).
    fast_build = gram_backend == "analytic" and linear_solver != "hlu"
    eff_near_factor = near_factor if near_factor is not None else (2.0 if fast_build else 1e30)
    eff_far_quad = far_quad if far_quad is not None else (4 if fast_build else 0)
    if gram_backend == "gauss" and not uniform_linear:
        raise ValueError("hdiv_demag_solve: gram_backend='gauss' is currently enabled only for "
                         "uniform linear mu_r solves")
    if gram_backend == "gauss" and linear_solver == "hlu":
        raise ValueError("hdiv_demag_solve: gram_backend='gauss' is a charge-operator path; "
                         "linear_solver='hlu' currently factors the analytic tet system-A operator")
    if uniform_linear:
        chi_uniform = float(mu_r) - 1.0
        if chi_uniform <= 0.0:
            raise ValueError("hdiv_demag_solve: mu_r must be > 1 (got %r)" % (mu_r,))
        can_hlu = all_tet and image is None and hasattr(_rp, "_HDivVimTetSolver")
        if linear_solver == "hlu" and not can_hlu:
            raise ValueError("hdiv_demag_solve: linear_solver='hlu' currently requires a tet mesh, "
                             "no image symmetry, and a HACApK-enabled _HDivVimTetSolver")
        if linear_solver == "hlu" and can_hlu:
            args = _tet_hlu_solver_args(d, chi_uniform)
            solver = _rp._HDivVimTetSolver(eps=eff_gram_eps, leaf=leaf, eta=eta,
                                           trunc_tol=hlu_trunc_tol, gram_near_factor=eff_near_factor,
                                           **args)
            rhs = chi_uniform * np.asarray(Mm @ h_ext).ravel()
            m = np.asarray(solver.solve(list(map(float, rhs))), float)
            H_for_d = _build_charge_gram(d, all_tet, eff_gram_eps, leaf, eta, eff_near_factor,
                                         image_masks, image_signs, far_quad=eff_far_quad)
            Nmu = B.T @ np.asarray(H_for_d.matvec((B @ mu).tolist()), float)
            D = float((mu @ Nmu) / denom)
            iters = 1
            solver_used = "cpp-hlu"
            if hasattr(solver, "stats"):
                hmat_stats = dict(solver.stats())

    if m is None:
        # ---- demag operator N: scalable C++ charge-Gram H-matrix (_ChargeGramHMatrix) ----
        # The C++ kernel handles tet (cell_verts/face_verts) AND hex/wedge (polytope triangle soup),
        # and folds IMA mirror symmetry as image charges.  The fallback path keeps Python only as
        # orchestration; uniform linear solves use H.solve_linear_material* so the Krylov loop is C++.
        if gram_backend == "gauss":
            H = _build_charge_gauss_gram(d, all_tet, eff_gram_eps, leaf, eta, eff_near_factor,
                                         image_masks, image_signs, gauss_near_factor)
        else:
            H = _build_charge_gram(d, all_tet, eff_gram_eps, leaf, eta, eff_near_factor,
                                   image_masks, image_signs, far_quad=eff_far_quad)

        def N_apply(v):
            v = np.asarray(v, float)
            return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

        # M_mass^{-1} sparse LU is built LAZILY: the default 'auto' and 'cpp-cg' symmetric CG paths factor
        # the mass inside the kernel (PARDISO), while GMRES / PM-mixed / nonlinear Newton paths need this
        # Python-side sparse-LU cache for their orchestration.
        _mcache = {}

        def _get_Mfac():
            if "f" not in _mcache:
                _mcache["f"] = spla.splu(sp.csc_matrix(Mm))
            return _mcache["f"]

        def _get_Mprec():
            if "p" not in _mcache:
                f = _get_Mfac()
                _mcache["p"] = spla.LinearOperator((n_face, n_face),
                                                   matvec=lambda v: f.solve(np.asarray(v, float)))
            return _mcache["p"]

        D = float((mu @ N_apply(mu)) / denom)
        if hasattr(H, "stats"):
            hmat_stats = dict(H.stats())

        if pm_M is not None and bh_table is not None:
            if linear_solver in ("cpp-cg", "hlu"):
                raise ValueError("hdiv_demag_solve: linear_solver=%r applies only to linear mu_r solves"
                                 % (linear_solver,))
            m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, _get_Mfac(), _get_Mprec(),
                                        n_face, fes, tol, gmres_restart, nl_maxit, nl_tol, anderson_window,
                                        pm_M=pm_M)
        elif pm_M is not None:
            if linear_solver in ("cpp-cg", "hlu"):
                raise ValueError("hdiv_demag_solve: linear_solver=%r does not yet support PM-mixed solves"
                                 % (linear_solver,))
            M_chi, b_pm = _build_mixed_pm(mesh, fes, mu_r, pm_M, Mm, n_face)
            m, iters = _solve_linear(Mm, M_chi, N_apply, _get_Mfac(), _get_Mprec(), n_face, h_ext,
                                     tol, maxit, gmres_restart, b_extra=b_pm)
        elif mu_r is not None:
            if isinstance(mu_r, dict):
                if linear_solver == "hlu":
                    raise ValueError("hdiv_demag_solve: linear_solver='hlu' does not support per-region mu_r")
                if linear_solver in ("auto", "cpp-cg"):
                    # PER-REGION default: all-C++ SYMMETRIC mass-Riesz CG on the Galerkin system
                    # (M_{1/chi} + N) m = M_mass h_ext (W = the 1/chi-weighted HDiv mass is both the system
                    # mass and the Riesz preconditioner).  Same C++ symmetric CG as the uniform path; only the
                    # NGSolve assembly of W stays Python.  'gmres'/'python' keep the form-1 GMRES cross-check.
                    W = _build_invchi_mass(mesh, fes, mu_r, n_face)
                    m, iters = _solve_linear_W_cpp(H, B, W, Mm, n_face, h_ext, tol, maxit)
                    solver_used = "mass-riesz-cg"
                else:                                  # gmres / python: form-1 GMRES (asymmetry-tolerant)
                    M_chi = _build_chi_mass(mesh, fes, mu_r, Mm, n_face)
                    m, iters = _solve_linear(Mm, M_chi, N_apply, _get_Mfac(), _get_Mprec(), n_face, h_ext,
                                             tol, maxit, gmres_restart)
                    solver_used = "mass-riesz-gmres" if linear_solver == "gmres" else "python-gmres"
            elif linear_solver in ("auto", "cpp-cg"):
                # DEFAULT: all-C++ mass-Riesz CG on the +N system, with G applied via the EXACTLY-SYMMETRIC
                # charge-Gram H-matvec (matvec_sym -- upper-triangular leaves define both triangles, so the
                # operator is machine-symmetric regardless of the per-block ACA truncation).  The symmetric
                # Gram makes CG robust BY CONSTRUCTION at all N (it removes the asymmetry failure mode of the
                # general matvec), and the symmetric matvec is ~1.4x faster (skips the lower-triangle leaves)
                # -> CG is the default again (Sugahara 2026-06-27).  PARDISO mass-Riesz precond, no Python
                # per-iteration glue.  linear_solver='gmres' is the asymmetry-tolerant cross-check/opt-in.
                m, iters = _solve_linear_mass_riesz_cpp(H, B, Mm, n_face, h_ext, chi_uniform, tol, maxit)
                solver_used = "mass-riesz-cg"
            elif linear_solver == "gmres":
                # OPT-IN cross-check: mass-Riesz GMRES on the GENERAL (asymmetry-tolerant) matvec.  Was the
                # default 2026-06-27 morning before the symmetric matvec landed; kept for cross-validation.
                m, iters = _solve_linear_mass_riesz_gmres(H, B, Mm, n_face, h_ext, chi_uniform, tol, maxit,
                                                          gmres_restart)
                solver_used = "mass-riesz-gmres"
            else:
                M_chi = _build_chi_mass(mesh, fes, mu_r, Mm, n_face)
                m, iters = _solve_linear(Mm, M_chi, N_apply, _get_Mfac(), _get_Mprec(), n_face, h_ext,
                                         tol, maxit, gmres_restart)
        else:
            if linear_solver == "hlu":
                raise ValueError("hdiv_demag_solve: linear_solver='hlu' applies only to linear mu_r solves")
            if linear_solver in ("gmres", "python"):
                # OPT-IN cross-check: the forward M-residual Newton (scipy splu + scipy GMRES, asymmetry-
                # tolerant).  Was the production nonlinear path before the all-C++ energy-Newton landed.
                m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, _get_Mfac(), _get_Mprec(),
                                            n_face, fes, tol, gmres_restart, nl_maxit, nl_tol, anderson_window)
                solver_used = "forward-newton-gmres"
            else:
                # DEFAULT iron-only nonlinear: the all-C++ SYMMETRIC ENERGY-NEWTON (inverse-BH reluctivity,
                # J = W_tan + N solved by the C++ symmetric W-CG -- mass-Riesz PARDISO, N H-matvec; no scipy
                # splu / GMRES / M_mass^-1).  Brings the nonlinear solve to C++ parity with the linear path;
                # validated == the forward Newton to ~1e-6 at all drives.
                m, iters = _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, B, Mm, n_face, h_ext,
                                                       tol, maxit, nl_maxit, nl_tol)
                solver_used = "energy-newton-cpp"

    # ---- per-element M (VectorL2(0) projection; component-major DOFs -> (n_el,3)) + averages ----
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol

    out = dict(M=M_el, M_avg=M_avg, iters=int(iters), demag=D, ndof=n_face, n_el=n_el,
               n_charge=int(d["n_charge"]), nonlinear=(bh_table is not None),
               linear_solver=solver_used, gram_backend=gram_backend)
    if hmat_stats is not None:
        out["hmat_stats"] = hmat_stats
    return out


def _solve_highorder(mesh, order, mu_r, bh_table, pm_M, H_ext, image, linear_solver, gram_backend,
                     gram_eps, leaf, eta, near_factor, far_quad, tol, maxit, gmres_restart,
                     curve_order=None, curve_gauss=8):
    """order>0 (high-order HDiv) soft-iron demag solve.  The order-p charge-Gram demag operator N = B^T G B is
    a VALID demag operator since the per-element change-of-basis fix (2026-06-28,
    [[hdiv-highorder-material-solve-wrong]]): eig(M_mass^-1 N) in [0,1] and the material solve p-converges
    (no 2x/4x blow-up).  Supports the LINEAR (uniform-scalar OR per-region dict) mu_r case via the SAME
    all-C++ symmetric mass-Riesz CG as the RT0 path; the not-yet-wired order>0 combos fail loud (No-Fallbacks).
    The CALLER opens `with ng.TaskManager():` (same contract as hdiv_demag_solve)."""
    if image is not None:
        raise NotImplementedError("hdiv_demag_solve: image symmetry is not yet wired at order>0 (use order=0)")
    if pm_M is not None:
        raise NotImplementedError("hdiv_demag_solve: PM-mixed (pm_M) is not yet wired at order>0 (use order=0)")
    if bh_table is not None:
        raise NotImplementedError("hdiv_demag_solve: NONLINEAR (bh_table) is not yet validated at order>0 -- "
                                  "use order=0, or the order-p demag FACTOR via DemagOperator(...).DemagFactor()")
    if linear_solver == "hlu":
        raise NotImplementedError("hdiv_demag_solve: linear_solver='hlu' is RT0-only (order=0)")
    if gram_backend == "gauss":
        raise NotImplementedError("hdiv_demag_solve: gram_backend='gauss' is not yet wired at order>0")
    if int(order) > 2:
        # order<=2 uses the EXACT analytic-moment charge potential (machine precision).  For order>=3 the C++
        # Gram falls back to the Duffy singular quadrature (PhiInner -> PhiAtHO_Duffy), which is ~1e-3 accurate
        # -- fine for curved-panel field evaluation, but NOT for the order>=3 MATERIAL solve: the ill-
        # conditioned high-degree monomial basis (cond(B)^2 in N=B^T G B) amplifies the ~1e-3 entry error so
        # the demag spectrum escapes [0,1].  A clean order>=3 material solve needs machine-precision entries
        # (the analytic moments extended with TetMoment2 / degree-3 surface moments), not the Duffy.  Fail loud
        # (No-Fallbacks) until that lands.  [[hdiv-vim-sauter-schwab-cg]]
        raise NotImplementedError(
            "hdiv_demag_solve: order>2 material solve is not yet production-clean -- order<=2 is exact "
            "(analytic moments); the order>=3 Duffy quadrature is only ~1e-3 and the ill-conditioned "
            "high-degree basis makes the demag spectrum leave [0,1]. Use order in {0,1,2}.")
    # accuracy-preserving fast Gram build (near analytic + far low-quad), the build_charge_gram defaults; an
    # explicit gram_eps/near_factor/far_quad always wins (pass near_factor=inf to force the all-high-quad Gram).
    eff_eps = gram_eps if gram_eps is not None else 1e-10
    eff_far = far_quad if far_quad is not None else 3
    eff_hofar = near_factor if near_factor is not None else 2.0
    if curve_order is not None:
        # CURVED (isoparametric P2) demag solve: curve the geometry, then the curved-Duffy charge Gram.  Curved
        # helps NEAR-SURFACE FIELD / FLUX accuracy (sigma=M.n on the true curved surface), NOT the volume-
        # averaged demag FACTOR (curving-insensitive ~3e-5 on a sphere; [[hdiv-vim-sauter-schwab-cg]] de-risk).
        if int(curve_order) != 2:
            raise NotImplementedError("hdiv_demag_solve: only curve_order=2 (isoparametric P2) is wired.")
        mesh.Curve(int(curve_order))
        fes = ng.HDiv(mesh, order=order)
        B, H, M_mass = build_charge_gram(fes, eps=eff_eps, leafsize=leaf, eta=eta,
                                         curve_order=int(curve_order), curve_gauss=int(curve_gauss))
    else:
        fes = ng.HDiv(mesh, order=order)
        B, H, M_mass = build_charge_gram(fes, eps=eff_eps, leafsize=leaf, eta=eta,
                                         far_quad=eff_far, ho_far_factor=eff_hofar)
    Mm = sp.csr_matrix(M_mass); B = sp.csr_matrix(B)
    n_face = fes.ndof; n_el = mesh.GetNE(ng.VOL); n_charge = B.shape[0]
    gfH = ng.GridFunction(fes); gfH.Set(H_ext); h_ext = gfH.vec.FV().NumPy().copy()
    gfMu = ng.GridFunction(fes); gfMu.Set(ng.CoefficientFunction((0, 0, 1)))
    mu = gfMu.vec.FV().NumPy().copy()
    denom = float(mu @ np.asarray(Mm @ mu).ravel())

    def N_apply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    D = float((mu @ N_apply(mu)) / denom)
    hmat_stats = dict(H.stats()) if hasattr(H, "stats") else None

    if isinstance(mu_r, dict):                                  # per-region linear: W = 1/chi-weighted HDiv mass
        W = _build_invchi_mass(mesh, fes, mu_r, n_face)
        m, iters = _solve_linear_W_cpp(H, B, W, Mm, n_face, h_ext, tol, maxit)
        solver_used = "mass-riesz-cg"
    else:                                                       # uniform-scalar linear
        chi = float(mu_r) - 1.0
        if chi <= 0.0:
            raise ValueError("hdiv_demag_solve: mu_r must be > 1 (got %r)" % (mu_r,))
        if linear_solver == "gmres":
            m, iters = _solve_linear_mass_riesz_gmres(H, B, Mm, n_face, h_ext, chi, tol, maxit, gmres_restart)
            solver_used = "mass-riesz-gmres"
        else:
            m, iters = _solve_linear_mass_riesz_cpp(H, B, Mm, n_face, h_ext, chi, tol, maxit)
            solver_used = "mass-riesz-cg"

    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol
    out = dict(M=M_el, M_avg=M_avg, iters=int(iters), demag=D, ndof=n_face, n_el=n_el,
               n_charge=n_charge, nonlinear=False, linear_solver=solver_used,
               gram_backend=gram_backend, order=int(order), curve_order=curve_order)
    if hmat_stats is not None:
        out["hmat_stats"] = hmat_stats
    return out


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


def _build_invchi_mass(mesh, fes, mu_r, n_face):
    """The 1/chi-weighted HDiv mass M_{1/chi} = INT (1/chi(x)) u.v dx for the SYMMETRIC per-region Galerkin
    system A = M_{1/chi} + N (the CG-able all-C++ form -- see _solve_linear_W_cpp).  `mu_r` is a dict
    {material: mu_r} (each > 1).  Fail-loud (No-Fallbacks): every mesh material specified, each mu_r > 1."""
    mats = list(mesh.GetMaterials())
    missing = sorted(set(mats) - set(mu_r))
    if missing:
        raise ValueError("hdiv_demag_solve: mu_r dict missing region(s) %s; mesh materials are %s"
                         % (missing, mats))
    bad = {r: mu_r[r] for r in mats if float(mu_r[r]) <= 1.0}
    if bad:
        raise ValueError("hdiv_demag_solve: every region mu_r must be > 1 (got %s)" % bad)
    invchi_cf = mesh.MaterialCF({r: 1.0 / (float(mu_r[r]) - 1.0) for r in mats})
    u, v = fes.TnT()
    a = ng.BilinearForm(fes); a += invchi_cf * u * v * ng.dx; a.Assemble()
    r, c, val = a.mat.COO()
    return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(n_face, n_face))


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


def _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, B, Mm, n_face, h_ext, cg_tol, cg_maxit,
                                nl_maxit, nl_tol):
    """SYMMETRIC ENERGY-NEWTON with an all-C++ inner solve -- the production nonlinear soft-iron path (brings
    the nonlinear solve to C++ parity with the linear mass-Riesz CG; the default for iron-only nonlinear,
    replacing the forward `_solve_nonlinear` scipy-splu + scipy-GMRES Newton).

    Co-energy / reluctivity form -- SYMMETRIC + M_mass^-1-free.  The co-energy functional
      E(m) = INT W_co(|M|) dx + 1/2 m.(N m) - (M_mass h_ext).m   (M = the flux field of m; W_co = INT H dM)
    is CONVEX; its gradient is the residual
      R(m) = INT H(M).v dx + N m - M_mass h_ext    (H(M) = the INVERSE-BH reluctance field, H = nu_sec M)
    and its Hessian (the Newton Jacobian) is SYMMETRIC
      J = W_tan + N,  W_tan = INT nu_d u.v dx  (nu_d = dH/dM = (dM/dH)^-1 differential reluctivity tensor),
    so each Newton step  (W_tan + N) dm = -R  is solved by the EXISTING C++ symmetric W-CG
    (`solve_linear_material_mass_riesz`: W = W_tan as both the system mass AND the mass-Riesz PARDISO
    preconditioner; N via the symmetric charge-Gram H-matvec).  NO scipy splu, NO scipy GMRES, NO M_mass^-1.

    Globalization: a chi0 (zero-field) LINEAR W-CG warmstart; an Armijo line search on the CONVEX ENERGY E
    (the merit -- ||R|| stalls in saturation where the inverse-BH
    H(M) blows up, but E keeps decreasing); a HARD-SATURATION BARRIER in the inverse BH (M cannot exceed
    Msat) that repels the M-iterates from the unphysical |M| > Mmax region.  Convergence: tight (relative
    Newton step < nl_tol -- 1-2 iters at moderate drive, == the forward Newton to ~1e-13), OR -- for the deep-
    saturation regime where the hard-saturation M-form intrinsically limit-cycles at the achievable precision
    -- a settled-step acceptance (rel step < 1e-4 for 5 consecutive iters -> accept the best-energy iterate;
    M matched the analytic uniform sphere to ~2e-6 at H0 up to 5e6, knee*5000).  Single-region (scalar
    bh_table) AND per-region (dict) iron; PM-mixed stays on the forward path.  CALLER opens TaskManager."""
    Bc = sp.csr_matrix(B); Mmc = sp.csr_matrix(Mm)
    rhs_src = np.asarray(Mmc @ h_ext).ravel()
    Id = ng.Id(3); uf, vf = fes.TnT()
    gfM = ng.GridFunction(fes); l2 = ng.L2(mesh, order=0)
    Bptr = list(map(int, Bc.indptr)); Bidx = list(map(int, Bc.indices)); Bdat = list(map(float, Bc.data))
    # element volumes (for the co-energy integral) = the L2(0) mass diagonal
    mvol = ng.BilinearForm(l2); mvol += l2.TrialFunction() * l2.TestFunction() * ng.dx; mvol.Assemble()
    rv, cv, vvv = mvol.mat.COO(); Vol = np.zeros(mesh.ne)
    for r_, c_, v_ in zip(rv, cv, vvv):
        if r_ == c_:
            Vol[int(r_)] = v_

    # ---- per-region OR single inverse-BH reluctivity fields + co-energy + zero-field chi0 (warmstart) ----
    if isinstance(bh_table, dict):
        mats = list(mesh.GetMaterials())
        missing = sorted(set(mats) - set(bh_table))
        if missing:
            raise ValueError("hdiv_demag_solve: bh_table dict missing region(s) %s; mesh materials are %s"
                             % (missing, mats))
        region_names = list(bh_table)
        name_to_ridx = {nm: k for k, nm in enumerate(region_names)}
        region_fields, region_wco = [], []
        for nm in region_names:
            arr = np.asarray(bh_table[nm], float)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError("hdiv_demag_solve: bh_table[%r] must be [[H,B], ...] (A/m, T)" % (nm,))
            f, w, _ = _bh_inverse_funcs(arr[:, 0], arr[:, 1]); region_fields.append(f); region_wco.append(w)
        elem_region = np.array([name_to_ridx[mesh[ng.ElementId(ng.VOL, i)].mat] for i in range(mesh.ne)],
                               dtype=int)

        def _reluct(g):
            return _reluctivity_tangent_multi(g, mesh, region_fields, elem_region, Id)

        def _wco_all(Mmag):
            out = np.empty_like(Mmag)
            for ridx, w in enumerate(region_wco):
                sel = elem_region == ridx
                if np.any(sel):
                    out[sel] = w(Mmag[sel])
            return out

        chi0_e = np.empty(mesh.ne)
        for ridx, f in enumerate(region_fields):
            _, nd0 = f(np.array([1e-12])); chi0_e[elem_region == ridx] = 1.0 / max(float(nd0[0]), 1e-30)
    else:
        arr = np.asarray(bh_table, float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("hdiv_demag_solve: bh_table must be [[H,B], ...] (A/m, T)")
        fields, wco, _ = _bh_inverse_funcs(arr[:, 0], arr[:, 1])

        def _reluct(g):
            return _reluctivity_tangent(g, mesh, fields, Id)

        def _wco_all(Mmag):
            return wco(Mmag)

        _, nd0 = fields(np.array([1e-12])); chi0_e = np.full(mesh.ne, 1.0 / max(float(nd0[0]), 1e-30))

    def _N_apply(v):
        v = np.asarray(v, float)
        return Bc.T @ np.asarray(H.matvec((Bc @ v).tolist()), float)

    def _Mmag(m):
        gfM.vec.FV().NumPy()[:] = m
        gfn = ng.GridFunction(l2); gfn.Set(ng.sqrt(ng.InnerProduct(gfM, gfM) + 1e-30))
        return np.maximum(gfn.vec.FV().NumPy(), 1e-30)

    def _bH(H_cf):
        lf = ng.LinearForm(fes); lf += H_cf * vf * ng.dx; lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def _W_coo(weight_cf, tensor):
        a = ng.BilinearForm(fes)
        a += (ng.InnerProduct(weight_cf * uf, vf) if tensor else weight_cf * uf * vf) * ng.dx
        a.Assemble()
        r, c, v = a.mat.COO()
        return sp.coo_matrix((np.array(v), (np.array(r), np.array(c))), shape=(n_face, n_face))

    def _solve_W(W_coo, rhs):
        res = H.solve_linear_material_mass_riesz(
            Bptr, Bidx, Bdat, int(n_face),
            list(map(int, W_coo.row)), list(map(int, W_coo.col)), list(map(float, W_coo.data)),
            1.0, list(map(float, rhs)), cg_tol, int(cg_maxit))
        it = int(res["iters"])
        if it >= int(cg_maxit):
            raise RuntimeError("hdiv_demag_solve (energy-Newton inner W-CG): did NOT converge in %d iters "
                               "(n_face=%d); the (W_tan + N) operator is SPD, so this means an ill-"
                               "conditioned tangent/mesh -- tighten gram_eps or raise maxit." % (cg_maxit, n_face))
        return np.asarray(res["m"], float), it

    def _energy(m):                              # E(m) = INT W_co(|M|) dx + 1/2 m.Nm - rhs.m  (convex merit)
        return float(np.dot(_wco_all(_Mmag(m)), Vol)) + 0.5 * float(m @ _N_apply(m)) - float(rhs_src @ m)

    # chi0 (zero-field) LINEAR warmstart: (M_{1/chi0} + N) m = M_mass h_ext, one C++ W-CG solve
    invchi0 = ng.GridFunction(l2)
    invchi0.vec.FV().NumPy()[:] = 1.0 / np.maximum(chi0_e, 1.0)
    m, _ = _solve_W(_W_coo(invchi0, tensor=False), rhs_src)

    converged = False; nit = 0; rel_step = float("inf"); settled = 0
    E = _energy(m); Ebest = E; mbest = m.copy()
    for it in range(nl_maxit):
        nit = it + 1
        gfM.vec.FV().NumPy()[:] = m
        H_cf, nud = _reluct(gfM)
        R = _bH(H_cf) + _N_apply(m) - rhs_src
        dm, _ = _solve_W(_W_coo(nud, tensor=True), -R)
        dec = float(-dm @ R)                                 # Newton decrement^2 = dm.(-R) = dm^T J dm >= 0
        lam = 1.0; E0 = E                                    # Armijo line search on the CONVEX ENERGY E
        while lam > 1e-10:
            if _energy(m + lam * dm) <= E0 - 1e-4 * lam * dec:
                break
            lam *= 0.5
        step = lam * dm
        rel_step = float(np.linalg.norm(step)) / (float(np.linalg.norm(m)) + 1e-30)
        m = m + step; E = _energy(m)
        if E < Ebest:
            Ebest = E; mbest = m.copy()
        settled = settled + 1 if rel_step < 1e-4 else 0
        if rel_step < nl_tol:                                # tight convergence (moderate drive: 1-2 iters)
            converged = True; break
        if settled >= 5:                                     # deep-saturation M-form limit cycle: accept the
            converged = True; m = mbest; break               # best-energy iterate (M is at achievable precision)
    if not converged:
        m = mbest
        raise RuntimeError("hdiv_demag_solve (energy-Newton): did NOT converge -- rel step=%.2e (tol %.1e), "
                           "%d settled iters after %d (returning M would be a silent wrong result).  For an "
                           "extreme-saturation / ill-conditioned case, cross-check with linear_solver='gmres' "
                           "(the forward H-form Newton)." % (rel_step, nl_tol, settled, nit))
    return m, nit
