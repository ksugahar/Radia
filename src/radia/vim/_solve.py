"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`vim.Solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r), and
`vim.Solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the FEEC/HDiv counterpart to the multipole-moment MMM MSC hex/wedge/pyramid soft-iron demag in `rad.Solve`.
Both modes take an ARBITRARY applied field `H_ext` (any NGSolve CoefficientFunction -- e.g. a coil's
Biot-Savart field `rad.RadiaField(coil,'h')`, the C-type electromagnet driver) and return per-element M.

## Formulation (verified-first, 2026-06-15)
ONE projected weak form everywhere -- the magnetization M is the RT1 primary, the constitutive law
M = M(H) is imposed in the L2 sense (M_mass m = INT M(H).v dx), and H = h_ext - M_mass^-1 N m is the
weak total field (N = B^T G B, h_ext = H_ext L2-projected onto HDiv order 1).  LINEAR soft iron is the
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
`((1/chi)M_mass + N) m = M_mass h_ext` by CG, preconditioned with the FULL RT1 H(div) mass inverse
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
`linear_solver="python"` is the form-1 GMRES + `M_mass^{-1}` sparse LU, used for per-region chi /
nonlinear Newton paths until their material-specific operators move to C++.  The old system-A H-LU path
is retired from the public HDiv-VIM backend.  (The mass Riesz makes the
operator well-conditioned by construction -- the earlier "h-explosion => need AMS" was a monopole-Gram
artifact; the accurate analytic Gram + mass Riesz needs no auxiliary-space preconditioner.)

The uniform-linear Krylov paths (default CG = auto/cpp-cg, and the gmres cross-check) build the analytic Gram
at the tight `gram_eps=1e-12`; per-region / nonlinear keep `1e-10`.  (With the symmetric matvec the CG
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
For RT0, an explicit `near_factor` / `far_quad` always wins (pass `near_factor=1e30` to force the
all-analytic Gram).  For high-order, use `ho_far_factor` for the separation threshold.

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
_LINEAR_SOLVERS = {"auto", "python", "cpp-cg", "gmres"}   # 'hlu' retired (RT0-only system-A H-LU)
_GRAM_BACKENDS = {"analytic"}                             # 'gauss' retired (RT0 build-speed experiment)


def _resolve_gram_params(*, order, gram_backend, linear_solver, uniform_linear, gram_eps,
                         near_factor, far_quad, ho_far_factor):
    """Resolve the charge-Gram BUILD defaults (ACA eps + near/far split) in ONE place.

    Single source of truth -- was inline + duplicated in vim.Solve (order=0) and _solve_highorder
    (order>0).  `gram_eps` / `far_quad` apply to both paths; `near_factor` is RT0-only and
    `ho_far_factor` is high-order-only.  Wrong-order knobs fail loud instead of being silently remapped.

    RT0 (order=0):
      eps      -- the uniform-linear symmetric-CG + gmres cross-check paths use the validated tight 1e-12
                  (kept for solution ACCURACY + golden stability even though the SYMMETRIC matvec no longer
                  NEEDS it for symmetry); per-region / nonlinear / H-LU keep 1e-10.
      near/far -- the Gram BUILD (per-pair analytic quadrature) dominates the cost (cube N=8: 47s build vs a
                  0.3s mass-riesz solve; nonlinear sphere nf=9403: 200s exact build vs ~1s/Newton-step solve).
                  N=B^T G B is GEOMETRY-ONLY, so the precision-preserving fast build near_factor=2 (near =
                  exact analytic) + far_quad=4 (far = a low-order double-quad of 1/r, O((size/r)^4))
                  REPRODUCES the all-analytic Gram (sphere transverse 7.26e-4 == exact 7.25e-4; nonlinear Mz
                  3e-7, same Newton iters) at ~monopole cost.  The bare centroid-monopole far (far_quad=0) is
                  equally fast but leaks ~0.12% transverse, so it is NEVER defaulted.  H-LU factors its own
                  system-A operator -> keep exact (near=1e30, far=0).  Gauss backend is a separate path.

    High-order (order>0): far_quad=3 (the low far-quad) + ho_far_factor=2.0 (the separation threshold).
      near_factor (RT0) and ho_far_factor (high-order) are ORDER-SPECIFIC and NOT interchangeable -- a
      wrong-order knob fails loud (No-Fallbacks), it is not silently remapped.
    """
    if int(order) == 0:
        if ho_far_factor is not None:
            raise ValueError("vim.Solve: ho_far_factor is an order>0 (high-order) parameter; "
                             "at order=0 (RT0) use near_factor for the near/far split")
        fast_uniform = uniform_linear and linear_solver in ("auto", "cpp-cg", "gmres")
        fast_build = gram_backend == "analytic" and linear_solver != "hlu"
        return {
            "eps": gram_eps if gram_eps is not None else (1e-12 if fast_uniform else 1e-10),
            "near_factor": near_factor if near_factor is not None else (2.0 if fast_build else 1e30),
            "far_quad": far_quad if far_quad is not None else (4 if fast_build else 0),
        }
    if near_factor is not None:
        raise ValueError("vim.Solve: near_factor is an order=0 (RT0) parameter; at order>0 "
                         "(high-order) use ho_far_factor for the near/far separation threshold")
    return {
        "eps": gram_eps if gram_eps is not None else 1e-10,
        "far_quad": far_quad if far_quad is not None else 3,
        "ho_far_factor": ho_far_factor if ho_far_factor is not None else 2.0,
    }


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
            "vim.Solve (symmetric mass-riesz CG): did NOT converge in %d iters (n_face=%d).  The "
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
        raise RuntimeError("vim.Solve (mass-riesz GMRES): did NOT converge (info=%d, n_face=%d, "
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
            "vim.Solve (per-region symmetric mass-riesz CG): did NOT converge in %d iters "
            "(n_face=%d).  The (M_{1/chi} + N) operator is SPD, so a non-convergence means an ill-"
            "conditioned material/mesh; tighten gram_eps or raise maxit, or cross-check with "
            "linear_solver='gmres' (the form-1 GMRES)." % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, bh_table=None, pm_M=None, magnets=None,
                     image=None, gram_eps=None, leaf=32, eta=2.0, near_factor=None, far_quad=None, tol=1e-8,
                     maxit=4000, gmres_restart=400, nl_maxit=300, nl_tol=1e-6, anderson_window=6,
                     linear_solver="auto", hlu_trunc_tol=1e-8,
                     gram_backend="analytic", gauss_near_factor=2.0, order=1,
                     curve_order=None, curve_gauss=8, ho_far_factor=None):
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
    near/far Gram-build tuning (ORDER-SPECIFIC, NOT interchangeable -- a wrong-order knob fails loud):
      near_factor   -- order=0 (RT0) ONLY: the near(=exact analytic)/far(=far_quad) distance boundary
                       (pass 1e30 to force the all-analytic Gram).
      ho_far_factor -- order>0 (high-order) ONLY: the near/far separation threshold (pass inf to force the
                       all-high-quad Gram).
    Returns dict: M (n_el,3) per-element magnetization, M_avg (3,), iters, demag (Rayleigh factor),
    ndof, n_el, n_charge, nonlinear(bool).  The caller must open `with ng.TaskManager():`.
    """
    if H_ext is None:
        raise ValueError("vim.Solve: H_ext (applied-field CoefficientFunction) is required")
    # ---- HDiv-VIM scope: RT1 (order 1) on a pure-TET, pure-HEX, or pure-WEDGE mesh ----
    # RT0 retired: per-element INACCURATE (the demag FACTOR is right ~1/3, but the per-element M leaks --
    # raising the solution order to RT1 is what fixes it); RT2+ retired: no per-element gain over RT1, slower.
    # TET, pure-HEX (Q1 volume charge + Q2 geometry), and pure-WEDGE solve LINEAR + NONLINEAR through
    # _solve_highorder here -- the C++ energy-Newton is Gram-agnostic (verified 2026-07-04: hex cube linear
    # & nonlinear match collocation MMMM to ~1%).  Pyramid / mixed + pm_M remain the collocation MMMM
    # backend's job (rad.Solve); curved IMA, 'gauss', and 'hlu' fail loud.
    if int(order) != 1:
        raise ValueError(
            "HDiv-VIM is RT1 (HDiv order=1) only.  RT0 (order=0) is RETIRED -- it is per-element INACCURATE "
            "(the demag factor is right ~1/3 but the per-element magnetization leaks; RT1 is what fixes it); "
            "use the collocation MMMM backend for a low-order surface-charge demag.  RT2+ is RETIRED too (no "
            "per-element gain over RT1, markedly slower).  Pass order=1 (the default).  (The geometry "
            "curve_order is a SEPARATE knob: curve_order=2 isoparametric P2 is still allowed.)")
    if pm_M is not None:
        raise NotImplementedError(
            "HDiv-VIM (RT1) is a SOFT-IRON (linear / nonlinear) demag solver and does not mix permanent "
            "magnets (pm_M).  Place permanent magnets as direct-M elements solved by collocation MMMM "
            "(rad.Solve), or as an applied-field source folded into H_ext.")
    # IMA mirror symmetry (image=) is wired for flat pure-TET / pure-HEX / pure-WEDGE RT1 (2026-07-04);
    # reduced curved / mixed / pyramid cases fail loud downstream instead of silently dropping the image.
    if mesh.dim == 2:
        if image is not None:
            raise NotImplementedError(
                "vim.Solve (2D): the planar HDiv-VIM layer does not support IMA image symmetry; "
                "use a full (un-reduced) 2D model.")
        # ---- PLANAR (2D motor cross-section) branch: the dense planar layer (_vim2d) ----
        # The 2D layer supports the core single-region surface: mu_r (linear) / bh_table
        # (nonlinear) + H_ext.  The 3D-only knobs must stay at their defaults -- fail loud.
        if gram_backend != "analytic":
            raise ValueError("vim.Solve (2D): gram_backend must be 'analytic' -- the planar "
                             "layer has one Gram (the C++ 2D log-kernel charge Gram)")
        if linear_solver != "auto":
            raise ValueError("vim.Solve (2D): linear_solver must be 'auto' (the planar layer "
                             "is dense; got %r)" % (linear_solver,))
        for _nm, _val in (("gram_eps", gram_eps), ("near_factor", near_factor),
                          ("far_quad", far_quad), ("ho_far_factor", ho_far_factor),
                          ("curve_order", curve_order)):
            if _val is not None:
                raise ValueError("vim.Solve (2D): %s is a 3D knob; the 2D Gram parameters "
                                 "are fixed by its own gates (got %r)" % (_nm, _val))
        from ._vim2d import solve_planar_demag
        return solve_planar_demag(mesh, mu_r=mu_r, H_ext=H_ext, bh_table=bh_table, magnets=magnets,
                                  eta=eta, nl_tol=nl_tol, nl_maxit=nl_maxit)
    if magnets is not None:
        raise NotImplementedError(
            "vim.Solve: magnets= (separate-body permanent-magnet source) is wired for the 2D "
            "planar layer only; in 3D place PMs as direct-M collocation MMMM elements or fold their "
            "field into H_ext.")
    _vtx = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if _vtx not in ({4}, {8}, {6}):
        raise ValueError(
            "HDiv-VIM (order 1) supports a pure-TET (4-vertex), pure-HEX (8-vertex), or pure-WEDGE/prism "
            "(6-vertex) mesh; got vertex counts %s.  pyramid / MIXED-element soft-iron demag uses the "
            "collocation MMMM backend (rad.Solve demag_backend='collocation_mmmm'), which rad.Solve's "
            "'auto' split uses for unsupported mesh-backed iron." % sorted(_vtx))
    if linear_solver not in _LINEAR_SOLVERS:
        raise ValueError("vim.Solve: linear_solver must be one of %s (got %r)"
                         % (sorted(_LINEAR_SOLVERS), linear_solver))
    if gram_backend not in _GRAM_BACKENDS:
        raise ValueError("vim.Solve: gram_backend must be one of %s (got %r)"
                         % (sorted(_GRAM_BACKENDS), gram_backend))
    if pm_M is None:
        if (mu_r is None) == (bh_table is None):
            raise ValueError("vim.Solve: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")
    else:
        if (mu_r is not None) and (bh_table is not None):
            raise ValueError("vim.Solve: with pm_M, give the iron as EITHER mu_r (linear) OR "
                             "bh_table (nonlinear), not both")

    # AUTO-MATCH: a CURVED mesh (mesh.GetCurveOrder()>=2) needs a Gram built on the SAME curved geometry as
    # B/M_mass, else N=B^T G B (straight Gram) is geometry-inconsistent and the demag DRIFTS with geometry order
    # (tet sphere: straight-Gram 0.336/0.308/0.279 at curve 1/2/3; matched curved Gram restores ~1/3 -- 0.338 at
    # curve 2).  Only curve_order=2 (isoparametric P2) is wired in build_charge_gram.  curve_order=0 forces the
    # straight Gram (a deliberate flat-Gram probe); an explicit int overrides the auto-match.  (TET enforced above.)
    if curve_order is None:
        _k = mesh.GetCurveOrder()
        if _k >= 2:
            curve_order = _k
    elif curve_order == 0:
        curve_order = None

    # ---- RT1 (HDiv order=1): the order-1 charge-Gram material solve (the SOLE HDiv-VIM solve path) ----
    # The per-element change-of-basis in `_vim._charge_basis` (2026-06-28, [[hdiv-highorder-material-solve-wrong]])
    # makes the order-1 demag operator N = B^T G B valid (eig(M_mass^-1 N) in [0,1]; per-element M p-converges).
    # LINEAR (uniform-scalar OR per-region dict) mu_r AND flat/curved NONLINEAR (bh_table) are wired via the same
    # all-C++ symmetric mass-Riesz CG / energy-Newton.  order is always 1 here (the tet/RT1 guard above), so this
    # is unconditional.  Golden: validation_test/feec/test_hdiv_vim_demag_solve*, test_hdiv_vim_curved_solve_nonlinear*.
    return _solve_highorder(mesh, int(order), mu_r, bh_table, pm_M, H_ext, image, linear_solver,
                            gram_backend, gram_eps, leaf, eta, near_factor, far_quad, tol, maxit,
                            gmres_restart, curve_order, curve_gauss, ho_far_factor, nl_maxit, nl_tol)


def _solve_highorder(mesh, order, mu_r, bh_table, pm_M, H_ext, image, linear_solver, gram_backend,
                     gram_eps, leaf, eta, near_factor, far_quad, tol, maxit, gmres_restart,
                     curve_order=None, curve_gauss=8, ho_far_factor=None, nl_maxit=300, nl_tol=1e-6):
    """order>0 (high-order HDiv) soft-iron demag solve.  The order-p charge-Gram demag operator N = B^T G B is
    a VALID demag operator since the per-element change-of-basis fix (2026-06-28,
    [[hdiv-highorder-material-solve-wrong]]): eig(M_mass^-1 N) in [0,1] and the material solve p-converges
    (no 2x/4x blow-up).  Supports the LINEAR (uniform-scalar OR per-region dict) mu_r case via the SAME
    all-C++ symmetric mass-Riesz CG as the RT0 path; the not-yet-wired order>0 combos fail loud (No-Fallbacks).
    The CALLER opens `with ng.TaskManager():` (same contract as vim.Solve)."""
    # IMA mirror symmetry: WIRED for the FLAT pure-TET (C++ highorder QuadDotRefl->PhiInner) AND pure-HEX /
    # pure-WEDGE (the C++ QuadBlockHex/Wedge(mask) reflected block) paths -- the Gram folds the mirror-image
    # charge interactions so a reduced 1/2,1/4,1/8 model reproduces the full model.  CURVED (curve_order) +
    # MIXED / pyramid fail loud (not yet wired) -> collocation MMMM.
    image_masks, image_signs = [], []
    if image is not None:
        if curve_order is not None:
            raise NotImplementedError(
                "vim.Solve: IMA image symmetry is not yet wired for the CURVED (curve_order) "
                "HDiv-VIM path -- use a flat tet/hex/wedge mesh, or collocation MMMM (rad.Solve).")
        _ivtx = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
        if _ivtx not in ({4}, {8}, {6}):
            raise NotImplementedError(
                "vim.Solve: IMA image symmetry is wired for the FLAT pure-TET / pure-HEX / pure-WEDGE "
                "RT1 Gram; MIXED / pyramid reduced models use collocation MMMM (rad.Solve "
                "demag_backend='collocation_mmmm', image=...).  Got vertex counts %s." % sorted(_ivtx))
        _planes = _tet.parse_image_string(image)
        # (2026-07-05) hex/wedge IMA now handles ANTISYMMETRIC (negative-sign, field-PERPENDICULAR) planes too:
        # the on-plane cut-face self-term is computed with the EXACT self-radial in the reflected block (the
        # QuadBlockHex/Wedge "R(host)==host -> self_pair" fix), so the large perpendicular cut-face charge
        # cancels exactly for sign -1 instead of the earlier ~1.5% hex / ~29% wedge quadrature residual.
        for axes, sign in _tet.image_group(_planes):
            image_masks.append(int(sum(1 << a for a in axes)))
            image_signs.append(float(sign))
    if pm_M is not None:
        raise NotImplementedError("vim.Solve: PM-mixed (pm_M) is not yet wired at order>0 (use order=0)")
    # flat (non-curved) order>0 NONLINEAR is now WIRED: the symmetric energy-Newton (_solve_nonlinear_energy_cpp)
    # is Gram-AGNOSTIC (it consumes only H.matvec + H.solve_linear_material_mass_riesz), so it runs on the flat
    # high-order Gram exactly as on the RT0 / curved Gram.  VERIFIED: flat RT1 nonlinear on a tet sphere matches
    # the RT0 nonlinear solve to ~7e-4 (golden: test_hdiv_vim_curved_solve_nonlinear::test_flat_rt1_nonlinear_*).
    if linear_solver == "hlu":
        raise NotImplementedError("vim.Solve: linear_solver='hlu' is RT0-only (order=0)")
    if gram_backend == "gauss":
        raise NotImplementedError("vim.Solve: gram_backend='gauss' is not yet wired at order>0")
    if int(order) > 2:
        # order<=2 uses the EXACT analytic-moment charge potential (machine precision).  For order>=3 the C++
        # Gram falls back to the Duffy singular quadrature (PhiInner -> PhiAtHO_Duffy), which is ~1e-3 accurate
        # -- fine for curved-panel field evaluation, but NOT for the order>=3 MATERIAL solve: the ill-
        # conditioned high-degree monomial basis (cond(B)^2 in N=B^T G B) amplifies the ~1e-3 entry error so
        # the demag spectrum escapes [0,1].  A clean order>=3 material solve needs machine-precision entries
        # (the analytic moments extended with TetMoment2 / degree-3 surface moments), not the Duffy.  Fail loud
        # (No-Fallbacks) until that lands.  [[hdiv-vim-sauter-schwab-cg]]
        raise NotImplementedError(
            "vim.Solve: order>2 material solve is not yet production-clean -- order<=2 is exact "
            "(analytic moments); the order>=3 Duffy quadrature is only ~1e-3 and the ill-conditioned "
            "high-degree basis makes the demag spectrum leave [0,1]. Use order in {0,1,2}.")
    # Gram-build defaults resolved in ONE place (_resolve_gram_params; rationale in its docstring).
    # _solve_highorder serves order>0 AND the order=0 curved (curve_order) path; BOTH resolve via the
    # HIGH-ORDER branch (which returns ho_far_factor).  max(order,1) routes order=0+curve_order there too
    # (eps 1e-10, far_quad 3, ho_far_factor 2.0 -- matching the pre-consolidation inline defaults; a plain
    # order=0 would wrongly hit the RT0 branch -> KeyError 'ho_far_factor' + the RT0 eps 1e-12).
    _gp = _resolve_gram_params(order=max(int(order), 1), gram_backend=gram_backend,
                               linear_solver=linear_solver, uniform_linear=False, gram_eps=gram_eps,
                               near_factor=near_factor, far_quad=far_quad, ho_far_factor=ho_far_factor)
    eff_eps = _gp["eps"]; eff_far = _gp["far_quad"]; eff_hofar = _gp["ho_far_factor"]
    if curve_order is not None:
        # CURVED (isoparametric P2) demag solve: curve the geometry, then the curved-Duffy charge Gram.  Curved
        # helps NEAR-SURFACE FIELD / FLUX accuracy (sigma=M.n on the true curved surface), NOT the volume-
        # averaged demag FACTOR (curving-insensitive ~3e-5 on a sphere; [[hdiv-vim-sauter-schwab-cg]] de-risk).
        if int(curve_order) != 2:
            raise NotImplementedError("vim.Solve: only curve_order=2 (isoparametric P2) is wired.")
        mesh.Curve(int(curve_order))
        fes = ng.HDiv(mesh, order=order)
        B, H, M_mass = build_charge_gram(fes, eps=eff_eps, leafsize=leaf, eta=eta,
                                         curve_order=int(curve_order), curve_gauss=int(curve_gauss),
                                         nonlinear=bh_table is not None)
    else:
        fes = ng.HDiv(mesh, order=order)
        B, H, M_mass = build_charge_gram(fes, eps=eff_eps, leafsize=leaf, eta=eta,
                                         far_quad=eff_far, ho_far_factor=eff_hofar,
                                         nonlinear=bh_table is not None,
                                         image_masks=image_masks, image_signs=image_signs)
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

    if bh_table is not None:                                    # CURVED nonlinear: energy-Newton on the curved Gram
        m, iters = _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, B, Mm, n_face, h_ext,
                                               tol, maxit, nl_maxit, nl_tol)
        solver_used = "energy-newton-cpp"
    elif isinstance(mu_r, dict):                                # per-region linear: W = 1/chi-weighted HDiv mass
        W = _build_invchi_mass(mesh, fes, mu_r, n_face)
        m, iters = _solve_linear_W_cpp(H, B, W, Mm, n_face, h_ext, tol, maxit)
        solver_used = "mass-riesz-cg"
    else:                                                       # uniform-scalar linear
        chi = float(mu_r) - 1.0
        if chi <= 0.0:
            raise ValueError("vim.Solve: mu_r must be > 1 (got %r)" % (mu_r,))
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
               n_charge=n_charge, nonlinear=bh_table is not None, linear_solver=solver_used,
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
            raise ValueError("vim.Solve: mu_r dict missing region(s) %s; mesh materials are %s"
                             % (missing, mats))
        bad = {r: mu_r[r] for r in mats if float(mu_r[r]) <= 1.0}
        if bad:
            raise ValueError("vim.Solve: every region mu_r must be > 1 (got %s)" % bad)
        chi_cf = mesh.MaterialCF({r: float(mu_r[r]) - 1.0 for r in mats})
        u, v = fes.TnT()
        a = ng.BilinearForm(fes); a += chi_cf * u * v * ng.dx; a.Assemble()
        r, c, val = a.mat.COO()
        return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(n_face, n_face))
    chi = float(mu_r) - 1.0
    if chi <= 0.0:
        raise ValueError("vim.Solve: mu_r must be > 1 (got %r)" % (mu_r,))
    return chi * Mm


def _build_invchi_mass(mesh, fes, mu_r, n_face):
    """The 1/chi-weighted HDiv mass M_{1/chi} = INT (1/chi(x)) u.v dx for the SYMMETRIC per-region Galerkin
    system A = M_{1/chi} + N (the CG-able all-C++ form -- see _solve_linear_W_cpp).  `mu_r` is a dict
    {material: mu_r} (each > 1).  Fail-loud (No-Fallbacks): every mesh material specified, each mu_r > 1."""
    mats = list(mesh.GetMaterials())
    missing = sorted(set(mats) - set(mu_r))
    if missing:
        raise ValueError("vim.Solve: mu_r dict missing region(s) %s; mesh materials are %s"
                         % (missing, mats))
    bad = {r: mu_r[r] for r in mats if float(mu_r[r]) <= 1.0}
    if bad:
        raise ValueError("vim.Solve: every region mu_r must be > 1 (got %s)" % bad)
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
        raise ValueError("vim.Solve: pm_M region(s) %s not in mesh materials %s" % (missing_pm, mats))
    if mu_r is None:
        iron_mu = {}
    elif isinstance(mu_r, dict):
        iron_mu = {r: float(mu_r[r]) for r in mu_r}
    else:
        iron_mu = {r: float(mu_r) for r in mats if r not in pm_regions}
    iron_regions = set(iron_mu)
    overlap = sorted(iron_regions & pm_regions)
    if overlap:
        raise ValueError("vim.Solve: region(s) %s are both iron (mu_r) and PM (pm_M)" % overlap)
    unspec = sorted(set(mats) - iron_regions - pm_regions)
    if unspec:
        raise ValueError("vim.Solve: region(s) %s are neither iron (mu_r) nor PM (pm_M); "
                         "mesh materials are %s" % (unspec, mats))
    bad = {r: iron_mu[r] for r in iron_regions if iron_mu[r] <= 1.0}
    if bad:
        raise ValueError("vim.Solve: every iron mu_r must be > 1 (got %s)" % bad)
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
            "vim.Solve: a permanent-magnet region directly TOUCHES a soft-iron region "
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
            raise ValueError("vim.Solve: pm_M region(s) %s not in mesh materials %s" % (missing_pm, mats))

    # ---- iron regions + per-region/single BH constitutive ----
    if isinstance(bh_table, dict):
        iron_regions = set(bh_table)
        if not pm_M:
            missing = sorted(set(mats) - iron_regions)
            if missing:
                raise ValueError("vim.Solve: bh_table dict missing region(s) %s; mesh materials are %s"
                                 % (missing, mats))
        region_names = list(bh_table)
        name_to_ridx = {nm: k for k, nm in enumerate(region_names)}
        region_funcs = []
        chi0_by_name = {}
        for nm in region_names:
            arr = np.asarray(bh_table[nm], float)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError("vim.Solve: bh_table[%r] must be [[H,B], ...] (A/m, T)" % (nm,))
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
            raise ValueError("vim.Solve: bh_table must be [[H,B], ...] (A/m, T)")
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
        raise ValueError("vim.Solve: region(s) %s are both iron (bh_table) and PM (pm_M)" % overlap)
    unspec = sorted(set(mats) - iron_regions - pm_regions)
    if unspec:
        raise ValueError("vim.Solve: region(s) %s are neither iron (bh_table) nor PM (pm_M); "
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
        raise RuntimeError("vim.Solve (linear): form-1 GMRES did not converge (info=%d, "
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
        raise RuntimeError("vim.Solve (nonlinear Newton): did NOT converge -- rel=%.2e > "
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
    -- a settled-step acceptance (rel step < 3e-4 for 5 consecutive iters -> accept the best-energy iterate;
    the RT1 limit-cycle plateau is ~1.5-1.9e-4, higher than RT0's <1e-4, so the floor is 3e-4; M matched the
    analytic uniform sphere to ~2e-3 at H0 up to 5e6, knee*5000).  Single-region (scalar
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
            raise ValueError("vim.Solve: bh_table dict missing region(s) %s; mesh materials are %s"
                             % (missing, mats))
        region_names = list(bh_table)
        name_to_ridx = {nm: k for k, nm in enumerate(region_names)}
        region_fields, region_wco = [], []
        for nm in region_names:
            arr = np.asarray(bh_table[nm], float)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError("vim.Solve: bh_table[%r] must be [[H,B], ...] (A/m, T)" % (nm,))
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
            raise ValueError("vim.Solve: bh_table must be [[H,B], ...] (A/m, T)")
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
            raise RuntimeError("vim.Solve (energy-Newton inner W-CG): did NOT converge in %d iters "
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
        # RT1's M-form limit cycle plateaus HIGHER than RT0's (~1.5-1.9e-4 rel step on a real BH table vs
        # <1e-4 at RT0 -- the larger RT1 charge system), so the settled-acceptance floor is 3e-4 (above the
        # observed RT1 plateau, below an actively-converging step).  The accepted M is the BEST-ENERGY iterate
        # (the energy minimum), VERIFIED on the sphere to ~2e-3 vs the analytic fixed point at H0 up to 5e6.
        settled = settled + 1 if rel_step < 3e-4 else 0
        if rel_step < nl_tol:                                # tight convergence (moderate drive: 1-2 iters)
            converged = True; break
        if settled >= 5:                                     # deep-saturation M-form limit cycle: accept the
            converged = True; m = mbest; break               # best-energy iterate (M is at achievable precision)
    if not converged:
        m = mbest
        raise RuntimeError("vim.Solve (energy-Newton): did NOT converge -- rel step=%.2e (tol %.1e), "
                           "%d settled iters after %d (returning M would be a silent wrong result).  For an "
                           "extreme-saturation / ill-conditioned case, cross-check with linear_solver='gmres' "
                           "(the forward H-form Newton)." % (rel_step, nl_tol, settled, nit))
    return m, nit
