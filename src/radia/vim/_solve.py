"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`vim.Solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r),
`vim.Solve(mesh, mu_r=.., B_r=.., H_ext=..)`   -- LINEAR recoil permanent magnet, and
`vim.Solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the production FEEC/HDiv soft-iron demag path used by `rad.Solve`.
Both modes take an ARBITRARY applied field `H_ext` (any NGSolve CoefficientFunction -- e.g. a coil's
Biot-Savart field `rad.RadiaField(coil,'h')`, the C-type electromagnet driver) and return per-element M.

## Formulation (verified-first, 2026-06-15)
ONE projected weak form everywhere -- the magnetization M is the BDM1/BDM2 primary, the constitutive law
M = M(H) is imposed in the L2 sense (M_mass m = INT M(H).v dx), and H = h_ext - M_mass^-1 N m is the
weak total field (N = B^T G B, h_ext = H_ext L2-projected onto the selected HDiv order).  LINEAR soft iron is the
CONSTANT-chi special case M = chi H, giving the form-1 system

    (M_mass + M_chi M_mass^-1 N) m = M_chi h_ext ,   M_chi = INT chi(x) u.v dx   (the chi-weighted mass),

solved +N (the -N system is mu_r-independent but NON-physical -- wrong sign).  For a SINGLE region
(uniform chi) this is identical bit-for-bit to ((1/chi) M_mass + N) m = M_mass h_ext (verified vs the
analytic sphere to 1.4e-4); for PER-REGION chi it is the consistent generalization, and it makes a
linear region agree with its nonlinear-table equivalent (the 1/chi-weighted form differs at material
interfaces by O(h)).  NONLINEAR iron solves the same projected residual by Newton (below).

NONLINEAR -- a **symmetric energy-Newton** method on the co-energy residual.  Its Hessian is
`W_tan + N`, with `W_tan` the differential-reluctivity HDiv mass and `N` the symmetric demag operator.
Each Newton step is one all-C++ preconditioned CG solve plus an Armijo line search, after a scalar-chi
linear warmstart.  WHY Newton, not the earlier Anderson-Hantila fixed point: Hantila/Picard's contraction
rho = (chi_max - chi_min)/(chi_max + chi_min) -> 1 as the iron saturates (chi_max ~ 1e4 unsaturated,
chi_min ~ 10 at the pole), so it STALLS on real silicon steel (measured: ~1e-2 residual after 300
iters, NOT converged) -- Newton's quadratic step is immune to the chi-range.  VERIFIED: the real CEFC
Si-steel C-type at 3000 AT converges relF 0.55 -> 8e-7 in ~24 iters, gap B within ~1% of the FEM truth
(Anderson-Hantila could not solve it at all).  N = B^T G B stays SPD-PSD (G the SPD Coulomb Gram), the
de-Rham loops sit in ker(N) and are carried by M_mass (no loop-star), and the scalable C++ charge-Gram
H-matvec is the only O(N log N) cost.

## Linear solve dispatch -- SYMMETRIC C++ CG (symmetric HACApK)
For the common scalar-mu_r case, `linear_solver="auto"` (the DEFAULT) solves the SPD +N system
`((1/chi)M_mass + N) m = M_mass h_ext` by CG, preconditioned with the FULL BDM1 H(div) mass inverse
`M_mass^{-1}` (the MASS RIESZ map), ENTIRELY in C++ (PARDISO mass factor + C++ Krylov, no Python glue).
CG is used because the charge-Gram is applied via the EXACTLY-SYMMETRIC H-matvec (`matvec_sym`): the
HACApK H-matrix stores both (I,J) and (J,I) leaf blocks but ACA-truncates them INDEPENDENTLY, so the
GENERAL matvec is only approximately symmetric; `matvec_sym` instead applies the UPPER-triangular leaves
only -- each upper leaf supplies its own block AND the mirror as its exact transpose -- so the operator is
machine-symmetric (||G - G^T|| == 0) regardless of the per-block truncation.  This makes CG robust BY
CONSTRUCTION at all N (it removes the asymmetry failure mode entirely), and the symmetric matvec is ~1.4x
FASTER than the general one (it touches half the leaves).  So CG is the default again (Sugahara 2026-06-27,
"対称HACApKを実装しよう。CGがいいね").  `linear_solver="cpp-cg"` is an explicit name for this symmetric C++ CG.
The old Python nonsymmetric Krylov, system-A H-LU, and Python sparse-LU paths are removed from HDiv-VIM.
(The mass Riesz makes the
operator well-conditioned by construction -- the earlier "h-explosion => need AMS" was a monopole-Gram
artifact; the accurate analytic Gram + mass Riesz needs no auxiliary-space preconditioner.)

The uniform-linear C++ CG path (auto/cpp-cg) builds the analytic Gram
at the tight `gram_eps=1e-12`; per-region / nonlinear keep `1e-10`.  (With the symmetric matvec the CG
no longer NEEDS 1e-12 for symmetry -- symmetry is now STRUCTURAL, independent of the ACA accuracy -- 1e-12 is
kept only for solution ACCURACY + golden stability.)  An explicit `gram_eps` always wins.  All material solve
paths are fail-loud: a non-converged solve RAISES (No-Fallbacks) rather than returning a wrong M.

ON THE EARLIER "+N CG SCALE WALL" (recorded honestly): an earlier nonsymmetric-solver retreat was motivated by a
report that the GENERAL-matvec CG diverges past nf ~ 20k (the spurious ACA antisymmetry growing with N).  When
the symmetric matvec was added, a re-measurement could NOT reproduce that divergence at HDiv scales: even with
the lossy monopole far (max asymmetry) + a distorted hex + mu_r up to 1e6, the GENERAL-matvec CG converges
fine through nf ~ 51k (the measured operator asymmetry stays ~1e-9, far below a CG-breaking level).  So the
symmetric CG default is a BY-CONSTRUCTION robustness guarantee + a speedup, not a fix for an actively
reproducing failure at these scales; the original retreat was conservative.

The Gram BUILD dominates the cost (the per-pair analytic quadrature; cube N=8 = 47 s all-analytic vs a
~0.3 s mass-riesz solve; nonlinear sphere nf=9403 = 200 s exact build vs ~1 s/Newton-step solve).  Because
N = B^T G B is GEOMETRY-ONLY (material-independent), the PRECISION-PRESERVING fast build is the default for
the analytic-Gram material paths: uniform-linear `auto` / `cpp-cg` (symmetric mass-Riesz CG, already
validated at tight Gram eps), plus per-region linear and nonlinear Newton.  Prescribed fixed-M sources
enter only through the assembled applied-field LinearForm; all material paths use the same symmetric
C++ operator on TET and polytope HEX/WEDGE.
`ho_far_factor=2` (near pairs = exact analytic) + `far_quad=4` (far pairs = a low-order double-quadrature of
 1/r, O((size/r)^4) -- degree-2 4-pt tet / 3-pt tri, or, for hex/wedge, the same degree-2 rule on the
centroid-fan sub-tets / sub-triangles).  This REPRODUCES the all-analytic Gram (uniform-linear sphere
transverse 7.26e-4 == exact 7.25e-4, demag identical; NONLINEAR sphere nf=9403 Mz agrees to 3e-7 with the
SAME 8 Newton iters) at ~4.5-9.4x faster build (linear cube N=8: 47 -> 10 s; nonlinear nf=9403: 200 -> 21 s).
The cheap centroid-monopole far (`far_quad=0`) is equally fast but leaks ~0.12% transverse (> the 1e-3
golden) -- so it is never defaulted; the low-quad far is what makes the fast build lossless.  Use
`ho_far_factor` for the separation threshold; `inf` forces the all-high-order reference build.

KELVIN-LESS: the 1/r charge Gram IS the open boundary; only
the iron is meshed -- no air box / Kelvin needed.  The NONLINEAR path uses the analytic charge Gram
(scalable `_ChargeGramHMatrix` at tight gram_eps), REQUIRED for div M != 0 (non-uniform M) bodies.

## Scope (M1)
Per-region soft iron, LINEAR (`mu_r` scalar or `{material: mu_r}` dict) AND NONLINEAR (`bh_table` one
[[H,B]] table or `{material: [[H,B]]}` dict).  N = B^T G B is geometry-only, so multi-grade iron enters
ONLY through the (1/chi)-weighted HDiv mass (linear) / the per-element constitutive law (nonlinear).
Permanent magnets use the documented four-level ladder: fixed/given MagnetizationSource, this solve's
linear-recoil B_r law, simplified Play, and full B-input EnergyStop.  Fixed PM + iron source coupling is
live.  `vim.SolveCoupled` supplies mutually coupled linear-recoil PM + nonlinear iron and segmented
linear-recoil bodies on independent spaces.  `vim.SolveCoupledHysteresis` supplies the stateful
EnergyStop/Play PM + nonlinear-iron history coupling.  This entry is also the soft-iron demag backend
used by `rad.Solve`.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- this library helper does NOT
open a TaskManager; the caller wraps the call in `with ng.TaskManager():`.
"""
import math
import os
import time
from math import pi as _PI

import numpy as np
import scipy.sparse as sp
import ngsolve as ng

import radia._radia_pybind as _rp
from . import _image
from ._vim import (
    _curve_mesh,
    build_charge_gram,
    _volume_vertex_counts,
)
from ._capabilities import validate_hdiv_configuration
from ._nonlinear import (_bh_table_funcs, _table_tensor_tangent, _table_tensor_tangent_multi,
                         _bh_inverse_funcs, _reluctivity_tangent, _reluctivity_tangent_multi,
                         _validate_bh_table)


class _OperatorBackedResult(dict):
    """Mapping result with a private, non-mapping geometry-cache attachment."""

    __slots__ = ("_operator_cache",)

    def __init__(self, *args, operator_cache=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._operator_cache = operator_cache


_MU0 = 4e-7 * _PI
_LINEAR_SOLVERS = {"auto", "cpp-cg"}
_NONLINEAR_SOLVERS = {"energy-newton", "picard-mass-riesz", "picard-energy"}
_PRECONDITIONERS = {"auto", "mass-riesz", "jacobi"}
_LAST_CPP_SOLVE_TIMINGS = {}
_LAST_NONLINEAR_SOLVE_STATS = {}
_AUTO_TET_JACOBI_NFACE_DEFAULT = 6000
_AUTO_TET_JACOBI_NFACE_ENV = "RADIA_HDIV_AUTO_JACOBI_TET_NFACE"


def _hdiv_space_with_image_constraints(
        mesh, order, image_planes, cyclic_periodic_boundaries=()):
    """Create the NGSolve HDiv space, removing symmetry-odd normal traces.

    An image plane with positive scalar-potential parity has odd normal
    magnetization.  Its normal trace is therefore exactly zero on the cut
    plane.  Leaving those trace unknowns in the reduced problem makes the
    result depend on tiny quadrature cross-couplings; an explicit full mesh
    has no corresponding internal-face charge.  NGSolve's compressed-space
    wrapper is the natural FEEC representation of this essential trace
    condition and keeps every assembled matrix on the native C++ path.
    """
    base = ng.HDiv(mesh, order=int(order))
    periodic_slave_dofs = ()
    if cyclic_periodic_boundaries:
        if image_planes:
            raise NotImplementedError(
                "vim.Solve: mirror image constraints and an azimuthal "
                "periodic HDiv trace cannot be combined")
        if not mesh.ngmesh.GetIdentifications():
            raise ValueError(
                "vim.Solve: cyclic_periodic_boundaries requires NGSolve "
                "PERIODIC point identifications between the two sector faces")
        available = tuple(mesh.GetBoundaries())
        missing = sorted(set(cyclic_periodic_boundaries) - set(available))
        if missing:
            raise ValueError(
                "vim.Solve: cyclic periodic boundary labels are missing: %s"
                % missing)
        boundary_counts = {name: 0 for name in cyclic_periodic_boundaries}
        boundary_vertices = {name: set() for name in cyclic_periodic_boundaries}
        for element in mesh.Elements(ng.BND):
            name = available[element.index]
            if name in boundary_counts:
                boundary_counts[name] += 1
                boundary_vertices[name].update(
                    int(vertex.nr) + 1 for vertex in element.vertices)
        for name in cyclic_periodic_boundaries:
            if boundary_counts[name] <= 0:
                raise ValueError(
                    "vim.Solve: cyclic periodic boundary %r has no facets"
                    % name)
        first, second = cyclic_periodic_boundaries
        if boundary_vertices[first] & boundary_vertices[second]:
            raise ValueError(
                "vim.Solve: cyclic periodic sector faces must have disjoint "
                "vertex sets")
        linked = {first: set(), second: set()}
        for endpoint_a, endpoint_b in mesh.ngmesh.GetIdentifications():
            endpoint_a = int(getattr(endpoint_a, "nr", endpoint_a))
            endpoint_b = int(getattr(endpoint_b, "nr", endpoint_b))
            if (endpoint_a in boundary_vertices[first]
                    and endpoint_b in boundary_vertices[second]):
                linked[first].add(endpoint_a)
                linked[second].add(endpoint_b)
            elif (endpoint_b in boundary_vertices[first]
                    and endpoint_a in boundary_vertices[second]):
                linked[first].add(endpoint_b)
                linked[second].add(endpoint_a)
        if any(linked[name] != boundary_vertices[name]
               for name in cyclic_periodic_boundaries):
            raise ValueError(
                "vim.Solve: NGSolve PERIODIC identifications must pair every "
                "vertex on the two cyclic_periodic_boundaries")
        periodic = ng.Periodic(base)
        active = periodic.FreeDofs()
        periodic_slave_dofs = tuple(
            index for index in range(periodic.ndof) if not active[index])
        if not periodic_slave_dofs:
            raise ValueError(
                "vim.Solve: mesh identifications did not identify any HDiv "
                "trace DoFs on the cyclic sector faces")
        base = ng.Compress(periodic, active)
    constrained_axes = tuple(sorted({int(axis) for axis, sign in image_planes
                                     if float(sign) > 0.0}))
    if not constrained_axes:
        return base, (), periodic_slave_dofs

    scale = max((abs(float(value))
                 for vertex in mesh.vertices for value in vertex.point), default=1.0)
    plane_tol = 128.0 * np.finfo(float).eps * max(scale, 1.0)
    axis_dofs = {axis: set() for axis in constrained_axes}
    for element in mesh.Elements(ng.BND):
        coordinates = np.asarray([mesh.vertices[v.nr].point for v in element.vertices], dtype=float)
        for axis in constrained_axes:
            if np.max(np.abs(coordinates[:, axis])) <= plane_tol:
                axis_dofs[axis].update(int(dof) for dof in base.GetDofNrs(element) if int(dof) >= 0)

    missing = ["xyz"[axis] for axis, dofs in axis_dofs.items() if not dofs]
    if missing:
        raise ValueError(
            "vim.Solve: image symmetry requires a boundary face on the coordinate plane(s) %s=0; "
            "no HDiv trace DOFs were found. The reduced mesh must end exactly on every image plane."
            % ",".join(missing))

    constrained = tuple(sorted(set().union(*axis_dofs.values())))
    active = ng.BitArray(base.ndof)
    active[:] = True
    for dof in constrained:
        active[dof] = False
    return ng.Compress(base, active), constrained, periodic_slave_dofs


def _auto_tet_jacobi_nface_threshold():
    raw = os.environ.get(_AUTO_TET_JACOBI_NFACE_ENV)
    if raw in (None, ""):
        return _AUTO_TET_JACOBI_NFACE_DEFAULT
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError("%s must be a positive integer (got %r)"
                         % (_AUTO_TET_JACOBI_NFACE_ENV, raw)) from exc
    if value <= 0:
        raise ValueError("%s must be a positive integer (got %r)"
                         % (_AUTO_TET_JACOBI_NFACE_ENV, raw))
    return value


def _resolve_highorder_preconditioner(preconditioner, *, nonlinear, nonlinear_solver, vertex_counts, n_face):
    """Resolve the production preconditioner policy for BDM1/BDM2 HDiv-VIM.

    The diagonal W+N preconditioner is dramatically faster for large hex/wedge nonlinear scaling runs because
    it avoids one PARDISO phase-33 mass solve per CG iteration.  Small tet problems still favor the exact
    mass-Riesz map.  Keep the policy local and explicit so the user can still force either branch.
    """
    if preconditioner != "auto":
        return preconditioner, "explicit:%s" % preconditioner
    if not nonlinear:
        return "mass-riesz", "auto:linear-mass-riesz"
    if nonlinear_solver != "energy-newton":
        return "mass-riesz", "auto:%s-mass-riesz" % nonlinear_solver
    if vertex_counts in ({6}, {8}):
        return "jacobi", "auto:hex-wedge-energy-newton-jacobi"
    threshold = _auto_tet_jacobi_nface_threshold()
    if int(n_face) >= threshold:
        return "jacobi", "auto:tet-energy-newton-jacobi-nface>=%d" % threshold
    return "mass-riesz", "auto:tet-energy-newton-mass-riesz-nface<%d" % threshold


def _clear_cpp_solve_timings():
    global _LAST_CPP_SOLVE_TIMINGS
    _LAST_CPP_SOLVE_TIMINGS = {}


def _capture_cpp_solve_timings(res):
    global _LAST_CPP_SOLVE_TIMINGS
    try:
        timings = dict(res.get("timings", {}) or {})
    except Exception:
        timings = {}
    clean = {}
    for key, value in timings.items():
        try:
            clean[str(key)] = float(value)
        except (TypeError, ValueError):
            pass
    if not _LAST_CPP_SOLVE_TIMINGS:
        _LAST_CPP_SOLVE_TIMINGS = clean
        return
    # Nonlinear solves call the C++ W-CG many times.  Sum timing/count-like
    # quantities so the final artifact reflects the full nonlinear solve, while
    # preserving state/identity diagnostics from the most recent inner solve.
    for key, value in clean.items():
        # "last_*" (like "hmatvec_last_*") are state of the most recent inner
        # solve -- converged flag, final residual -- not additive quantities.
        if (key.startswith("hmatvec_last_") or key.startswith("last_") or
                key == "mass_riesz_geometry_preconditioner"):
            _LAST_CPP_SOLVE_TIMINGS[key] = value
        else:
            _LAST_CPP_SOLVE_TIMINGS[key] = _LAST_CPP_SOLVE_TIMINGS.get(key, 0.0) + value


def _clear_nonlinear_solve_stats():
    global _LAST_NONLINEAR_SOLVE_STATS
    _LAST_NONLINEAR_SOLVE_STATS = {}


def _capture_nonlinear_solve_stats(stats):
    global _LAST_NONLINEAR_SOLVE_STATS
    _LAST_NONLINEAR_SOLVE_STATS = dict(stats or {})


def _i32(a):
    return np.ascontiguousarray(a, dtype=np.int32)


def _f64(a):
    return np.ascontiguousarray(a, dtype=np.float64)


def _configure_cpp_mass(H, mass, n_face):
    """Keep NGSolve sparse matrices native; accept SciPy only for external callers."""
    if sp.issparse(mass):
        coo = sp.coo_matrix(mass)
        H.configure_mass_matrix(
            _i32(coo.row), _i32(coo.col), _f64(coo.data), int(n_face))
    else:
        H.configure_mass_matrix_ngsolve(mass)


def _geometry_mass_apply(H, vector):
    """Apply the immutable NGSolve HDiv mass through the persistent C++ operator."""
    return np.asarray(H.apply_configured_geometry_mass(_f64(vector)), dtype=float)


def _h_solve_auto_prec(H, mass, n_face, inv_chi, rhs, tol, maxit, *, x0=None):
    if mass is not None:
        _configure_cpp_mass(H, mass, n_face)
    kwargs = {"x0": None if x0 is None else _f64(x0)}
    return H.solve_configured_linear_material_auto_prec(
        float(inv_chi), _f64(rhs), float(tol), int(maxit), **kwargs)


def _h_solve_mass_riesz(H, mass, n_face, inv_chi, rhs, tol, maxit, *,
                        symmetric=True, x0=None):
    if mass is not None:
        _configure_cpp_mass(H, mass, n_face)
    kwargs = {"x0": None if x0 is None else _f64(x0)}
    return H.solve_configured_linear_material_mass_riesz(
        float(inv_chi), _f64(rhs), float(tol), int(maxit), bool(symmetric), **kwargs)


def _resolve_gram_params(*, gram_eps, far_quad, ho_far_factor):
    """Resolve the production HDiv BDM1/BDM2 charge-Gram build defaults.

    ``far_quad`` controls the smooth, well-separated rule and ``ho_far_factor``
    controls the near/far boundary.  ``inf`` forces the all-high-order
    reference build.  There is no lower-order solver branch or legacy
    near-factor knob.
    """
    return {
        "eps": gram_eps if gram_eps is not None else 1e-10,
        "far_quad": far_quad if far_quad is not None else 3,
        "ho_far_factor": ho_far_factor if ho_far_factor is not None else 2.0,
    }


def _solve_linear_mass_riesz_cpp(H, n_face, h_ext, chi, tol, maxit):
    """DEFAULT uniform-chi linear solve: mass-Riesz-preconditioned CG ENTIRELY in C++ on the SPD +N system
    ((1/chi)M_mass + B^T G B) m = M_mass h_ext, with G applied via the EXACTLY-SYMMETRIC charge-Gram
    H-matvec (`matvec_sym`: the upper-triangular leaves define both triangles, so the operator is
    machine-symmetric -- ||G - G^T|| == 0 -- regardless of the per-block ACA truncation).  This is why CG
    is the default again (Sugahara 2026-06-27, "対称HACApKを実装しよう。CGがいいね"): the symmetric Gram
    makes CG robust BY CONSTRUCTION at all N -- it removes the asymmetry failure mode that motivated the
    earlier solver retreat (the spurious antisymmetric part of the GENERAL ACA matvec).  Mass-Riesz precond
    via a PARDISO SPD factor of the HDiv mass (eigenvalues vs M_mass are (1/chi)+d, d in [0,1], bounded ->
    ~3-5x fewer iters than diagonal Jacobi).  The whole Krylov loop (O(N log N) symmetric H-matvec +
    per-iteration PARDISO mass solve + vector ops) runs in C++ -- no Python per-iteration glue, no splu.
    The symmetric matvec is also ~1.4x FASTER than the general one (it skips the lower-triangle leaves).
    `H.solve_configured_linear_material_mass_riesz(..., symmetric=True)` is the default; pass
    symmetric=False only to
    cross-check against the general (asymmetric) matvec."""
    inv_chi = 1.0 / float(chi)
    rhs = _geometry_mass_apply(H, h_ext)
    res = _h_solve_mass_riesz(
        H, None, int(n_face), inv_chi, rhs, tol, int(maxit))
    _capture_cpp_solve_timings(res)
    iters = int(res["iters"])
    if iters >= int(maxit):                      # fail-loud (No-Fallbacks): never return a non-converged M
        raise RuntimeError(
            "vim.Solve (symmetric mass-riesz CG): did NOT converge in %d iters (n_face=%d).  The "
            "operator is the EXACTLY-symmetric SPD +N system, so CG should converge -- a non-convergence "
            "here means an ill-conditioned material/mesh. Tighten gram_eps or raise maxit."
            % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def _solve_linear_jacobi_cpp(H, system_mass, n_face, h_ext, inv_chi, tol, maxit):
    """Diagnostic linear solve: exact Jacobi diagonal of (inv_chi*system_mass + N), then C++ CG.

    This is intentionally not the default.  It gives the hex/H-matrix benchmark a direct apples-to-apples
    comparison against the production mass-Riesz preconditioner while keeping the same symmetric charge-Gram
    H-matvec and C++ Krylov loop.
    """
    rhs = _geometry_mass_apply(H, h_ext)
    res = _h_solve_auto_prec(
        H, system_mass, int(n_face), float(inv_chi), rhs, tol, int(maxit))
    _capture_cpp_solve_timings(res)
    iters = int(res["iters"])
    if iters >= int(maxit):
        raise RuntimeError(
            "vim.Solve (symmetric Jacobi-CG): did NOT converge in %d iters (n_face=%d).  "
            "Use preconditioner='mass-riesz' for the production Riesz-map preconditioner."
            % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def _solve_linear_W_cpp(H, W, n_face, h_ext, tol, maxit):
    """PER-REGION linear solve ENTIRELY in C++: the SYMMETRIC Galerkin system (M_{1/chi} + N) m = M_mass h_ext
    by symmetric mass-Riesz CG.  W = M_{1/chi} = INT (1/chi(x)) u.v dx is the SYSTEM mass; the immutable
    geometry mass M_mass is the Riesz preconditioner.  Passing W as the 'mass' COO with inv_chi=1.0 makes
    the C++ kernel (`solve_configured_linear_material_mass_riesz`, symmetric=True) compute
    A = W + B^T G B while preconditioning with M_mass^{-1} (PARDISO).  Keeping the Riesz map material-
    independent is the same stable C++ CG contract as the uniform path and lets nonlinear tangent updates
    reuse its factor.  The whole Krylov loop runs in C++ (symmetric charge-Gram H-matvec + PARDISO solve
    + vector ops).
    Python declares W and the geometric mass as NGSolve forms; NGSolve assembles them in C++, pybind extracts
    their native sparse matrices directly, and the persistent C++ operator applies the geometric-mass RHS.
    The
    symmetric 1/chi-weighted Galerkin form is CG-able (W, N both symmetric), and for a UNIFORM region it is
    bit-identical to the scalar +N system."""
    rhs = _geometry_mass_apply(H, h_ext)
    res = _h_solve_mass_riesz(
        H, W, int(n_face), 1.0, rhs, tol, int(maxit))
    _capture_cpp_solve_timings(res)
    iters = int(res["iters"])
    if iters >= int(maxit):                      # fail-loud (No-Fallbacks)
        raise RuntimeError(
            "vim.Solve (per-region symmetric mass-riesz CG): did NOT converge in %d iters "
            "(n_face=%d).  The (M_{1/chi} + N) operator is SPD, so a non-convergence means an ill-"
            "conditioned material/mesh; tighten gram_eps or raise maxit." % (maxit, n_face))
    return np.asarray(res["m"], float), iters


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, B_r=None, bh_table=None,
                     magnetization_sources=None, magnets=None,
                     image=None, image_cyclic=None, image_cyclic_alternating=False,
                     cyclic_periodic_boundaries=None,
                     gram_eps=None, leaf=64, eta=2.0, far_quad=None, tol=1e-8,
                     maxit=4000, nl_maxit=300, nl_tol=1e-6,
                     nonlinear_solver="energy-newton",
                     preconditioner="auto",
                     linear_solver="auto", order=1,
                     curve_order=None, curve_gauss=8, ho_far_factor=None,
                     newton_inner_tol="auto", newton_warmstart="linear",
                     newton_continuation=1, newton_reuse_tangent_steps=1,
                     newton_cg_x0=False, gram_backend="hmat",
                     exact_dense_memory_mb=None, _operator_cache=None):
    """HDiv-type VIM soft-iron demag solve (the +N physical material system).

    ``leaf=64`` is the production ChargeGram H-matrix default. It is
    calibrated for high-order charge rows sharing clustering centres; a
    smaller explicit value requires raw-Gram and PSD validation for the
    selected mesh family.

    Material spec (EXACTLY ONE):
      mu_r     : float > 1 -> LINEAR isotropic soft iron (ONE region), OR a dict {material_name: mu_r}
                 for PER-REGION linear soft iron (multiple iron grades; each mu_r > 1).  N = B^T G B is
                 geometry-only, so per-region enters ONLY through the chi-weighted HDiv mass (form 1).
      B_r      : optional vector remanent flux density (T), supplied with scalar ``mu_r``.  This selects
                 the linear-recoil permanent-magnet law B = mu0*mu_r*H + B_r.  ``B_r`` may be a constant
                 vector or a spatial NGSolve CoefficientFunction within one conforming magnet body.
                 A normal-discontinuous segmented magnet requires separate body spaces so its internal
                 surface charge is retained.  ``mu_r`` is the recoil permeability and must be > 1; the
                 rigid mu_r=1 limit is ``vim.MagnetizationSource``.
      bh_table : [[H,B], ..] (A/m, T; the MatSatIsoTab data) -> NONLINEAR isotropic soft iron (ONE
                 region), OR a dict {material_name: [[H,B], ..]} for PER-REGION nonlinear soft iron.
                 N = B^T G B is geometry-only, so per-region nonlinear enters ONLY through the per-element
                 constitutive law + warmstart.
      magnetization_sources : one or more ``vim.MagnetizationSource`` objects.  Each fixed/given
                  magnetization is L2-projected in its OWN HDiv space and contributes its immutable C++
                  field to the iron weak form.  The PM coefficients are never included among the unknowns;
                  separate PM and iron spaces preserve their physical normal jump, including touching
                  interfaces.
    H_ext      : NGSolve CoefficientFunction, the applied field (A/m) -- uniform, analytic, or a coil's
                 Biot-Savart field rad.RadiaField(coil,'h').  Required unless ``B_r``, a magnetization
                  source, or the planar magnets path supplies the drive.
    cyclic_periodic_boundaries : the two named azimuthal cut faces of a
                 connected pure-HEX sector.  The mesh must carry NGSolve
                 PERIODIC point identifications and ``image_cyclic=N`` must
                 be set.  The solve uses ``Compress(Periodic(HDiv))`` and
                 removes the paired seam faces from the physical charge skin;
                 omitting this argument on a mesh with periodic/cyclic labels
                 fails loudly instead of leaving artificial seam charge.
    near/far Gram-build tuning:
      ho_far_factor -- the HDiv near/far separation threshold (pass inf to force the all-high-order
                       reference build).
    Returns dict: M (n_el,3) per-element magnetization, M_avg (3,), iters, demag (Rayleigh factor),
    ndof, n_el, n_charge, nonlinear(bool).  The caller must open `with ng.TaskManager():`.
    """
    from ._magnetization_source import MagnetizationSource
    if magnetization_sources is None:
        sources = ()
    elif isinstance(magnetization_sources, MagnetizationSource):
        sources = (magnetization_sources,)
    else:
        try:
            sources = tuple(magnetization_sources)
        except TypeError as exc:
            raise TypeError("vim.Solve: magnetization_sources must be a MagnetizationSource or iterable") from exc
        if not all(isinstance(source, MagnetizationSource) for source in sources):
            raise TypeError("vim.Solve: every magnetization_sources item must be a vim.MagnetizationSource")
    linear_recoil_pm = B_r is not None
    if linear_recoil_pm:
        if bh_table is not None:
            raise ValueError("vim.Solve: B_r selects the linear-recoil permanent-magnet law and cannot "
                             "be combined with bh_table")
        if isinstance(mu_r, dict):
            raise NotImplementedError(
                "vim.Solve: B_r currently requires one scalar recoil mu_r; per-region recoil materials "
                "must be solved as separate permanent-magnet bodies")
        try:
            recoil_mu_r = float(mu_r)
        except (TypeError, ValueError) as exc:
            raise ValueError("vim.Solve: B_r requires scalar recoil mu_r > 1") from exc
        if not np.isfinite(recoil_mu_r) or recoil_mu_r <= 1.0:
            raise ValueError(
                "vim.Solve: linear-recoil B_r requires mu_r > 1; use vim.MagnetizationSource for the "
                "rigid mu_r=1 fixed-magnetization limit")
        dim = int(mesh.dim)
        try:
            uniform_B_r = np.asarray(B_r, dtype=float)
        except (TypeError, ValueError):
            uniform_B_r = None
        if (uniform_B_r is not None and uniform_B_r.shape == (dim,)
                and np.isfinite(uniform_B_r).all()):
            B_r_cf = ng.CoefficientFunction(tuple(float(value) for value in uniform_B_r))
        elif uniform_B_r is not None and uniform_B_r.shape == (dim,):
            raise ValueError("vim.Solve: B_r vector must contain only finite values")
        elif getattr(B_r, "dim", None) == dim:
            B_r_cf = B_r
        else:
            raise ValueError(
                "vim.Solve: B_r must be a length-%d vector or a %d-component NGSolve "
                "CoefficientFunction (got %r)" % (dim, dim, B_r))
    if H_ext is None:
        if linear_recoil_pm or sources or magnets is not None:
            H_ext = ng.CoefficientFunction((0.0,) * int(mesh.dim))
        else:
            raise ValueError("vim.Solve: H_ext (applied-field CoefficientFunction) is required when no "
                             "magnetization source is supplied")
    if getattr(H_ext, "dim", None) != int(mesh.dim):
        raise ValueError("vim.Solve: H_ext must be a %d-component NGSolve CoefficientFunction"
                         % int(mesh.dim))
    if linear_recoil_pm:
        # M = chi*H + B_r/mu0 is exactly the existing symmetric linear HDiv system with
        # H_ext shifted by B_r/(mu0*chi).  The unknown and returned field remain total M.
        H_ext = H_ext + B_r_cf / (_MU0 * (recoil_mu_r - 1.0))
    # ---- HDiv-VIM scope: BDM1/BDM2 on pure TET/HEX/WEDGE/planar meshes ----
    # The tetrahedral monomial-charge Gram, IMA fold, and field evaluator are exact
    # for orders 1 and 2.  The planar log kernel also supports BDM1/Q2 and BDM2/Q3;
    # specialized HEX/WEDGE kernels use the same order without fallback.
    # Mapped HEX BDM2 uses the cancellation-preserving composite charge rule
    # documented on the production high-order path below.
    order = int(order)
    if curve_order is None and mesh.dim == 3 and mesh.GetCurveOrder() >= 2:
        curve_order = int(mesh.GetCurveOrder())
    _vtx = _volume_vertex_counts(mesh)
    validate_hdiv_configuration(mesh.dim, _vtx, order, mesh.GetCurveOrder())
    # IMA mirror symmetry is wired for flat/curved pure-TET/HEX/WEDGE BDM1/BDM2.
    # Mixed and pyramid cases fail loud downstream instead of silently dropping the image.
    if mesh.dim == 2:
        if sources:
            raise NotImplementedError(
                "vim.Solve (2D): use the established magnets=[(mesh, M), ...] planar source path; "
                "vim.MagnetizationSource is the 3D prescribed-source API.")
        image_masks, image_signs = [], []
        if image is not None:
            for axes, sign in _image.image_group(_image.parse_image_string(image)):
                if any(axis >= 2 for axis in axes):
                    raise ValueError("vim.Solve (2D): image may contain x/y planes only.")
                image_masks.append(int(sum(1 << axis for axis in axes)))
                image_signs.append(float(sign))
        # ---- PLANAR (2D motor cross-section) branch: matrix-free C++ charge Gram + mass-Riesz CG ----
        # The 2D layer supports the core single-region surface: mu_r (linear) / bh_table
        # (nonlinear) + H_ext.  The 3D-only knobs must stay at their defaults -- fail loud.
        if linear_solver != "auto":
            raise ValueError("vim.Solve (2D): linear_solver must be 'auto' (the planar layer "
                             "uses C++ mass-Riesz CG; got %r)" % (linear_solver,))
        if preconditioner == "auto":
            preconditioner = "mass-riesz"
        if preconditioner != "mass-riesz":
            raise ValueError("vim.Solve (2D): preconditioner must be 'mass-riesz' (got %r)"
                             % (preconditioner,))
        for _nm, _val in (("gram_eps", gram_eps), ("far_quad", far_quad), ("ho_far_factor", ho_far_factor),
                          ("curve_order", curve_order), ("exact_dense_memory_mb", exact_dense_memory_mb)):
            if _val is not None:
                raise ValueError("vim.Solve (2D): %s is a 3D knob; the 2D Gram parameters "
                                 "are fixed by its own gates (got %r)" % (_nm, _val))
        if gram_backend != "hmat":
            raise ValueError("vim.Solve (2D): gram_backend must be 'hmat' (got %r)"
                             % (gram_backend,))
        from ._vim2d import solve_planar_demag
        result = solve_planar_demag(
            mesh, mu_r=mu_r, H_ext=H_ext, bh_table=bh_table, magnets=magnets,
            order=order, eta=eta, cg_tol=tol, cg_maxit=maxit,
            nl_tol=nl_tol, nl_maxit=nl_maxit,
            image_masks=image_masks, image_signs=image_signs)
        result["image"] = image
        if linear_recoil_pm:
            result["permanent_magnet_model"] = "linear-recoil"
            result["permanent_magnet_level"] = 2
            result["recoil_mu_r"] = recoil_mu_r
            result["B_r_supplied"] = True
            result["_B_r"] = B_r_cf
        return result
    if magnets is not None:
        raise NotImplementedError(
            "vim.Solve: magnets= is the 2D planar source API; in 3D construct "
            "vim.MagnetizationSource and pass magnetization_sources=[...].")
    for source in sources:
        if source.mesh is mesh:
            raise ValueError(
                "vim.Solve: a prescribed magnetization source and soft iron must use separate mesh "
                "objects/HDiv spaces; using the same mesh would erase the PM/iron normal jump.")
        H_ext = H_ext + source.field_cf
    if linear_solver not in _LINEAR_SOLVERS:
        raise ValueError("vim.Solve: linear_solver must be one of %s (got %r)"
                         % (sorted(_LINEAR_SOLVERS), linear_solver))
    if nonlinear_solver not in _NONLINEAR_SOLVERS:
        raise ValueError("vim.Solve: nonlinear_solver must be one of %s (got %r)"
                         % (sorted(_NONLINEAR_SOLVERS), nonlinear_solver))
    if preconditioner not in _PRECONDITIONERS:
        raise ValueError("vim.Solve: preconditioner must be one of %s (got %r)"
                         % (sorted(_PRECONDITIONERS), preconditioner))
    if (mu_r is None) == (bh_table is None):
        raise ValueError("vim.Solve: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")
    if bh_table is not None:
        _validate_bh_table(bh_table)

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

    # ---- BDM1/BDM2 material solve ----
    # The per-element change-of-basis in `_vim._charge_basis` (2026-06-28, [[hdiv-highorder-material-solve-wrong]])
    # makes the order-p demag operator N = B^T G B valid (eig(M_mass^-1 N) in [0,1]; per-element M p-converges).
    # LINEAR (uniform-scalar OR per-region dict) mu_r AND flat/curved NONLINEAR (bh_table) are wired via the same
    # all-C++ symmetric mass-Riesz CG / energy-Newton on each supported pure topology.
    # Pure-HEX BDM1 and BDM2 material solves cover affine and mapped Q2
    # geometry.  Mapped BDM2 shape derivatives remain a separate fail-loud
    # limitation until the composite rule is differentiated consistently.
    # Golden: validation_test/feec/test_hdiv_vim_demag_solve*, test_hdiv_vim_highorder_cpp.py.
    result = _solve_highorder(mesh, int(order), mu_r, bh_table, H_ext, image, linear_solver,
                              gram_eps, leaf, eta, far_quad, tol, maxit,
                              curve_order, curve_gauss, ho_far_factor, nl_maxit, nl_tol,
                              nonlinear_solver, preconditioner, newton_inner_tol, newton_warmstart,
                              newton_continuation, newton_reuse_tangent_steps, newton_cg_x0,
                              vertex_counts=_vtx, magnetization_sources=sources,
                              gram_backend=gram_backend,
                              exact_dense_memory_mb=exact_dense_memory_mb,
                              operator_cache=_operator_cache,
                              image_cyclic=image_cyclic,
                              image_cyclic_alternating=image_cyclic_alternating,
                              cyclic_periodic_boundaries=cyclic_periodic_boundaries)
    if linear_recoil_pm:
        result["permanent_magnet_model"] = "linear-recoil"
        result["permanent_magnet_level"] = 2
        result["recoil_mu_r"] = recoil_mu_r
        result["B_r_supplied"] = True
        result["_B_r"] = B_r_cf
    return result


def _solve_highorder(mesh, order, mu_r, bh_table, H_ext, image, linear_solver,
                     gram_eps, leaf, eta, far_quad, tol, maxit,
                     curve_order=None, curve_gauss=8, ho_far_factor=None, nl_maxit=300, nl_tol=1e-6,
                      nonlinear_solver="energy-newton", preconditioner="auto",
                      newton_inner_tol="auto", newton_warmstart="linear", newton_continuation=1,
                      newton_reuse_tangent_steps=1, newton_cg_x0=False, vertex_counts=None,
                      magnetization_sources=(), gram_backend="hmat",
                      exact_dense_memory_mb=None, operator_cache=None,
                      image_cyclic=None, image_cyclic_alternating=False,
                      cyclic_periodic_boundaries=None):
    """BDM1/BDM2 HDiv soft-iron demag solve.  The order-p charge-Gram demag operator N = B^T G B is
    a VALID demag operator since the per-element change-of-basis fix (2026-06-28,
    [[hdiv-highorder-material-solve-wrong]]): eig(M_mass^-1 N) in [0,1] and the material solve p-converges
    (no 2x/4x blow-up).  Supports the LINEAR (uniform-scalar OR per-region dict) mu_r case via the SAME
    all-C++ symmetric mass-Riesz CG as the production path; unsupported combinations fail loud (No-Fallbacks).
    Mapped/non-affine HEX BDM2 uses a complete-host tensor source rule for
    smooth pairs and target-anchored radial Duffy integration for adjacent
    volume-boundary pairs.  The charge quotient remains B^T G B, so loop
    nullspaces are annihilated to roundoff on affine and mapped meshes alike.
    The CALLER opens `with ng.TaskManager():` (same contract as vim.Solve)."""
    t_total = time.perf_counter()
    vertex_counts = (_volume_vertex_counts(mesh) if vertex_counts is None
                     else frozenset(vertex_counts))
    # IMA mirror symmetry: wired for flat/Curve(2) pure-TET (C++ highorder QuadDotRefl->PhiInner) AND pure-HEX /
    # pure-WEDGE (the C++ QuadBlockHex/Wedge(mask) reflected block) paths -- the Gram folds the mirror-image
    # charge interactions so a reduced 1/2,1/4,1/8 model reproduces the full model.  The same fold is used
    # on P2-curved TET/HEX/WEDGE geometry; MIXED / pyramid still fail loud.
    image_masks, image_signs = [], []
    image_planes = []
    if image is not None:
        _ivtx = vertex_counts
        if _ivtx not in ({4}, {8}, {6}):
            raise NotImplementedError(
                "vim.Solve: IMA image symmetry is wired for flat/Curve(2) pure-TET / pure-HEX / pure-WEDGE "
                "BDM1/BDM2 Gram; MIXED / pyramid reduced models are not supported.  Got vertex counts %s."
                % sorted(_ivtx))
        _planes = _image.parse_image_string(image)
        image_planes = list(_planes)
        # (2026-07-05) hex/wedge IMA now handles ANTISYMMETRIC (negative-sign, field-PERPENDICULAR) planes too:
        # the on-plane cut-face self-term is computed with the EXACT self-radial in the reflected block (the
        # QuadBlockHex/Wedge "R(host)==host -> self_pair" fix), so the large perpendicular cut-face charge
        # cancels exactly for sign -1 instead of the earlier ~1.5% hex / ~29% wedge quadrature residual.
        for axes, sign in _image.image_group(_planes):
            image_masks.append(int(sum(1 << a for a in axes)))
            image_signs.append(float(sign))
    # CYCLIC (N-fold rotational) reduction: solve ONE sector and let the Gram fold in the other N-1 poles as
    # rotated images about +z.  Unlike a mirror this is not an involution, so the C++ side maps eval points
    # through the INVERSE rotation; and unlike an infinite translational array the finite rotational sum is
    # unconditionally well posed (an infinite dipole lattice sum is only conditionally convergent -- that
    # shape dependence is the demagnetizing-factor phenomenon).  Charges are SCALARS under a rotation that
    # carries the magnetization with the geometry, so alternating N/S poles are just signs (-1)^k.
    image_rot_angle = []
    if isinstance(cyclic_periodic_boundaries, str):
        raise TypeError(
            "vim.Solve: cyclic_periodic_boundaries must be a pair of labels, "
            "not one string")
    cyclic_periodic_boundaries = (() if cyclic_periodic_boundaries is None
                                  else tuple(str(name) for name in
                                             cyclic_periodic_boundaries))
    labeled_periodic_boundaries = tuple(
        name for name in mesh.GetBoundaries()
        if "periodic" in name.lower() or "cyclic" in name.lower())
    if image_cyclic is not None:
        n_fold = int(image_cyclic)
        if n_fold < 2:
            raise ValueError("vim.Solve: image_cyclic must be an N-fold count >= 2; got %r" % (image_cyclic,))
        if image_cyclic_alternating and n_fold % 2:
            raise ValueError(
                "vim.Solve: image_cyclic_alternating needs an EVEN pole count (the (-1)^k pattern must "
                "close around the ring); got N=%d" % n_fold)
        if image_masks:
            raise NotImplementedError(
                "vim.Solve: combining image= mirror planes with image_cyclic= is not wired yet; "
                "pass one or the other")
        for k in range(1, n_fold):
            image_masks.append(0)                                  # pure rotation, no mirror
            image_signs.append(-1.0 if (image_cyclic_alternating and k % 2) else 1.0)
            image_rot_angle.append(2.0*math.pi*k/n_fold)
    if cyclic_periodic_boundaries:
        if image_cyclic is None:
            raise ValueError(
                "vim.Solve: cyclic_periodic_boundaries requires image_cyclic=N")
        if len(cyclic_periodic_boundaries) != 2 or (
                cyclic_periodic_boundaries[0] == cyclic_periodic_boundaries[1]):
            raise ValueError(
                "vim.Solve: cyclic_periodic_boundaries must name two distinct "
                "rotation-related sector faces")
        if image_cyclic_alternating:
            raise NotImplementedError(
                "vim.Solve: connected antiperiodic HDiv sectors are not yet "
                "wired; FFAG guide-field sectors use periodic traces")
        if vertex_counts != {8}:
            raise NotImplementedError(
                "vim.Solve: connected cyclic periodic traces are currently "
                "wired for pure HEX BDM1/BDM2")
    elif image_cyclic is not None and labeled_periodic_boundaries:
        raise ValueError(
            "vim.Solve: image_cyclic mesh has azimuthal periodic/cyclic "
            "boundary labels %s; pass cyclic_periodic_boundaries=(min,max) "
            "so the HDiv traces are identified and the seam charge is removed"
            % (labeled_periodic_boundaries,))
    # The flat nonlinear path uses the same high-order Gram as the curved path.  The symmetric energy-Newton
    # is Gram-agnostic and consumes only H.matvec plus the C++ mass-Riesz solve.
    if int(order) > 2:
        # order<=2 uses the EXACT analytic-moment charge potential (machine precision).  For order>=3 the C++
        # Gram falls back to the Duffy singular quadrature (PhiInner -> PhiAtHO_Duffy), which is ~1e-3 accurate
        # -- fine for curved-panel field evaluation, but NOT for the order>=3 MATERIAL solve: the ill-
        # conditioned high-degree monomial basis (cond(B)^2 in N=B^T G B) amplifies the ~1e-3 entry error so
        # the demag spectrum escapes [0,1].  A clean order>=3 material solve needs machine-precision entries
        # (the analytic moments extended with TetMoment2 / degree-3 surface moments), not the Duffy.  Fail loud
        # (No-Fallbacks) until that lands.  [[hdiv-vim-sauter-schwab-cg]]
        raise NotImplementedError(
            "vim.Solve: order>2 material solve is not yet production-clean -- order in {1,2} is exact "
            "(analytic moments); the order>=3 Duffy quadrature is only ~1e-3 and the ill-conditioned "
            "high-degree basis makes the demag spectrum leave [0,1]. Use order in {1,2}.")
    _gp = _resolve_gram_params(gram_eps=gram_eps, far_quad=far_quad, ho_far_factor=ho_far_factor)
    eff_eps = _gp["eps"]; eff_far = _gp["far_quad"]; eff_hofar = _gp["ho_far_factor"]
    t_before_fes = time.perf_counter()
    operator_reused = operator_cache is not None
    if operator_reused:
        if not isinstance(operator_cache, dict):
            raise TypeError("vim.Solve: internal operator cache must be owned by vim.HDivSolver")
        expected = dict(mesh=mesh, order=int(order), curve_order=curve_order, image=image,
                        image_cyclic=image_cyclic,
                        image_cyclic_alternating=bool(image_cyclic_alternating),
                        cyclic_periodic_boundaries=cyclic_periodic_boundaries,
                        vertex_counts=frozenset(vertex_counts), gram_backend=gram_backend,
                        exact_dense_memory_mb=exact_dense_memory_mb,
                        nonlinear=bool(bh_table is not None))
        for key, value in expected.items():
            cached = operator_cache.get(key)
            matches = (cached is value) if key == "mesh" else (cached == value)
            if not matches:
                raise ValueError(
                    "vim.Solve: prepared operator mismatch for %s (cached=%r, requested=%r)"
                    % (key, cached, value))
        fes = operator_cache["fes"]
        H = operator_cache["charge_gram"]
        n_face = int(operator_cache["n_face"])
        n_el = int(operator_cache["n_el"])
        n_charge = int(operator_cache["n_charge"])
        symmetry_constrained_dofs = tuple(
            operator_cache.get("symmetry_constrained_dofs", ()))
        periodic_slave_dofs = tuple(
            operator_cache.get("periodic_slave_dofs", ()))
        t_before_charge_gram = t_before_fes
        t_after_charge_gram = t_before_fes
        charge_build_timings = {}
    else:
        if curve_order is not None:
            # CURVED (isoparametric P2) demag solve: curve the geometry, then the curved-Duffy charge Gram.
            if int(curve_order) != 2:
                raise NotImplementedError("vim.Solve: only curve_order=2 (isoparametric P2) is wired.")
            _curve_mesh(mesh, int(curve_order))
        (fes, symmetry_constrained_dofs,
         periodic_slave_dofs) = _hdiv_space_with_image_constraints(
            mesh, order, image_planes, cyclic_periodic_boundaries)
        t_before_charge_gram = time.perf_counter()
        B, H, M_mass = build_charge_gram(
            fes, eps=eff_eps, leafsize=leaf, eta=eta,
            far_quad=eff_far, ho_far_factor=eff_hofar,
            curve_order=(int(curve_order) if curve_order is not None else None),
            curve_gauss=int(curve_gauss), nonlinear=bh_table is not None,
            image_masks=image_masks, image_signs=image_signs,
            image_rot_angle=image_rot_angle,
            excluded_boundaries=cyclic_periodic_boundaries,
            gram_backend=gram_backend,
            exact_dense_memory_mb=exact_dense_memory_mb,
            _materialize_mass=False)
        t_after_charge_gram = time.perf_counter()
        charge_build_timings = dict(getattr(build_charge_gram, "last_timings", {}) or {})
        n_face = fes.ndof
        n_el = mesh.GetNE(ng.VOL)
        n_charge = B.shape[0]
        del B, M_mass
    preconditioner_requested = preconditioner
    preconditioner, preconditioner_policy = _resolve_highorder_preconditioner(
        preconditioner,
        nonlinear=bh_table is not None,
        nonlinear_solver=nonlinear_solver,
        vertex_counts=vertex_counts,
        n_face=n_face,
    )
    # True weak L2 projection of arbitrary applied fields, including the native
    # C++ fields of prescribed HDiv magnetization sources.  GridFunction.Set is
    # interpolation; assembling the LinearForm and applying the already-pinned
    # mass Riesz map preserves the NGSolve Galerkin contract at PM/iron gaps and
    # touching interfaces.
    source_rhs = ng.LinearForm(fes)
    source_rhs += H_ext * fes.TestFunction() * ng.dx
    source_rhs.Assemble()
    h_ext = np.asarray(H.apply_configured_mass_riesz(
        _f64(source_rhs.vec.FV().NumPy())), dtype=float)
    t_after_projection = time.perf_counter()

    def N_apply(v):
        return np.asarray(H.apply_configured_demag(_f64(v), True), float)

    if operator_reused:
        D = float(operator_cache["demag"])
        cached_stats = operator_cache.get("hmat_stats")
        hmat_stats = None if cached_stats is None else dict(cached_stats)
    else:
        gfMu = ng.GridFunction(fes)
        gfMu.Set(ng.CoefficientFunction((0, 0, 1)))
        mu = gfMu.vec.FV().NumPy().copy()
        denom = float(mu @ _geometry_mass_apply(H, mu))
        D = float((mu @ N_apply(mu)) / denom)
        hmat_stats = (dict(H.stats()) if gram_backend == "hmat" and hasattr(H, "stats")
                      else None)
        if hmat_stats is not None and hasattr(H, "hex_state_breakdown"):
            try:
                _hex_diag = H.hex_state_breakdown()
                if "hexUniformAffineCells" in _hex_diag:
                    hmat_stats["hex_uniform_affine_cells"] = bool(_hex_diag["hexUniformAffineCells"])
                if "hexUniformTransHosts" in _hex_diag:
                    hmat_stats["hex_uniform_trans_hosts"] = bool(_hex_diag["hexUniformTransHosts"])
            except Exception:
                pass
        operator_cache = dict(
            mesh=mesh, order=int(order), curve_order=curve_order, image=image,
            image_cyclic=image_cyclic,
            image_cyclic_alternating=bool(image_cyclic_alternating),
            cyclic_periodic_boundaries=cyclic_periodic_boundaries,
            vertex_counts=frozenset(vertex_counts), gram_backend=gram_backend,
            exact_dense_memory_mb=exact_dense_memory_mb,
            nonlinear=bool(bh_table is not None),
            fes=fes, charge_gram=H, n_face=int(n_face), n_el=int(n_el),
            n_charge=int(n_charge), demag=float(D),
            hmat_stats=(None if hmat_stats is None else dict(hmat_stats)),
            symmetry_constrained_dofs=tuple(symmetry_constrained_dofs),
            periodic_slave_dofs=tuple(periodic_slave_dofs),
            build_timings=dict(charge_build_timings))
    t_after_demag_probe = time.perf_counter()

    t_solve = time.perf_counter()
    _clear_cpp_solve_timings()
    _clear_nonlinear_solve_stats()
    setup_wall_s = t_solve - t_total
    if bh_table is None and not isinstance(mu_r, dict):
        # Per-region and nonlinear solves replace the mutable material mass on
        # the persistent C++ operator.  A later scalar-mu load case needs the
        # immutable geometry mass again.  The C++ call is a no-op while the
        # geometry mass is already active, preserving factor reuse for normal
        # uniform load sweeps.
        H.restore_geometry_mass_matrix()
    if bh_table is not None:
        if preconditioner != "mass-riesz" and nonlinear_solver != "energy-newton":
            raise NotImplementedError("vim.Solve: nonlinear_solver=%r is wired only with "
                                      "preconditioner='mass-riesz' for now" % (nonlinear_solver,))
        if nonlinear_solver == "picard-mass-riesz":
            m, iters = _solve_nonlinear_picard_mass_riesz_cpp(mesh, fes, bh_table, H, n_face, h_ext,
                                                              tol, maxit, nl_maxit, nl_tol)
            solver_used = "picard-mass-riesz-cpp"
        elif nonlinear_solver == "picard-energy":
            m0, it0 = _solve_nonlinear_picard_mass_riesz_cpp(mesh, fes, bh_table, H, n_face, h_ext,
                                                             tol, maxit, min(nl_maxit, 12),
                                                             max(nl_tol, 1e-2),
                                                             require_convergence=False)
            m, it1 = _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, n_face, h_ext,
                                                 tol, maxit, nl_maxit, nl_tol, m0=m0,
                                                 inner_tol=newton_inner_tol,
                                                 warmstart="linear",
                                                 continuation_steps=newton_continuation,
                                                 reuse_tangent_steps=newton_reuse_tangent_steps,
                                                 cg_x0=bool(newton_cg_x0),
                                                 inner_preconditioner=preconditioner)
            iters = int(it0) + int(it1)
            solver_used = "picard-energy-cpp"
        else:                                                   # energy-Newton on the same charge Gram
            if newton_warmstart not in ("linear", "picard", "none"):
                raise ValueError("vim.Solve: newton_warmstart must be 'linear', 'picard', or 'none' "
                                 "(got %r)" % (newton_warmstart,))
            m0 = None
            warm_iters = 0
            if newton_warmstart == "picard":
                m0, warm_iters = _solve_nonlinear_picard_mass_riesz_cpp(
                    mesh, fes, bh_table, H, n_face, h_ext, tol, maxit,
                    min(nl_maxit, 8), max(nl_tol, 5e-3), require_convergence=False)
            m, iters = _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, n_face, h_ext,
                                                   tol, maxit, nl_maxit, nl_tol, m0=m0,
                                                   inner_tol=newton_inner_tol,
                                                   warmstart=newton_warmstart,
                                                   continuation_steps=newton_continuation,
                                                   reuse_tangent_steps=newton_reuse_tangent_steps,
                                                   cg_x0=bool(newton_cg_x0),
                                                   inner_preconditioner=preconditioner)
            iters = int(warm_iters) + int(iters)
            solver_used = "energy-newton-cpp"
    elif isinstance(mu_r, dict):                                # per-region linear: W = 1/chi-weighted HDiv mass
        W = _build_invchi_mass(mesh, fes, mu_r, n_face)
        if preconditioner == "jacobi":
            m, iters = _solve_linear_jacobi_cpp(H, W, n_face, h_ext, 1.0, tol, maxit)
            solver_used = "jacobi-cg"
        else:
            m, iters = _solve_linear_W_cpp(H, W, n_face, h_ext, tol, maxit)
            solver_used = "mass-riesz-cg"
    else:                                                       # uniform-scalar linear
        chi = float(mu_r) - 1.0
        if chi <= 0.0:
            raise ValueError("vim.Solve: mu_r must be > 1 (got %r)" % (mu_r,))
        if preconditioner == "jacobi":
            m, iters = _solve_linear_jacobi_cpp(H, None, n_face, h_ext, 1.0 / chi, tol, maxit)
            solver_used = "jacobi-cg"
        else:
            m, iters = _solve_linear_mass_riesz_cpp(H, n_face, h_ext, chi, tol, maxit)
            solver_used = "mass-riesz-cg"
    solve_wall_s = time.perf_counter() - t_solve
    cpp_solve_timings = dict(_LAST_CPP_SOLVE_TIMINGS)
    t_post = time.perf_counter()

    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)
    M_el = np.vstack([
        np.asarray(ng.Integrate(gfM[i], mesh, element_wise=True), float) / vol_el
        for i in range(3)
    ]).T.copy()
    vol = float(np.sum(vol_el))
    M_avg_reduced = np.average(M_el, axis=0, weights=vol_el)   # raw REDUCED-domain average
    # With image= (a symmetry-reduced solve), the raw average above is over the REDUCED domain only.  A
    # magnetization component that is ODD across an image mirror plane (e.g. Mx across a '+x' plane) integrates
    # to a nonzero ONE-SIDED value over the reduced half/quarter/octant, but that value CANCELS exactly against
    # its mirror image in the full domain -- so comparing the raw reduced average to a full-model average is a
    # category error (it is the source of the apparent "+x-z 15% transverse error", which is NOT a Gram bug:
    # the demag/energy matches the full model to ~5e-6).  Report the physical FULL-DOMAIN average as `M_avg`:
    # component c survives iff it is EVEN under every image plane a (normal-to-plane parity = -sign_a,
    # tangential parity = +sign_a; sign +1 symmetric/field-parallel, -1 antisymmetric/field-perpendicular),
    # otherwise its full-domain average is exactly 0.  The raw reduced average stays available as M_avg_reduced.
    M_avg = M_avg_reduced.copy()
    if image is not None:
        for _c in range(3):
            for _a, _s in _image.parse_image_string(image):
                if ((-_s) if _c == _a else _s) < 0:      # component c odd across plane a -> cancels in full domain
                    M_avg[_c] = 0.0
                    break
    out = _OperatorBackedResult(
               M=M_el, M_avg=M_avg, gfM=gfM, iters=int(iters), demag=D, ndof=n_face, n_el=n_el,
               n_charge=n_charge, nonlinear=bh_table is not None, linear_solver=solver_used,
               preconditioner=preconditioner, preconditioner_requested=preconditioner_requested,
               preconditioner_policy=preconditioner_policy,
               order=int(order), curve_order=curve_order, image=image,
               gram_backend=gram_backend,
               exact_dense_normalized_gram=bool(
                   getattr(H, "uses_exact_dense_normalized_gram", False)),
               operator_reused=bool(operator_reused),
               symmetry_constrained_dofs=len(symmetry_constrained_dofs), setup_wall_s=setup_wall_s,
               periodic_slave_dofs=len(periodic_slave_dofs),
               cyclic_periodic_boundaries=cyclic_periodic_boundaries,
               solve_wall_s=solve_wall_s, post_wall_s=time.perf_counter() - t_post,
               total_wall_s_internal=time.perf_counter() - t_total,
               fes_wall_s=t_before_charge_gram - t_before_fes,
               charge_gram_wall_s=t_after_charge_gram - t_before_charge_gram,
               charge_basis_wall_s=charge_build_timings.get("charge_basis_wall_s"),
               charge_gram_cpp_wall_s=charge_build_timings.get("charge_gram_cpp_wall_s"),
               hex_state_check_wall_s=charge_build_timings.get("hex_state_check_wall_s"),
               projection_wall_s=t_after_projection - t_after_charge_gram,
               demag_probe_wall_s=t_after_demag_probe - t_after_projection,
               operator_cache=operator_cache)
    out["_charge_gram"] = H
    out["_m_coefficients"] = np.ascontiguousarray(m, dtype=np.float64)
    out["_magnetization_sources"] = tuple(magnetization_sources)
    out["magnetization_source_count"] = len(magnetization_sources)
    for _k, _v in charge_build_timings.items():
        out.setdefault(_k, _v)
    if cpp_solve_timings:
        out["cpp_solve_timings"] = cpp_solve_timings
        for _k, _v in cpp_solve_timings.items():
            out.setdefault(_k, _v)
    if _LAST_NONLINEAR_SOLVE_STATS:
        out["nonlinear_solve_stats"] = dict(_LAST_NONLINEAR_SOLVE_STATS)
        for _k, _v in _LAST_NONLINEAR_SOLVE_STATS.items():
            out.setdefault(_k, _v)
    if image is not None:
        out["M_avg_reduced"] = M_avg_reduced
    if hmat_stats is not None:
        out["hmat_stats"] = hmat_stats
    if int(order) in (1, 2):
        # rad.Fld is part of the solved-object contract.  Materialize its
        # immutable C++ charge source and source tree now, while gfM/mesh are
        # hot, instead of charging the first observation request and repeating
        # Python source packing on every subsequent request.
        from ._field_batch import _materialize_field_evaluator
        _materialize_field_evaluator(out)
    out["post_wall_s"] = time.perf_counter() - t_post
    out["total_wall_s_internal"] = time.perf_counter() - t_total
    return out


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
    return a.mat


def _solve_nonlinear_picard_mass_riesz_cpp(mesh, fes, bh_table, H, n_face, h_ext, cg_tol, cg_maxit,
                                           nl_maxit, nl_tol, require_convergence=True):
    """Secant-reluctivity Picard with the same C++ mass-Riesz W-CG as the linear HDiv path.

    This is the production nonlinear path inherited from the tet/mass-Riesz line: freeze the elementwise
    secant reluctivity nu_sec(|M|)=|H(M)|/|M|, solve the SPD Galerkin problem

        ( INT nu_sec(M_old) u.v dx + B^T G B ) m_new = M_mass h_ext

    by `solve_configured_linear_material_mass_riesz`, then refresh nu_sec from the inverse BH table.  It keeps the
    weak HDiv / NGSolve formulation (no scalar global chi shortcut) but avoids the energy-Newton tensor
    tangent and line search.  `energy-newton` remains available for hard saturation / robustness checks."""
    rhs_src = _geometry_mass_apply(H, h_ext)
    uf, vf = fes.TnT()
    gfM = ng.GridFunction(fes); l2 = ng.L2(mesh, order=0); gfNu = ng.GridFunction(l2)

    if isinstance(bh_table, dict):
        mats = list(mesh.GetMaterials())
        missing = sorted(set(mats) - set(bh_table))
        if missing:
            raise ValueError("vim.Solve: bh_table dict missing region(s) %s; mesh materials are %s"
                             % (missing, mats))
        region_names = list(bh_table)
        name_to_ridx = {nm: k for k, nm in enumerate(region_names)}
        region_fields = []
        elem_region = np.array([name_to_ridx[mesh[ng.ElementId(ng.VOL, i)].mat] for i in range(mesh.ne)],
                               dtype=int)
        invchi0_e = np.empty(mesh.ne)
        for ridx, nm in enumerate(region_names):
            arr = np.asarray(bh_table[nm], float)
            if arr.ndim != 2 or arr.shape[1] != 2:
                raise ValueError("vim.Solve: bh_table[%r] must be [[H,B], ...] (A/m, T)" % (nm,))
            f, _, _ = _bh_inverse_funcs(arr[:, 0], arr[:, 1])
            region_fields.append(f)
            _, nd0 = f(np.array([1e-12]))
            invchi0_e[elem_region == ridx] = max(float(nd0[0]), 1e-30)

        def _nu_sec_all(Mmag):
            out = np.empty_like(Mmag)
            for ridx, f in enumerate(region_fields):
                sel = elem_region == ridx
                if np.any(sel):
                    out[sel] = f(Mmag[sel])[0]
            return out
    else:
        arr = np.asarray(bh_table, float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("vim.Solve: bh_table must be [[H,B], ...] (A/m, T)")
        fields, _, _ = _bh_inverse_funcs(arr[:, 0], arr[:, 1])
        _, nd0 = fields(np.array([1e-12]))
        invchi0_e = np.full(mesh.ne, max(float(nd0[0]), 1e-30))

        def _nu_sec_all(Mmag):
            return fields(Mmag)[0]

    def _Mmag(m):
        gfM.vec.FV().NumPy()[:] = m
        gfn = ng.GridFunction(l2)
        gfn.Set(ng.sqrt(ng.InnerProduct(gfM, gfM) + 1e-30))
        return np.maximum(gfn.vec.FV().NumPy(), 1e-30)

    def _W_matrix(nu_vals):
        gfNu.vec.FV().NumPy()[:] = np.maximum(np.asarray(nu_vals, float), 1e-30)
        a = ng.BilinearForm(fes)
        a += gfNu * uf * vf * ng.dx
        a.Assemble()
        return a.mat

    def _solve_W(W_matrix, x0=None):
        res = _h_solve_mass_riesz(
            H, W_matrix, int(n_face), 1.0, rhs_src, cg_tol, int(cg_maxit), x0=x0)
        _capture_cpp_solve_timings(res)
        it = int(res["iters"])
        if it >= int(cg_maxit):
            raise RuntimeError("vim.Solve (Picard mass-Riesz inner W-CG): did NOT converge in %d iters "
                               "(n_face=%d); tighten gram_eps or raise maxit." % (cg_maxit, n_face))
        return np.asarray(res["m"], float), it

    nu = np.maximum(invchi0_e, 1e-30)
    m = np.zeros(n_face, dtype=float)
    rel_step = float("inf")
    nit = 0
    stats = {
        "nonlinear_solver": "picard-mass-riesz",
        "nonlinear_picard_iters": 0,
        "nonlinear_material_relaxation": 0.7,
        "nonlinear_linear_inner_iters": 0,
    }
    # Damping the material coefficient rather than the solved field keeps every outer step an SPD Galerkin
    # solve while avoiding saturation ping-pong on steep tables.
    relax = 0.7
    for it in range(int(nl_maxit)):
        nit = it + 1
        m_new, inner_iterations = _solve_W(
            _W_matrix(nu), x0=m if it > 0 else None
        )
        stats["nonlinear_picard_iters"] += 1
        stats["nonlinear_linear_inner_iters"] += int(inner_iterations)
        rel_step = float(np.linalg.norm(m_new - m)) / (float(np.linalg.norm(m_new)) + 1e-30)
        m = m_new
        nu_new = np.maximum(_nu_sec_all(_Mmag(m)), 1e-30)
        if rel_step < nl_tol:
            stats["nonlinear_final_rel_step"] = float(rel_step)
            stats["nonlinear_converged_final_stage"] = True
            _capture_nonlinear_solve_stats(stats)
            return m, nit
        nu = relax * nu_new + (1.0 - relax) * nu
    if not require_convergence:
        stats["nonlinear_final_rel_step"] = float(rel_step)
        stats["nonlinear_converged_final_stage"] = False
        _capture_nonlinear_solve_stats(stats)
        return m, nit
    stats["nonlinear_final_rel_step"] = float(rel_step)
    stats["nonlinear_converged_final_stage"] = False
    _capture_nonlinear_solve_stats(stats)
    raise RuntimeError("vim.Solve (Picard mass-Riesz): did NOT converge -- rel step=%.2e > "
                       "nl_tol=%.1e after %d iters.  Try nonlinear_solver='energy-newton' for the robust "
                       "co-energy Newton path." % (rel_step, nl_tol, nit))


def _solve_nonlinear_energy_cpp(mesh, fes, bh_table, H, n_face, h_ext, cg_tol, cg_maxit,
                                nl_maxit, nl_tol, m0=None, *, inner_tol="auto", warmstart="linear",
                                continuation_steps=1, reuse_tangent_steps=1, cg_x0=False,
                                inner_preconditioner="mass-riesz"):
    """SYMMETRIC ENERGY-NEWTON with an all-C++ inner solve -- the production nonlinear soft-iron path (brings
    the nonlinear solve to C++ parity with the linear mass-Riesz CG; the default for iron-only nonlinear,
    replacing the former Python nonlinear Newton path).

    Co-energy / reluctivity form -- SYMMETRIC + M_mass^-1-free.  The co-energy functional
      E(m) = INT W_co(|M|) dx + 1/2 m.(N m) - (M_mass h_ext).m   (M = the flux field of m; W_co = INT H dM)
    is CONVEX; its gradient is the residual
      R(m) = INT H(M).v dx + N m - M_mass h_ext    (H(M) = the INVERSE-BH reluctance field, H = nu_sec M)
    and its Hessian (the Newton Jacobian) is SYMMETRIC
      J = W_tan + N,  W_tan = INT nu_d u.v dx  (nu_d = dH/dM = (dM/dH)^-1 differential reluctivity tensor),
    so each Newton step  (W_tan + N) dm = -R  is solved by the EXISTING C++ symmetric W-CG
    (`solve_configured_linear_material_mass_riesz`: W = W_tan as both the system mass AND the mass-Riesz PARDISO
    preconditioner; N via the symmetric charge-Gram H-matvec).  For large exploratory scaling runs,
    `inner_preconditioner="jacobi"` switches only the inner W-CG preconditioner to the exact diagonal of
    (W+N), avoiding per-CG PARDISO phase-33 solves.  No SciPy linear solve and no M_mass^-1 are used.

    Globalization: a chi0 (zero-field) LINEAR W-CG warmstart; an Armijo line search on the CONVEX ENERGY E
    (the merit -- ||R|| stalls in saturation where the inverse-BH
    H(M) blows up, but E keeps decreasing); a HARD-SATURATION BARRIER in the inverse BH (M cannot exceed
    Msat) that repels the M-iterates from the unphysical |M| > Mmax region.  Convergence: tight (relative
    Newton step < nl_tol -- 1-2 iters at moderate drive, == the forward Newton to ~1e-13), OR -- for the deep-
    saturation regime where the hard-saturation M-form intrinsically limit-cycles at the achievable precision
    -- a settled-step acceptance (rel step < 3e-4 for 5 consecutive iters -> accept the best-energy iterate;
    the production order-1 limit-cycle plateau is ~1.5-1.9e-4, so the floor is 3e-4; M matched the
    analytic uniform sphere to ~2e-3 at H0 up to 5e6, knee*5000).  Single-region (scalar
    bh_table) AND per-region (dict) iron; prescribed fixed-M sources enter through h_ext before this
    constitutive solve.  CALLER opens TaskManager."""
    rhs_src = _geometry_mass_apply(H, h_ext)
    Id = ng.Id(3); uf, vf = fes.TnT()
    gfM = ng.GridFunction(fes); l2 = ng.L2(mesh, order=0)
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
        return np.asarray(H.apply_configured_demag(_f64(v), True), float)

    def _Mmag(m):
        gfM.vec.FV().NumPy()[:] = m
        gfn = ng.GridFunction(l2); gfn.Set(ng.sqrt(ng.InnerProduct(gfM, gfM) + 1e-30))
        return np.maximum(gfn.vec.FV().NumPy(), 1e-30)

    def _bH(H_cf):
        lf = ng.LinearForm(fes); lf += H_cf * vf * ng.dx; lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def _W_matrix(weight_cf, tensor):
        a = ng.BilinearForm(fes)
        a += (ng.InnerProduct(weight_cf * uf, vf) if tensor else weight_cf * uf * vf) * ng.dx
        a.Assemble()
        return a.mat

    if inner_preconditioner not in ("mass-riesz", "jacobi"):
        raise ValueError("vim.Solve: energy-newton inner_preconditioner must be 'mass-riesz' or 'jacobi' "
                         "(got %r)" % (inner_preconditioner,))

    def _solve_W(W_matrix, rhs, *, tol_override=None, x0=None):
        solve_tol = float(cg_tol if tol_override is None else max(float(cg_tol), float(tol_override)))
        if inner_preconditioner == "jacobi":
            res = _h_solve_auto_prec(
                H, W_matrix, int(n_face), 1.0, rhs, solve_tol, int(cg_maxit), x0=x0)
        else:
            res = _h_solve_mass_riesz(
                H, W_matrix, int(n_face), 1.0, rhs, solve_tol, int(cg_maxit), x0=x0)
        _capture_cpp_solve_timings(res)
        it = int(res["iters"])
        if it >= int(cg_maxit):
            raise RuntimeError("vim.Solve (energy-Newton inner W-CG): did NOT converge in %d iters "
                               "(n_face=%d); the (W_tan + N) operator is SPD, so this means an ill-"
                               "conditioned tangent/mesh -- tighten gram_eps or raise maxit." % (cg_maxit, n_face))
        return np.asarray(res["m"], float), it

    def _energy(m, rhs):                        # E(m) = INT W_co(|M|) dx + 1/2 m.Nm - rhs.m
        return float(np.dot(_wco_all(_Mmag(m)), Vol)) + 0.5 * float(m @ _N_apply(m)) - float(rhs @ m)

    def _forcing_tol(prev_rel_step, stage_final):
        if inner_tol in (None, "fixed"):
            return float(cg_tol)
        if inner_tol == "auto":
            if not stage_final:
                return max(float(cg_tol), 1e-3)
            if not np.isfinite(prev_rel_step):
                return max(float(cg_tol), 1e-3)
            return max(float(cg_tol), min(1e-3, 0.25 * max(float(prev_rel_step), float(nl_tol))))
        return max(float(cg_tol), float(inner_tol))

    nstage = max(1, int(continuation_steps))
    reuse_steps = max(1, int(reuse_tangent_steps))
    if warmstart not in ("linear", "picard", "none"):
        raise ValueError("vim.Solve: warmstart must be 'linear', 'picard', or 'none' (got %r)" % (warmstart,))
    stats = {
        "nonlinear_inner_preconditioner": inner_preconditioner,
        "nonlinear_inner_tol": inner_tol,
        "nonlinear_continuation_steps": int(nstage),
        "nonlinear_reuse_tangent_steps": int(reuse_steps),
        "nonlinear_cg_x0": bool(cg_x0),
        "nonlinear_newton_iters": 0,
        "nonlinear_warmstart_solves": 0,
        "nonlinear_linear_inner_iters": 0,
        "nonlinear_line_search_backtracks": 0,
        "nonlinear_tangent_assemblies": 0,
        "nonlinear_tangent_reuses": 0,
        "nonlinear_fresh_tangent_retries": 0,
    }
    alphas = np.linspace(1.0 / nstage, 1.0, nstage)
    invchi0 = None
    m0_provided = m0 is not None
    m = np.asarray(m0, float).copy() if m0_provided else np.zeros(n_face, dtype=float)
    if m0_provided and nstage > 1:
        m *= float(alphas[0])
    dm_prev = None
    total_nit = 0
    final_rel_step = float("inf")
    final_settled = 0
    converged_final = False

    for istage, alpha in enumerate(alphas):
        rhs_stage = float(alpha) * rhs_src
        stage_final = istage == len(alphas) - 1
        if not m0_provided and istage == 0 and warmstart == "linear":
            # chi0 (zero-field) LINEAR warmstart: (M_{1/chi0} + N) m = M_mass h_ext.
            invchi0 = ng.GridFunction(l2)
            invchi0.vec.FV().NumPy()[:] = 1.0 / np.maximum(chi0_e, 1.0)
            m, it0 = _solve_W(
                _W_matrix(invchi0, tensor=False), rhs_stage, tol_override=max(cg_tol, 1e-6))
            stats["nonlinear_warmstart_solves"] += 1
            stats["nonlinear_linear_inner_iters"] += int(it0)

        converged = False
        settled = 0
        rel_step = float("inf")
        E = _energy(m, rhs_stage)
        Ebest = E
        mbest = m.copy()
        cached_W = None
        cached_it = -10**9
        for it in range(int(nl_maxit)):
            total_nit += 1
            stats["nonlinear_newton_iters"] += 1
            gfM.vec.FV().NumPy()[:] = m
            H_cf, nud = _reluct(gfM)
            R = _bH(H_cf) + _N_apply(m) - rhs_stage
            rebuild = cached_W is None or (it - cached_it) >= reuse_steps
            if rebuild:
                cached_W = _W_matrix(nud, tensor=True)
                cached_it = it
                stats["nonlinear_tangent_assemblies"] += 1
            else:
                stats["nonlinear_tangent_reuses"] += 1
            solve_tol = _forcing_tol(rel_step, stage_final)
            dm, itlin = _solve_W(cached_W, -R, tol_override=solve_tol,
                                 x0=dm_prev if (cg_x0 and dm_prev is not None) else None)
            stats["nonlinear_linear_inner_iters"] += int(itlin)
            dec = float(-dm @ R)                             # dm.(-R) = dm^T J dm >= 0
            lam = 1.0
            E0 = E
            backtracks = 0
            while lam > 1e-10:
                if _energy(m + lam * dm, rhs_stage) <= E0 - 1e-4 * lam * dec:
                    break
                lam *= 0.5
                backtracks += 1
            if lam <= 1e-10 and not rebuild:
                # A reused/chord tangent can occasionally lose descent.  Pay for one fresh tangent before
                # accepting a microscopic step.
                cached_W = _W_matrix(nud, tensor=True)
                cached_it = it
                stats["nonlinear_tangent_assemblies"] += 1
                stats["nonlinear_fresh_tangent_retries"] += 1
                dm, itlin = _solve_W(cached_W, -R, tol_override=max(cg_tol, min(1e-4, solve_tol)),
                                     x0=dm_prev if (cg_x0 and dm_prev is not None) else None)
                stats["nonlinear_linear_inner_iters"] += int(itlin)
                dec = float(-dm @ R)
                lam = 1.0
                backtracks = 0
                while lam > 1e-10:
                    if _energy(m + lam * dm, rhs_stage) <= E0 - 1e-4 * lam * dec:
                        break
                    lam *= 0.5
                    backtracks += 1
            stats["nonlinear_line_search_backtracks"] += int(backtracks)
            step = lam * dm
            rel_step = float(np.linalg.norm(step)) / (float(np.linalg.norm(m)) + 1e-30)
            m = m + step
            dm_prev = dm
            E = _energy(m, rhs_stage)
            if E < Ebest:
                Ebest = E
                mbest = m.copy()
            settled = settled + 1 if rel_step < 3e-4 else 0
            if rel_step < nl_tol:
                converged = True
                break
            if settled >= 5:
                converged = True
                m = mbest
                break
        if not converged:
            m = mbest
            _capture_nonlinear_solve_stats(stats)
            raise RuntimeError("vim.Solve (energy-Newton): did NOT converge -- rel step=%.2e (tol %.1e), "
                               "%d settled iters after %d (returning M would be a silent wrong result).  For an "
                               "extreme-saturation / ill-conditioned case, increase continuation or use the "
                               "mass-Riesz preconditioner." % (rel_step, nl_tol, settled, total_nit))
        final_rel_step = rel_step
        final_settled = settled
        converged_final = stage_final

    stats["nonlinear_final_rel_step"] = float(final_rel_step)
    stats["nonlinear_final_settled_iters"] = int(final_settled)
    stats["nonlinear_converged_final_stage"] = bool(converged_final)
    _capture_nonlinear_solve_stats(stats)
    return m, total_nit
