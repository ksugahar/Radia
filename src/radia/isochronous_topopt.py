"""Per-element density adjoint for HDiv-MMM topology optimization.

Application-layer module for isochronous-magnet (and general accelerator
pole/yoke) density topology optimization on the HDiv-MMM forward engine
(design record: ``docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md``).
It composes :class:`radia.vim.DemagOperator` (build-once geometry operator
``N = B^T G B``) with the per-element design variable ``s_e = 1/chi_e``
(an ``L2(order=0)`` mass weight) into the adjoint gradient route:

    (M_s + N) m      = f_state        state solve (native H-matrix Jacobi-PCG)
    (M_s + N) lambda = f_adjoint      adjoint (same operator, same factor)
    J(s)             = f_adjoint^T m  (the design-dependent objective part)
    dJ/ds_e          = -lambda^T M^(e) m
                     = -Integrate(lambda . m, element_wise=True)[e]

The last identity is exact because ``A(s) = sum_e s_e M^(e) + N`` is affine in
``s`` and the element mass contraction of two HDiv fields IS the element-wise
integral of their inner product -- all sensitivities come from ONE
``Integrate`` call, with no per-mode machinery and no dense Jacobian.

The field functional ``F^T c`` (weighted point evaluations of the
demagnetizing field, e.g. orbit-quadrature ``dB_z/dx`` for the focusing
integral ``gL``) is realized by kernel reciprocity: the adjoint load is the
mass projection of the analytic H field of a point-DIPOLE array at the
evaluation points.  The dense magnetization-to-field matrix ``F`` is never
formed; derivative functionals enter as +/- point pairs
(:func:`gradient_pair_points`).

Verified 2026-07-28 on the promoted unit-ball gate (270 tets, log-uniform
``s`` in [1e-2, 1]): adjoint gradient matches
central finite differences to 8.1e-10 (directional, all elements) and 3.9e-7
worst-case per element (the smallest-gradient element, FD noise floor); the
reciprocity load matches the independent C++ analytic charge evaluator
(``vim.FieldFromSolution``) to 1.1e-10 at ``bonus_intorder=10``.  Locked by
``tests/test_isochronous_topopt.py``.

The caller wraps all NGSolve work (construction and solves) in
``with TaskManager():`` per the repository TaskManager policy.  The LP/OC
design-update drivers live in :mod:`radia.topology_optimization`; this module
supplies the objective/gradient callback layer above the same operators.
"""
import time
from dataclasses import dataclass
from math import pi

import numpy as np
import ngsolve as ng
from ngsolve.krylovspace import CGSolver

from .vim import DemagOperator, FieldFromSolution

MU0 = 4.0e-7 * pi
#: Validated ersatz-void susceptibility floor: with the s-weighted mass-Riesz
#: preconditioner, CG iteration counts DECREASE as chi_void -> 0 (measured
#: 74 -> 46 from 1e-1 to 1e-6), and the void response stays exactly the
#: physical chi*H over five decades.  See the design document, Sec. 3.
CHI_MIN = 1.0e-6

__all__ = [
    "MU0", "CHI_MIN", "AdjointGradientResult", "DensityAdjointVIM",
    "FunctionalLinearization", "DensityDesignResult", "HelmholtzFilter",
    "HeavisideProjection", "density_discreteness",
    "iron_only_verification_ready",
    "IronOnlyVerification", "density_to_s",
    "density_gradient_from_s_gradient", "gradient_pair_points",
    "dipole_array_field_cf", "field_functional_load", "uniform_field_load",
    "demag_field_from_solution", "orbit_arc_points", "optimize_density",
    "iron_only_mesh", "verify_design_iron_only",
]


# --------------------------------------------------------------------------
# density <-> s mapping
# --------------------------------------------------------------------------
def density_to_s(density, chi_iron, chi_min=CHI_MIN, penalty=1.0):
    """Map a per-element density ``rho in [0, 1]`` to ``s = 1/chi(rho)``.

    SIMP-style interpolation in susceptibility, ``chi(rho) = chi_min +
    (chi_iron - chi_min) rho^penalty``.  ``penalty=1`` (default) is the
    linear map; ``penalty=3`` makes intermediate densities material-
    inefficient and drives designs toward 0/1 (use for runs whose result
    will be thresholded/manufactured).  The filtering layer of the design
    loop operates on ``rho`` before this map.  ``chi_min`` defaults to the
    validated ersatz floor :data:`CHI_MIN`.
    """
    rho = np.asarray(density, dtype=float)
    chi_iron = float(chi_iron)
    chi_min = float(chi_min)
    penalty = float(penalty)
    if not (np.isfinite(chi_iron) and np.isfinite(chi_min)
            and chi_iron > chi_min > 0.0):
        raise ValueError(
            "density_to_s: need chi_iron > chi_min > 0 (got chi_iron=%r, chi_min=%r)"
            % (chi_iron, chi_min))
    if not np.isfinite(penalty) or penalty < 1.0:
        raise ValueError("density_to_s: penalty must be >= 1")
    if not np.all(np.isfinite(rho)):
        raise ValueError("density_to_s: density must be finite")
    if np.any(rho < -1e-9) or np.any(rho > 1.0 + 1e-9):
        raise ValueError(
            "density_to_s: density must lie in [0, 1] (got min=%r, max=%r)"
            % (float(rho.min()), float(rho.max())))
    chi = chi_min + (chi_iron - chi_min) * np.clip(rho, 0.0, 1.0) ** penalty
    return 1.0 / chi


def density_gradient_from_s_gradient(density, s_gradient, chi_iron,
                                     chi_min=CHI_MIN, penalty=1.0):
    """Chain rule ``dJ/drho_e = dJ/ds_e * ds/drho_e`` for :func:`density_to_s`.

    ``ds/drho = -(chi_iron - chi_min) penalty rho^(penalty-1) / chi(rho)^2``.
    """
    rho = np.asarray(density, dtype=float)
    grad_s = np.asarray(s_gradient, dtype=float)
    if rho.shape != grad_s.shape:
        raise ValueError("density_gradient_from_s_gradient: shape mismatch %r vs %r"
                         % (rho.shape, grad_s.shape))
    chi_iron = float(chi_iron)
    chi_min = float(chi_min)
    penalty = float(penalty)
    if not (np.isfinite(chi_iron) and np.isfinite(chi_min)
            and chi_iron > chi_min > 0.0):
        raise ValueError("density_gradient_from_s_gradient: need "
                         "chi_iron > chi_min > 0")
    if not np.isfinite(penalty) or penalty < 1.0:
        raise ValueError("density_gradient_from_s_gradient: penalty must be >= 1")
    if (not np.all(np.isfinite(rho)) or not np.all(np.isfinite(grad_s))
            or np.any(rho < 0.0) or np.any(rho > 1.0)):
        raise ValueError("density_gradient_from_s_gradient: density and "
                         "gradient must be finite, with density in [0, 1]")
    chi = chi_min + (chi_iron - chi_min) * rho ** penalty
    dchi = ((chi_iron - chi_min) * penalty * rho ** (penalty - 1.0)
            if penalty != 1.0 else (chi_iron - chi_min) * np.ones_like(rho))
    return grad_s * (-dchi / (chi * chi))


# --------------------------------------------------------------------------
# field-functional loads (F^T c by kernel reciprocity)
# --------------------------------------------------------------------------
def gradient_pair_points(points, weights, delta, axis=0, direction=None):
    """Realize ``sum_i w_i d(field)/d(direction) (x_i)`` as +/- evaluation pairs.

    Returns ``(points_out, weights_out)`` with ``2 N`` rows: each input point
    splits into ``x_i +/- (delta/2) d_i`` carrying weights ``+/- w_i/delta``
    (central finite difference of the point functional; ``delta`` is a
    physical stencil length on the orbit, not a numerical epsilon).
    ``direction`` is a single (3,) vector or per-point (N, 3) vectors
    (normalized here); when omitted the coordinate ``axis`` is used --
    per-point radial directions realize the ``d B_z/d r`` orbit functionals.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    wts = np.asarray(weights, dtype=float).ravel()
    if len(pts) != len(wts):
        raise ValueError("gradient_pair_points: %d points but %d weights"
                         % (len(pts), len(wts)))
    if not delta > 0.0:
        raise ValueError("gradient_pair_points: delta must be positive")
    if direction is None:
        if axis not in (0, 1, 2):
            raise ValueError("gradient_pair_points: axis must be 0, 1, or 2")
        direction = np.zeros(3)
        direction[axis] = 1.0
    d = np.asarray(direction, dtype=float)
    d = np.broadcast_to(d.reshape(-1, 3), pts.shape).copy()
    norms = np.linalg.norm(d, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError("gradient_pair_points: zero direction vector")
    d /= norms[:, None]
    shift = 0.5 * float(delta) * d
    return (np.concatenate([pts + shift, pts - shift]),
            np.concatenate([wts / delta, -wts / delta]))


def orbit_arc_points(radius, z, n, span=(0.0, 2.0 * pi)):
    """Equispaced points on a circular arc (radius, height ``z``), with the
    per-point radial unit vectors for :func:`gradient_pair_points`.

    Returns ``(points (n,3), radial_directions (n,3))``.  A full circle
    (span of 2 pi) omits the duplicate endpoint.
    """
    lo, hi = float(span[0]), float(span[1])
    full = abs((hi - lo) - 2.0 * pi) < 1e-12
    theta = np.linspace(lo, hi, int(n), endpoint=not full)
    pts = np.stack([radius * np.cos(theta), radius * np.sin(theta),
                    float(z) * np.ones_like(theta)], axis=1)
    radial = np.stack([np.cos(theta), np.sin(theta),
                       np.zeros_like(theta)], axis=1)
    return pts, radial


def dipole_array_field_cf(points, moments):
    """Analytic H field CoefficientFunction of a point-dipole array.

    ``H(y) = (1/4pi) sum_i [3 (m_i . r) r/|r|^5 - m_i/|r|^3]``, ``r = y - x_i``.
    Singular at the dipole points -- integrate it only over regions that
    exclude them (:func:`field_functional_load` enforces this).
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    mom = np.asarray(moments, dtype=float).reshape(-1, 3)
    if len(pts) != len(mom):
        raise ValueError("dipole_array_field_cf: %d points but %d moments"
                         % (len(pts), len(mom)))
    if not len(pts):
        raise ValueError("dipole_array_field_cf: empty dipole array")
    total = None
    for (px, py, pz), (mx, my, mz) in zip(pts, mom):
        rx, ry, rz = ng.x - px, ng.y - py, ng.z - pz
        r2 = rx * rx + ry * ry + rz * rz
        r1 = ng.sqrt(r2)
        r3 = r2 * r1
        r5 = r2 * r3
        mdotr = mx * rx + my * ry + mz * rz
        term = ng.CoefficientFunction((3.0 * mdotr * rx / r5 - mx / r3,
                                       3.0 * mdotr * ry / r5 - my / r3,
                                       3.0 * mdotr * rz / r5 - mz / r3))
        total = term if total is None else total + term
    return total * (1.0 / (4.0 * pi))


def _require_points_outside(mesh, pts, who):
    """Fail loud when an evaluation point lies inside the meshed iron."""
    mapped = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
    if not (isinstance(mapped, np.ndarray) and mapped.dtype.names
            and "nr" in mapped.dtype.names):
        raise RuntimeError(
            "%s: NGSolve's vectorized MeshPoint array API is required for the "
            "inside-mesh check" % who)
    inside = np.flatnonzero(np.asarray(mapped["nr"]) >= 0)
    if len(inside):
        raise ValueError(
            "%s: %d evaluation/dipole point(s) lie INSIDE the meshed body "
            "(first: %s) -- the dipole reciprocity load is singular there; "
            "orbit points must be outside the iron mesh"
            % (who, len(inside), pts[inside[0]].tolist()))


def field_functional_load(fes, points, weights, axis=2, scale=MU0,
                          bonus_intorder=10):
    """Assembled adjoint load ``f`` with ``f^T m = scale * sum_i w_i H_d[m]_axis(x_i)``.

    Kernel reciprocity: the load is the mass projection of the analytic field
    of dipoles ``w_i e_axis`` at ``x_i``.  With ``scale=MU0`` (default),
    ``f^T m`` is the magnetization-dependent part of
    ``sum_i w_i B_axis(x_i)`` in Tesla; the design-independent coil part
    ``scale * sum_i w_i H_ext_axis(x_i)`` is the caller's constant.
    ``bonus_intorder=10`` reproduces the independent C++ charge evaluator to
    ~1e-10 relative on the verification geometry (points >= half an element
    size outside the body); the returned value is exact FD-consistent for the
    adjoint regardless of quadrature, since J is DEFINED as ``f^T m``.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    wts = np.asarray(weights, dtype=float).ravel()
    if len(pts) != len(wts):
        raise ValueError("field_functional_load: %d points but %d weights"
                         % (len(pts), len(wts)))
    if axis not in (0, 1, 2):
        raise ValueError("field_functional_load: axis must be 0, 1, or 2")
    _require_points_outside(fes.mesh, pts, "field_functional_load")
    moments = np.zeros((len(wts), 3))
    moments[:, axis] = wts
    cf = dipole_array_field_cf(pts, moments)
    v = fes.TestFunction()
    f = ng.LinearForm(fes)
    f += float(scale) * ng.InnerProduct(cf, v) * ng.dx(
        bonus_intorder=int(bonus_intorder))
    f.Assemble()
    return f


def uniform_field_load(fes, H_ext):
    """Assembled state load for a uniform external field ``H_ext`` (A/m)."""
    H = tuple(float(c) for c in H_ext)
    if len(H) != 3:
        raise ValueError("uniform_field_load: H_ext must have 3 components")
    v = fes.TestFunction()
    f = ng.LinearForm(fes)
    f += ng.CoefficientFunction(H) * v * ng.dx
    f.Assemble()
    return f


def demag_field_from_solution(demag, gfM, points, algorithm="direct"):
    """Demagnetizing field H_d (A/m) of a solved magnetization at ``points``.

    Independent verification/reporting path: evaluates the exact analytic
    field of the BDM charge representation through the C++ evaluator of
    ``vim.FieldFromSolution``.  This is the ONE sanctioned coupling point to
    the vim result-dict internals (the same keys ``vim.Solve`` stores).
    """
    if not isinstance(demag, DemagOperator):
        raise TypeError("demag_field_from_solution: demag must be a vim.DemagOperator")
    if gfM.space is not demag.space:
        raise ValueError("demag_field_from_solution: gfM does not live on the "
                         "operator's HDiv space")
    res = {"gfM": gfM, "order": int(demag.space.globalorder),
           "_charge_gram": demag._G,
           "_m_coefficients": np.ascontiguousarray(
               gfM.vec.FV().NumPy(), dtype=np.float64)}
    return FieldFromSolution(res, points, algorithm=algorithm)


# --------------------------------------------------------------------------
# density filter
# --------------------------------------------------------------------------
class HelmholtzFilter:
    """Element-density Helmholtz filter ``rho_f = D_V^-1 B^T K^-1 B rho``.

    ``K = (u v + radius^2 grad u . grad v)`` on H1 order 1 (SPD, factored
    once); ``B_je = int_e phi_j``; ``D_V`` = element volumes.  ``radius`` is
    the minimum-feature length (a manufacturability input).  The K^-1 core is
    self-adjoint, so the gradient chain is the transpose map
    ``dJ/drho = B^T K^-1 B (D_V^-1 dJ/drho_f)`` -- :meth:`chain` -- exact to
    solver precision (locked by the filter FD test).

    The caller wraps construction and every apply in ``with TaskManager():``.
    """

    def __init__(self, mesh, radius):
        if not radius > 0.0:
            raise ValueError("HelmholtzFilter: radius must be positive")
        self.mesh = mesh
        self.radius = float(radius)
        fes = ng.H1(mesh, order=1)
        u, v = fes.TnT()
        a = ng.BilinearForm(fes, symmetric=True)
        a += (u * v + self.radius ** 2 * ng.grad(u) * ng.grad(v)) * ng.dx
        a.Assemble()
        self._fes = fes
        self._inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        self._volumes = np.asarray(
            ng.Integrate(ng.CoefficientFunction(1.0), mesh,
                         element_wise=True), dtype=float)
        self._l2 = ng.L2(mesh, order=0)
        if self._l2.ndof != int(mesh.ne):
            raise RuntimeError("HelmholtzFilter: L2(order=0) dof count %d != "
                               "element count %d" % (self._l2.ndof, mesh.ne))
        self._rho_gf = ng.GridFunction(self._l2)
        self._u_gf = ng.GridFunction(fes)

    def _bt_kinv_b(self, values):
        values = np.ascontiguousarray(values, dtype=float).ravel()
        if values.size != self._volumes.size:
            raise ValueError("HelmholtzFilter: vector has %d entries, mesh "
                             "has %d elements" % (values.size,
                                                  self._volumes.size))
        self._rho_gf.vec.FV().NumPy()[:] = values
        f = ng.LinearForm(self._fes)
        f += self._rho_gf * self._fes.TestFunction() * ng.dx
        f.Assemble()
        self._u_gf.vec.data = self._inv * f.vec
        return np.asarray(ng.Integrate(self._u_gf, self.mesh,
                                       element_wise=True), dtype=float)

    def apply(self, density):
        """Filtered element densities (same shape, smoothed over ``radius``)."""
        return self._bt_kinv_b(density) / self._volumes

    def chain(self, gradient_filtered):
        """Transpose map: gradient w.r.t. raw density from the filtered one."""
        return self._bt_kinv_b(
            np.asarray(gradient_filtered, dtype=float) / self._volumes)


class HeavisideProjection:
    """Smooth density projection with an exact analytic chain rule.

    ``beta`` controls sharpness and ``eta`` is the transition density.  Use a
    continuation sequence (for example 1, 2, 4, 8) across design phases; a
    large beta from a flat start suppresses useful topology gradients.
    """

    def __init__(self, beta, eta=0.5):
        self.beta = float(beta)
        self.eta = float(eta)
        if not np.isfinite(self.beta) or self.beta <= 0.0:
            raise ValueError("HeavisideProjection: beta must be positive")
        if not np.isfinite(self.eta) or not 0.0 < self.eta < 1.0:
            raise ValueError("HeavisideProjection: eta must lie in (0, 1)")
        self._denominator = (np.tanh(self.beta*self.eta)
                             + np.tanh(self.beta*(1.0-self.eta)))

    def apply(self, density):
        rho = np.asarray(density, dtype=float)
        if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
                or np.any(rho > 1.0)):
            raise ValueError("HeavisideProjection: density must lie in [0, 1]")
        return ((np.tanh(self.beta*self.eta)
                 + np.tanh(self.beta*(rho-self.eta))) / self._denominator)

    def derivative(self, density):
        rho = np.asarray(density, dtype=float)
        if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
                or np.any(rho > 1.0)):
            raise ValueError("HeavisideProjection: density must lie in [0, 1]")
        return (self.beta / self._denominator
                / np.cosh(self.beta*(rho-self.eta))**2)

    def chain(self, density, gradient_projected):
        gradient = np.asarray(gradient_projected, dtype=float)
        rho = np.asarray(density, dtype=float)
        if gradient.shape != rho.shape:
            raise ValueError("HeavisideProjection.chain: shape mismatch")
        return gradient*self.derivative(rho)


def density_discreteness(density, *, lower=0.1, upper=0.9):
    """Return stable grayness diagnostics for a continuous density design."""
    rho = np.asarray(density, dtype=float).reshape(-1)
    if (rho.size == 0 or not np.all(np.isfinite(rho))
            or np.any(rho < 0.0) or np.any(rho > 1.0)):
        raise ValueError("density_discreteness: density must be non-empty in [0,1]")
    lower, upper = float(lower), float(upper)
    if not 0.0 <= lower < upper <= 1.0:
        raise ValueError("density_discreteness: need 0 <= lower < upper <= 1")
    intermediate = (rho > lower) & (rho < upper)
    return {
        "intermediate_fraction": float(np.mean(intermediate)),
        "binary_distance_mean": float(np.mean(np.minimum(rho, 1.0-rho))),
        "binary_distance_max": float(np.max(np.minimum(rho, 1.0-rho))),
    }


def iron_only_verification_ready(density, *, maximum_intermediate_fraction=0.05,
                                 lower=0.1, upper=0.9):
    """Whether thresholded exact-void verification is scientifically useful."""
    limit = float(maximum_intermediate_fraction)
    if not np.isfinite(limit) or not 0.0 <= limit <= 1.0:
        raise ValueError("maximum_intermediate_fraction must lie in [0,1]")
    metrics = density_discreteness(density, lower=lower, upper=upper)
    return metrics["intermediate_fraction"] <= limit, metrics


# --------------------------------------------------------------------------
# adjoint solver
# --------------------------------------------------------------------------
@dataclass
class AdjointGradientResult:
    """State + adjoint solution pair with objective and per-element gradient."""
    objective: float          # f_adjoint^T m  (design-dependent part of J)
    gradient: np.ndarray      # dJ/ds_e, shape (n_el,)
    gfM: "ng.GridFunction"    # state magnetization m
    gfLambda: "ng.GridFunction"  # adjoint field lambda
    state_iterations: int
    adjoint_iterations: int


@dataclass
class FunctionalLinearization:
    """State + one adjoint per functional, all on one configured operator.

    ``values[k] = f_k^T m`` and ``jacobians[k, e] = d(f_k^T m)/ds_e`` for the
    functional loads passed to :meth:`DensityAdjointVIM.linearize` --
    functional 0 is conventionally the objective, the rest are constraints.
    """
    values: np.ndarray        # (n_loads,)
    jacobians: np.ndarray     # (n_loads, n_el)
    gfM: "ng.GridFunction"
    gfLambdas: tuple          # one adjoint GridFunction per load
    state_iterations: int
    adjoint_iterations: tuple


class DensityAdjointVIM:
    """Build-once HDiv-MMM operator with per-element density adjoint gradients.

    ``fes`` is the HDiv space on the design mesh (order 1 or 2).  The geometry
    operator ``N`` is built once (``demag=None``) or shared from an existing
    :class:`radia.vim.DemagOperator`; every design iterate then costs one
    weighted-mass assembly and the state/adjoint CG solves.  The caller
    wraps construction and every solve in ``with TaskManager():``.
    """

    def __init__(self, fes, demag=None, **gram_kwargs):
        if demag is None:
            demag = DemagOperator(fes, **gram_kwargs)
        else:
            if gram_kwargs:
                raise ValueError(
                    "DensityAdjointVIM: gram kwargs are only accepted when the "
                    "operator is built here (demag=None)")
            if demag.space is not fes:
                raise ValueError(
                    "DensityAdjointVIM: the supplied DemagOperator lives on a "
                    "different HDiv space")
        self.fes = fes
        self.mesh = fes.mesh
        self.demag = demag
        self.n_el = int(self.mesh.ne)
        self._l2 = ng.L2(self.mesh, order=0)
        if self._l2.ndof != self.n_el:
            raise RuntimeError(
                "DensityAdjointVIM: L2(order=0) dof count %d != element count %d"
                % (self._l2.ndof, self.n_el))
        self._s_gf = ng.GridFunction(self._l2)
        self._volumes = None

    @property
    def element_volumes(self):
        """Per-element volumes (the LP volume-budget row), shape (n_el,)."""
        if self._volumes is None:
            self._volumes = np.asarray(
                ng.Integrate(ng.CoefficientFunction(1.0), self.mesh,
                             element_wise=True), dtype=float)
        return self._volumes

    # -------------------------------------------------------------- internals
    def _system(self, s, *, need_ngsolve_solver=True):
        s = np.ascontiguousarray(s, dtype=float).ravel()
        if s.size != self.n_el:
            raise ValueError("DensityAdjointVIM: s has %d entries, mesh has %d "
                             "elements" % (s.size, self.n_el))
        if not np.all(s > 0.0):
            raise ValueError(
                "DensityAdjointVIM: s = 1/chi must be strictly positive "
                "(min=%r); map densities through density_to_s" % float(s.min()))
        self._s_gf.vec.FV().NumPy()[:] = s
        u, v = self.fes.TnT()
        mass = ng.BilinearForm(self.fes, symmetric=True)
        mass += self._s_gf * u * v * ng.dx
        mass.Assemble()
        if not need_ngsolve_solver:
            return mass, None, None
        A = mass.mat + self.demag.mat
        pre = mass.mat.Inverse(self.fes.FreeDofs(), inverse="sparsecholesky")
        return mass, A, pre

    def _cg(self, A, pre, rhs_vec, gf, tol, maxiter, warm_vec=None):
        # krylovspace's relative tol is anchored to the FIRST residual of the
        # run, which for a warm start is already tiny -- that would make warm
        # restarts HARDER, not cheaper.  Anchor the stopping criterion to the
        # x0=0 first residual ||rhs||_pre instead, so warm and cold solves
        # share one absolute criterion (cold behavior is unchanged).
        work = rhs_vec.CreateVector()
        work.data = pre * rhs_vec
        reference = abs(ng.InnerProduct(work, rhs_vec)) ** 0.5
        if reference == 0.0:
            gf.vec[:] = 0.0
            return 0
        inv = CGSolver(A, pre=pre, tol=None, atol=tol * reference,
                       maxiter=int(maxiter), printrates=False)
        if warm_vec is not None:
            gf.vec.data = warm_vec
            inv.Solve(rhs=rhs_vec, sol=gf.vec, initialize=False)
        else:
            inv.Solve(rhs=rhs_vec, sol=gf.vec, initialize=True)
        iters = int(inv.iterations)
        if iters >= int(maxiter):
            raise RuntimeError(
                "DensityAdjointVIM: CG did not converge within %d iterations "
                "(tol=%g); raise maxiter or inspect the design state"
                % (int(maxiter), tol))
        return iters

    def _native_solve_many(self, mass, rhs_vecs, tol, maxiter,
                           warm_vecs=None, mass_riesz=True):
        """Solve shared-matrix VIM systems in the configured C++ kernel.

        The weighted HDiv mass is registered once.  H-matrix application,
        Krylov updates, and true-residual stopping tests remain inside C++.
        ``mass_riesz=False`` uses the inexpensive exact system diagonal and is
        the study-scale default; ``True`` reuses the persistent PARDISO mass
        factor.  Python crosses the boundary only once per right-hand side.
        """
        gram = getattr(self.demag, "_G", None)
        solve_name = ("solve_configured_linear_material_mass_riesz"
                      if mass_riesz else
                      "solve_configured_linear_material_auto_prec")
        required = ("configure_mass_matrix_ngsolve", solve_name)
        if gram is None or any(not hasattr(gram, name) for name in required):
            raise RuntimeError(
                "DensityAdjointVIM: native H-matrix CG is unavailable "
                "on this DemagOperator; pass solver='ngsolve-cg' only for "
                "the explicit reference path")
        gram.configure_mass_matrix_ngsolve(mass.mat)
        rhs_vecs = list(rhs_vecs)
        if warm_vecs is None:
            warm_vecs = [None] * len(rhs_vecs)
        else:
            warm_vecs = list(warm_vecs)
            if len(warm_vecs) != len(rhs_vecs):
                raise ValueError(
                    "DensityAdjointVIM: native warm-start count mismatch")
        fields, iterations = [], []
        for rhs, warm_vec in zip(rhs_vecs, warm_vecs):
            rhs_array = np.ascontiguousarray(rhs.FV().NumPy(), dtype=float)
            x0 = None if warm_vec is None else np.ascontiguousarray(
                warm_vec.FV().NumPy(), dtype=float)
            if mass_riesz:
                result = gram.solve_configured_linear_material_mass_riesz(
                    1.0, rhs_array, float(tol), int(maxiter), True, x0=x0)
            else:
                result = gram.solve_configured_linear_material_auto_prec(
                    1.0, rhs_array, float(tol), int(maxiter), x0=x0)
            iters = int(result["iters"])
            if iters >= int(maxiter):
                raise RuntimeError(
                    "DensityAdjointVIM: native %s CG did not converge "
                    "within %d iterations (tol=%g)"
                    % ("mass-Riesz" if mass_riesz else "Jacobi",
                       int(maxiter), tol))
            gf = ng.GridFunction(self.fes)
            gf.vec.FV().NumPy()[:] = np.asarray(result["m"], dtype=float)
            fields.append(gf)
            iterations.append(iters)
        return fields, iterations

    # ------------------------------------------------------------ public API
    def solve(self, s, load, tol=1e-12, maxiter=20000, warm=None,
              solver="native-jacobi"):
        """One SPD solve ``(M_s + N) x = load.vec``.

        ``warm`` is an optional GridFunction whose vector seeds CG (warm start
        across design iterates).  Returns ``(GridFunction, iterations)``.
        """
        mass, A, pre = self._system(
            s, need_ngsolve_solver=(solver == "ngsolve-cg"))
        if solver in {"native", "native-jacobi"}:
            fields, iterations = self._native_solve_many(
                mass, [load.vec], tol, maxiter,
                warm_vecs=None if warm is None else [warm.vec],
                mass_riesz=(solver == "native"))
            return fields[0], iterations[0]
        if solver != "ngsolve-cg":
            raise ValueError(
                "DensityAdjointVIM.solve: solver must be 'native-jacobi', "
                "'native', or 'ngsolve-cg'")
        gf = ng.GridFunction(self.fes)
        iters = self._cg(A, pre, load.vec, gf, tol, maxiter,
                         warm_vec=None if warm is None else warm.vec)
        return gf, iters

    def linearize(self, s, state_load, functional_loads,
                  tol=1e-12, maxiter=20000, warm=None,
                  solver="native-jacobi"):
        """State + one adjoint per functional, sharing one configured matrix.

        The self-adjoint operator serves every solve; ``warm`` is an optional
        :class:`FunctionalLinearization` from the previous design iterate
        whose fields seed all CG solves.  Sensitivities come from one
        element-wise ``Integrate`` per functional:
        ``jacobians[k] = -Integrate(lambda_k . m, element_wise=True)``.
        """
        loads = list(functional_loads)
        if not loads:
            raise ValueError("DensityAdjointVIM.linearize: no functional loads")
        if warm is not None and len(warm.gfLambdas) != len(loads):
            raise ValueError(
                "DensityAdjointVIM.linearize: warm start carries %d adjoints "
                "for %d loads" % (len(warm.gfLambdas), len(loads)))
        mass, A, pre = self._system(
            s, need_ngsolve_solver=(solver == "ngsolve-cg"))
        if solver in {"native", "native-jacobi"}:
            rhs = [state_load.vec] + [load.vec for load in loads]
            warm_vecs = None if warm is None else [warm.gfM.vec] + [
                value.vec for value in warm.gfLambdas]
            fields, native_iterations = self._native_solve_many(
                mass, rhs, tol, maxiter, warm_vecs=warm_vecs,
                mass_riesz=(solver == "native"))
            gfm, gfls = fields[0], fields[1:]
            it_m, its = native_iterations[0], native_iterations[1:]
        elif solver == "ngsolve-cg":
            gfm = ng.GridFunction(self.fes)
            it_m = self._cg(A, pre, state_load.vec, gfm, tol, maxiter,
                            warm_vec=None if warm is None else warm.gfM.vec)
            gfls, its = [], []
            for k, load in enumerate(loads):
                gfl = ng.GridFunction(self.fes)
                wv = None if warm is None else warm.gfLambdas[k].vec
                its.append(self._cg(A, pre, load.vec, gfl, tol, maxiter,
                                    warm_vec=wv))
                gfls.append(gfl)
        else:
            raise ValueError(
                "DensityAdjointVIM.linearize: solver must be "
                "'native-jacobi', 'native', or 'ngsolve-cg'")
        values = np.array([float(ng.InnerProduct(load.vec, gfm.vec))
                           for load in loads])
        jacobians = np.stack([
            -np.asarray(ng.Integrate(gfl * gfm, self.mesh,
                                     element_wise=True), dtype=float)
            for gfl in gfls])
        return FunctionalLinearization(
            values=values, jacobians=jacobians, gfM=gfm,
            gfLambdas=tuple(gfls), state_iterations=it_m,
            adjoint_iterations=tuple(its))

    def objective_and_gradient(self, s, state_load, adjoint_load,
                               tol=1e-12, maxiter=20000, warm=None,
                               solver="native-jacobi"):
        """State + adjoint solve and the full per-element gradient.

        One-functional convenience over :meth:`linearize` (same configured
        operator).  ``warm`` is an optional
        :class:`AdjointGradientResult` from the previous design iterate; its
        fields seed both CG solves.  The objective is the design-dependent
        part ``adjoint_load^T m``; the gradient is
        ``dJ/ds_e = -Integrate(lambda . m, element_wise=True)``.
        """
        warm_lin = None
        if warm is not None:
            warm_lin = FunctionalLinearization(
                values=np.array([warm.objective]),
                jacobians=warm.gradient.reshape(1, -1), gfM=warm.gfM,
                gfLambdas=(warm.gfLambda,),
                state_iterations=warm.state_iterations,
                adjoint_iterations=(warm.adjoint_iterations,))
        lin = self.linearize(s, state_load, [adjoint_load], tol=tol,
                             maxiter=maxiter, warm=warm_lin, solver=solver)
        return AdjointGradientResult(
            objective=float(lin.values[0]), gradient=lin.jacobians[0],
            gfM=lin.gfM, gfLambda=lin.gfLambdas[0],
            state_iterations=lin.state_iterations,
            adjoint_iterations=lin.adjoint_iterations[0])


# --------------------------------------------------------------------------
# constrained design loop (trust-region SLP)
# --------------------------------------------------------------------------
@dataclass
class DensityDesignResult:
    """Converged/terminal state of :func:`optimize_density`."""
    density: np.ndarray       # raw (unfiltered) element densities in [0, 1]
    history: tuple            # one dict per ACCEPTED iterate
    converged: bool           # move limit collapsed below move_min
    final_move: float
    solves: int               # linearizations evaluated (incl. rejected trials)


def optimize_density(problem, state_load, objective_load, constraint_loads=(),
                     targets=(), *, chi_iron, volume_fraction,
                     density_filter=None, density_projection=None,
                     initial_density=None, penalty=1.0,
                     move_limit=0.1, max_iterations=30,
                     band_relative=5e-3, band_floor=None,
                     move_min=1e-3, objective_slack=1e-6,
                     tol=1e-10, cg_maxiter=20000,
                     linear_solver="native-jacobi", callback=None,
                     checkpoint_callback=None, initial_warm=None,
                     evaluation_callback=None):
    """MAXIMIZE ``J = objective_load^T m`` under linear-functional equality
    bands, an iron volume budget, and box/move limits (trust-region SLP).

    Constraints are engineering bands: ``|constraint_k - target_k| <=
    band_k`` with ``band_k = band_relative * |target_k|`` (or explicit
    ``band_floor``, absolute, per constraint).  Each iterate costs one
    weighted-mass assembly and ``1 + n_constraints``
    warm-started CG solves (:meth:`DensityAdjointVIM.linearize`), one
    element-wise Integrate per functional, and a milliseconds HiGHS LP
    (:func:`radia.topology_optimization.solve_lp_update` with the constraint
    rows in its ``A_ub`` slot, normalized to O(1) -- Tesla-scale rows sit
    below HiGHS's absolute feasibility tolerance and read as noise).

    Two-phase trust-region SLP.  While every violation is inside its band,
    the LP maximizes J and acceptance keeps the ascent MONOTONE and the
    violations BOUNDED: a trial is accepted only if J does not decrease
    (beyond ``objective_slack`` relative) and every violation is either
    under the absolute cap ``1.25 * band`` or in strict geometric decrease.
    While any violation sits OUTSIDE its band (e.g. an infeasible profile
    start), the loop switches to RESTORATION: the LP objective becomes the
    steepest combined violation descent (J free), the rows hold the
    non-worsening guard ``max(band, |viol|)`` (always feasible at the
    current point -- the LP can never be infeasible), and acceptance
    requires the band-weighted total violation to decrease (individual
    violations may not blow up while others improve).  Rejected trials
    halve the move limit against the SAME linearization.  When a constraint
    is active at the optimum the design rides its band edge; measured
    behavior on the verification cases (2026-07-28): strictly monotone J
    (+1.3 % ball / +19 % sector surrogate), violations peaking at 1.24 x
    band and riding at ~1.05 x band.

    ``callback(entry)`` receives each accepted history dict.  The caller
    wraps the whole call in ``with TaskManager():``.  Informal per-iterate
    timings are recorded in the history; benchmark-grade timings belong on
    the quiet compute hosts, per the repository benchmark policy.
    """
    from .topology_optimization import solve_lp_update

    constraint_loads = list(constraint_loads)
    targets = np.asarray(targets, dtype=float).ravel()
    if len(constraint_loads) != targets.size:
        raise ValueError("optimize_density: %d constraint loads but %d targets"
                         % (len(constraint_loads), targets.size))
    if not 0.0 < volume_fraction <= 1.0:
        raise ValueError("optimize_density: volume_fraction must be in (0, 1]")
    move_limit = float(move_limit)
    move_min = float(move_min)
    if (not np.isfinite(move_limit) or not np.isfinite(move_min)
            or not 0.0 < move_min <= move_limit <= 1.0):
        raise ValueError("optimize_density: need 0 < move_min <= "
                         "move_limit <= 1")
    max_iterations_i = int(max_iterations)
    if max_iterations_i < 1 or max_iterations_i != max_iterations:
        raise ValueError("optimize_density: max_iterations must be a positive integer")
    objective_slack = float(objective_slack)
    if not np.isfinite(objective_slack) or objective_slack < 0.0:
        raise ValueError("optimize_density: objective_slack must be non-negative")
    if not np.all(np.isfinite(targets)):
        raise ValueError("optimize_density: targets must be finite")
    if linear_solver not in {"native-jacobi", "native", "ngsolve-cg"}:
        raise ValueError(
            "optimize_density: linear_solver must be 'native-jacobi', "
            "'native', or 'ngsolve-cg'")
    volumes = problem.element_volumes
    volume_max = float(volume_fraction * volumes.sum())
    if initial_density is None:
        rho = np.full(problem.n_el, float(volume_fraction))
    else:
        rho = np.asarray(initial_density, dtype=float).copy()
        if rho.shape != (problem.n_el,):
            raise ValueError("optimize_density: initial_density shape %r != "
                             "(%d,)" % (rho.shape, problem.n_el))
    if (not np.all(np.isfinite(rho)) or np.any(rho < 0.0)
            or np.any(rho > 1.0)):
        raise ValueError("optimize_density: initial_density must be finite "
                         "and lie in [0, 1]")
    initial_volume = float(volumes @ rho)
    if initial_volume > volume_max + 1e-12 * max(1.0, abs(volume_max)):
        raise ValueError("optimize_density: initial_density exceeds the "
                         "iron volume budget")
    if band_floor is None:
        if not np.isfinite(band_relative) or band_relative <= 0.0:
            raise ValueError("optimize_density: band_relative must be positive")
        band = float(band_relative) * np.maximum(np.abs(targets), 1e-300)
    else:
        band = np.broadcast_to(np.asarray(band_floor, dtype=float),
                               targets.shape).astype(float)
    if targets.size and not np.all(band > 0.0):
        raise ValueError("optimize_density: constraint bands must be positive")
    loads = [objective_load] + constraint_loads

    evaluation_index = 0

    def evaluate(rho_vec, warm):
        nonlocal evaluation_index
        t_evaluate = time.perf_counter()
        # The P1 Helmholtz filter under/overshoots at bang-bang transitions
        # (measured -1.2e-2 on a coarse ball); clip the FILTERED density to
        # [0, 1] with the exact piecewise chain rule (zero derivative on
        # clipped elements) before the material map.
        if density_filter is not None:
            rho_f_raw = density_filter.apply(rho_vec)
            rho_f = np.clip(rho_f_raw, 0.0, 1.0)
            unclipped = (rho_f_raw > 0.0) & (rho_f_raw < 1.0)
        else:
            rho_f = rho_vec
            unclipped = None
        if density_projection is not None:
            rho_material = density_projection.apply(rho_f)
        else:
            rho_material = rho_f
        lin = problem.linearize(density_to_s(
            rho_material, chi_iron, penalty=penalty),
                                state_load, loads, tol=tol,
                                maxiter=cg_maxiter, warm=warm,
                                solver=linear_solver)

        def to_rho(g_s):
            g_rf = density_gradient_from_s_gradient(
                rho_material, g_s, chi_iron, penalty=penalty)
            if density_projection is not None:
                g_rf = density_projection.chain(rho_f, g_rf)
            if density_filter is None:
                return g_rf
            return density_filter.chain(np.where(unclipped, g_rf, 0.0))

        gJ = to_rho(lin.jacobians[0])
        gks = [to_rho(lin.jacobians[1 + k])
               for k in range(len(constraint_loads))]
        if evaluation_callback is not None:
            evaluation_callback(dict(
                evaluation=evaluation_index,
                elapsed_s=time.perf_counter() - t_evaluate,
                warm_start=warm is not None,
                state_iterations=int(lin.state_iterations),
                adjoint_iterations=[int(v) for v in lin.adjoint_iterations],
                values=np.asarray(lin.values, dtype=float).tolist()))
        evaluation_index += 1
        return lin, gJ, gks

    lin, gJ, gks = evaluate(rho, initial_warm)
    n_solves = 1
    history = []
    move = move_limit
    converged = False
    for iteration in range(max_iterations_i):
        t_iter = time.perf_counter()
        J = float(lin.values[0])
        viol = lin.values[1:] - targets
        # Three-tier SLP.  DEEP RESTORATION (any violation beyond the 1.25*
        # band acceptance cap, e.g. an infeasible profile start): the LP
        # OBJECTIVE becomes the steepest combined violation descent and J is
        # free -- putting the shrink demand in the constraint ROWS instead
        # goes infeasible once the move box shrinks and dead-ends (measured:
        # zero accepted iterates on an isochronous-profile start).  Inside
        # the cap, the LP maximizes J with the tested band rows: violations
        # in the (band, 1.25 band] transition zone get a geometric pull-back
        # row (fall back to a non-worsening hold row when that is
        # unreachable in the move box), and acceptance keeps J MONOTONE.
        deep = bool(targets.size) and bool(np.any(np.abs(viol)
                                                  > 1.25 * band))

        def lp_rows(bands_eff):
            # rows normalized to O(1) for HiGHS's absolute tolerances
            A_ub, b_ub = [], []
            for G, v, b in zip(gks, viol, bands_eff):
                scale = 1.0 / max(np.abs(G).max(), 1e-300)
                base = float(G @ rho)
                A_ub.append(G * scale)
                b_ub.append((b - v + base) * scale)
                A_ub.append(-G * scale)
                b_ub.append((b + v - base) * scale)
            if not A_ub:
                return None, None
            return np.array(A_ub), np.array(b_ub)

        if deep:
            direction = np.zeros_like(gJ)
            for G, v, b in zip(gks, viol, band):
                if abs(v) > 1.25 * b:
                    direction += (np.sign(v) / max(np.abs(G).max(), 1e-300)) * G
            lp_objective = direction / max(np.abs(direction).max(), 1e-300)
        else:
            lp_objective = -gJ / max(np.abs(gJ).max(), 1e-300)

        accepted = False
        trials = 0
        band_mode = "deep" if deep else "band"
        while move >= move_min:
            trials += 1
            if deep:
                bands_eff = np.maximum(band, np.abs(viol))  # always feasible
                A_ub, b_ub = lp_rows(bands_eff)
                update = solve_lp_update(rho, lp_objective, volumes,
                                         volume_max, move_limit=move,
                                         A_ub=A_ub, b_ub=b_ub)
            else:
                bands_eff = np.where(np.abs(viol) > band,
                                     np.maximum(band, 0.9 * np.abs(viol)),
                                     band)
                band_mode = ("band" if not targets.size
                             or np.all(np.abs(viol) <= band) else "restore")
                A_ub, b_ub = lp_rows(bands_eff)
                try:
                    update = solve_lp_update(rho, lp_objective, volumes,
                                             volume_max, move_limit=move,
                                             A_ub=A_ub, b_ub=b_ub)
                except RuntimeError:
                    bands_eff = np.maximum(band, np.abs(viol))
                    band_mode = "hold"
                    A_ub, b_ub = lp_rows(bands_eff)
                    update = solve_lp_update(rho, lp_objective, volumes,
                                             volume_max, move_limit=move,
                                             A_ub=A_ub, b_ub=b_ub)
            lin_new, gJ_new, gks_new = evaluate(update.density, warm=lin)
            n_solves += 1
            J_new = float(lin_new.values[0])
            viol_new = np.abs(lin_new.values[1:] - targets)
            if deep:
                # accept on feasibility progress: the band-weighted total
                # violation decreases, or full cap entry; individual
                # violations must not blow up while others improve.
                total = float(np.sum(np.abs(viol) / band))
                total_new = float(np.sum(viol_new / band))
                ok = ((total_new < 0.995 * total
                       or bool(np.all(viol_new <= 1.25 * band)))
                      and bool(np.all(viol_new
                                      <= np.maximum(1.25 * band,
                                                    1.05 * np.abs(viol)))))
            else:
                ok_J = J_new >= J - objective_slack * abs(J)
                # absolute cap 1.25*band at all times OR strict geometric
                # decrease while outside -- no ratchet path exists.
                ok_g = np.all((viol_new <= 1.25 * band)
                              | (viol_new <= 0.97 * np.abs(viol)))
                ok = ok_J and bool(ok_g)
            if ok:
                accepted = True
                break
            move *= 0.5
        if not accepted:
            converged = True
            break
        change = float(np.max(np.abs(update.delta)))
        rho, lin, gJ, gks = update.density, lin_new, gJ_new, gks_new
        move = min(float(move_limit), 1.5 * move)
        entry = dict(iteration=iteration, objective=float(lin.values[0]),
                     constraints=lin.values[1:].tolist(),
                     violation=viol_new.tolist(), band=band.tolist(),
                     volume=float(volumes @ rho), max_density_change=change,
                     move=move, trials=trials, band_mode=band_mode,
                     t_iter_s=time.perf_counter() - t_iter,
                     state_iterations=lin.state_iterations,
                     adjoint_iterations=list(lin.adjoint_iterations))
        if density_filter is not None:
            rho_report = np.clip(density_filter.apply(rho), 0.0, 1.0)
        else:
            rho_report = rho
        if density_projection is not None:
            rho_report = density_projection.apply(rho_report)
        entry.update(density_discreteness(rho_report))
        history.append(entry)
        if callback is not None:
            callback(entry)
        if checkpoint_callback is not None:
            checkpoint_callback(entry, rho.copy())
    return DensityDesignResult(density=rho, history=tuple(history),
                               converged=converged, final_move=move,
                               solves=n_solves)


# --------------------------------------------------------------------------
# Stage-3 verification protocol: exact-void iron-only re-evaluation
# --------------------------------------------------------------------------
# Boundary triangles of a netgen tet, ordered so the right-hand-rule normal
# points OUT of the element -- netgen stores surface elements outward for
# FaceDescriptor(domin=1, domout=0) (probed on an OCC sphere mesh).  The
# opposite handedness flips every surface charge and turns the demag solve
# into runaway magnetization (measured <Mz> 18-34 instead of ~2.2).
_TET_BOUNDARY_FACES = ((0, 3, 1), (1, 3, 2), (2, 3, 0), (0, 1, 2))


def iron_only_mesh(mesh, keep):
    """New straight-tet netgen/NGSolve mesh from the kept VOL elements.

    Exact void removal: vertices of the kept set are copied, kept tets are
    re-added as one ``iron`` material, and the boundary facets of the kept
    set (facets owned by exactly one kept element) become the new exterior
    surface.  ``keep`` is a boolean mask in NGSolve VOL element numbering;
    the numbering correspondence with the netgen element list is verified
    element-by-element (fail loud on mismatch).
    """
    import netgen.meshing as nm

    keep = np.asarray(keep, dtype=bool).ravel()
    if mesh.dim != 3:
        raise NotImplementedError("iron_only_mesh: 3D meshes only")
    if mesh.GetCurveOrder() >= 2:
        raise NotImplementedError(
            "iron_only_mesh: curved parent meshes are not supported (the "
            "extraction copies vertices only); pass the straight design mesh")
    src = mesh.ngmesh
    points = list(src.Points())
    elements = list(src.Elements3D())
    if keep.size != len(elements):
        raise ValueError("iron_only_mesh: mask has %d entries, mesh has %d "
                         "volume elements" % (keep.size, len(elements)))
    if not keep.any() or keep.all():
        raise ValueError("iron_only_mesh: the kept set must be a proper "
                         "non-empty subset (got %d of %d elements)"
                         % (int(keep.sum()), keep.size))
    # the mask lives in NGSolve VOL numbering -- verify it matches the
    # netgen element list order before using it
    for ng_el, nm_el in zip(mesh.Elements(ng.VOL), elements):
        if (tuple(sorted(v.nr for v in ng_el.vertices))
                != tuple(sorted(v.nr - 1 for v in nm_el.vertices))):
            raise RuntimeError(
                "iron_only_mesh: NGSolve VOL element numbering does not "
                "match the netgen Elements3D order on this mesh")
    new = nm.Mesh(3)
    new.SetMaterial(1, "iron")
    pmap = {}

    def pid(nr):
        if nr not in pmap:
            p = points[nr - 1].p
            pmap[nr] = new.Add(nm.MeshPoint(nm.Pnt(p[0], p[1], p[2])))
        return pmap[nr]

    facets = {}
    for index in np.flatnonzero(keep):
        vs = [v.nr for v in elements[index].vertices]
        if len(vs) != 4:
            raise NotImplementedError("iron_only_mesh: TET meshes only")
        new.Add(nm.Element3D(1, [pid(v) for v in vs]))
        for fa in _TET_BOUNDARY_FACES:
            tri = (vs[fa[0]], vs[fa[1]], vs[fa[2]])
            facets.setdefault(tuple(sorted(tri)), []).append(tri)
    descriptor = new.Add(nm.FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    for occurrences in facets.values():
        if len(occurrences) == 1:
            new.Add(nm.Element2D(descriptor,
                                 [pid(v) for v in occurrences[0]]))
        elif len(occurrences) != 2:
            raise RuntimeError("iron_only_mesh: facet shared by %d kept "
                               "elements" % len(occurrences))
    return ng.Mesh(new)


@dataclass
class IronOnlyVerification:
    """Exact-void re-evaluation of a thresholded design (Stage-3 protocol).

    ``values_*[k]`` are the functional values (state solve + load inner
    product) for the supplied builders; ``bands[k] = (embedded - iron_only)
    / |iron_only|`` quantifies the ersatz-void error of the embedded model
    honestly, per functional, at matched 0/1 representation.
    """
    keep: np.ndarray               # thresholded iron mask (parent numbering)
    iron_mesh: "ng.Mesh"
    values_embedded: np.ndarray    # 0/1 ersatz void on the parent mesh
    values_iron_only: np.ndarray   # exact void on the extracted mesh
    bands: np.ndarray
    embedded_iterations: int
    iron_only_iterations: int


def verify_design_iron_only(problem, density, state_load_builder,
                            functional_builders, *, chi_iron, threshold=0.5,
                            density_filter=None, chi_min=CHI_MIN,
                            tol=1e-10, cg_maxiter=5000, gram_kwargs=None):
    """Stage-3 final-verification protocol for a converged density design.

    Thresholds the (filtered) density at ``threshold``, evaluates every
    functional on (a) the PARENT mesh with the 0/1 ersatz void
    (``chi_min``) and (b) a NEW iron-only mesh with the void REMOVED
    (:func:`iron_only_mesh` -- the exact-void gold standard), and reports
    the per-functional ersatz band.  Loads are mesh-bound, so the caller
    passes BUILDERS: ``state_load_builder(fes)`` and
    ``functional_builders = [f(fes) -> assembled load, ...]`` are invoked
    on both spaces (orbit points and weights are geometry-fixed).

    The comparison is at MATCHED 0/1 representation: it quantifies the
    ersatz-void/interface-charge error alone.  The separate gap between
    the CONTINUOUS design objective and the thresholded one (large for
    gray designs; drive it down with ``penalty``/projection before
    manufacturing conclusions) is the caller's report from the design
    history.  Never report final design numbers from the embedded model --
    this function exists to replace them.
    """
    density = np.asarray(density, dtype=float).ravel()
    if density.shape != (problem.n_el,):
        raise ValueError("verify_design_iron_only: density shape %r != (%d,)"
                         % (density.shape, problem.n_el))
    rho_f = (np.clip(density_filter.apply(density), 0.0, 1.0)
             if density_filter is not None else density)
    keep = rho_f >= float(threshold)
    loads = [b(problem.fes) for b in functional_builders]
    s_bin = np.where(keep, 1.0 / chi_iron, 1.0 / chi_min)
    gf_emb, it_emb = problem.solve(s_bin, state_load_builder(problem.fes),
                                   tol=tol, maxiter=cg_maxiter)
    values_emb = np.array([float(ng.InnerProduct(load.vec, gf_emb.vec))
                           for load in loads])
    mesh_iron = iron_only_mesh(problem.mesh, keep)
    fes_iron = ng.HDiv(mesh_iron, order=int(problem.fes.globalorder))
    problem_iron = DensityAdjointVIM(fes_iron, **(gram_kwargs or {}))
    loads_iron = [b(fes_iron) for b in functional_builders]
    gf_iron, it_iron = problem_iron.solve(
        np.full(problem_iron.n_el, 1.0 / chi_iron),
        state_load_builder(fes_iron), tol=tol, maxiter=cg_maxiter)
    values_iron = np.array([float(ng.InnerProduct(load.vec, gf_iron.vec))
                            for load in loads_iron])
    bands = (values_emb - values_iron) / np.maximum(np.abs(values_iron),
                                                    1e-300)
    return IronOnlyVerification(
        keep=keep, iron_mesh=mesh_iron, values_embedded=values_emb,
        values_iron_only=values_iron, bands=bands,
        embedded_iterations=it_emb, iron_only_iterations=it_iron)
