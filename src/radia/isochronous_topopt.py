"""Per-element density adjoint for HDiv-MMM topology optimization.

Application-layer module for isochronous-magnet (and general accelerator
pole/yoke) density topology optimization on the HDiv-MMM forward engine
(design record: ``docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md``).
It composes :class:`radia.vim.DemagOperator` (build-once geometry operator
``N = B^T G B``) with the per-element design variable ``s_e = 1/chi_e``
(an ``L2(order=0)`` mass weight) into the adjoint gradient route:

    (M_s + N) m      = f_state        state solve (SPD; CG + mass-Riesz)
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

Verified 2026-07-28 (research run ``C:/temp/vim_topopt/stage1_adjoint_gate.py``,
unit ball, 270 tets, log-uniform ``s`` in [1e-2, 1]): adjoint gradient matches
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
    "density_to_s", "density_gradient_from_s_gradient",
    "gradient_pair_points", "dipole_array_field_cf",
    "field_functional_load", "uniform_field_load",
    "demag_field_from_solution", "orbit_arc_points", "optimize_density",
]


# --------------------------------------------------------------------------
# density <-> s mapping
# --------------------------------------------------------------------------
def density_to_s(density, chi_iron, chi_min=CHI_MIN):
    """Map a per-element density ``rho in [0, 1]`` to ``s = 1/chi(rho)``.

    Linear interpolation in susceptibility, ``chi(rho) = chi_min +
    (chi_iron - chi_min) rho`` -- the penalization/filtering layer of the
    design loop operates on ``rho`` before this map.  ``chi_min`` defaults to
    the validated ersatz floor :data:`CHI_MIN`.
    """
    rho = np.asarray(density, dtype=float)
    if not chi_iron > chi_min > 0.0:
        raise ValueError(
            "density_to_s: need chi_iron > chi_min > 0 (got chi_iron=%r, chi_min=%r)"
            % (chi_iron, chi_min))
    if np.any(rho < -1e-9) or np.any(rho > 1.0 + 1e-9):
        raise ValueError(
            "density_to_s: density must lie in [0, 1] (got min=%r, max=%r)"
            % (float(rho.min()), float(rho.max())))
    chi = chi_min + (chi_iron - chi_min) * np.clip(rho, 0.0, 1.0)
    return 1.0 / chi


def density_gradient_from_s_gradient(density, s_gradient, chi_iron,
                                     chi_min=CHI_MIN):
    """Chain rule ``dJ/drho_e = dJ/ds_e * ds/drho_e`` for :func:`density_to_s`.

    ``ds/drho = -(chi_iron - chi_min)/chi(rho)^2``.
    """
    rho = np.clip(np.asarray(density, dtype=float), 0.0, 1.0)
    grad_s = np.asarray(s_gradient, dtype=float)
    if rho.shape != grad_s.shape:
        raise ValueError("density_gradient_from_s_gradient: shape mismatch %r vs %r"
                         % (rho.shape, grad_s.shape))
    chi = chi_min + (chi_iron - chi_min) * rho
    return grad_s * (-(chi_iron - chi_min) / (chi * chi))


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
    """State + one adjoint per functional, all on one factorization.

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
    weighted-mass assembly + factorization and two CG solves.  The caller
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
    def _system(self, s):
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

    # ------------------------------------------------------------ public API
    def solve(self, s, load, tol=1e-12, maxiter=5000, warm=None):
        """One SPD solve ``(M_s + N) x = load.vec``.

        ``warm`` is an optional GridFunction whose vector seeds CG (warm start
        across design iterates).  Returns ``(GridFunction, iterations)``.
        """
        mass, A, pre = self._system(s)
        gf = ng.GridFunction(self.fes)
        iters = self._cg(A, pre, load.vec, gf, tol, maxiter,
                         warm_vec=None if warm is None else warm.vec)
        return gf, iters

    def linearize(self, s, state_load, functional_loads,
                  tol=1e-12, maxiter=5000, warm=None):
        """State + one adjoint per functional, sharing one factorization.

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
        mass, A, pre = self._system(s)
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
                               tol=1e-12, maxiter=5000, warm=None):
        """State + adjoint solve and the full per-element gradient.

        One-functional convenience over :meth:`linearize` (same operator,
        same factorization).  ``warm`` is an optional
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
                             maxiter=maxiter, warm=warm_lin)
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
                     density_filter=None, initial_density=None,
                     move_limit=0.1, max_iterations=30,
                     band_relative=5e-3, band_floor=None, restore_shrink=0.1,
                     move_min=1e-3, objective_slack=1e-6,
                     tol=1e-10, cg_maxiter=5000, callback=None):
    """MAXIMIZE ``J = objective_load^T m`` under linear-functional equality
    bands, an iron volume budget, and box/move limits (trust-region SLP).

    Constraints are engineering bands: ``|constraint_k - target_k| <=
    band_k`` with ``band_k = band_relative * |target_k|`` (or explicit
    ``band_floor``, absolute, per constraint).  Each iterate costs one
    weighted-mass assembly + factorization and ``1 + n_constraints``
    warm-started CG solves (:meth:`DensityAdjointVIM.linearize`), one
    element-wise Integrate per functional, and a milliseconds HiGHS LP
    (:func:`radia.topology_optimization.solve_lp_update` with the constraint
    rows in its ``A_ub`` slot, normalized to O(1) -- Tesla-scale rows sit
    below HiGHS's absolute feasibility tolerance and read as noise).

    Trust-region acceptance makes the ascent MONOTONE and the violations
    BOUNDED: a trial is accepted only if J does not decrease (beyond
    ``objective_slack`` relative) and every violation is either under the
    absolute cap ``1.25 * band`` or in strict geometric decrease (outside
    the band); rejected trials halve the move limit against the SAME
    linearization.  When a constraint is active at the optimum the design
    rides its band edge; measured behavior on the verification cases
    (2026-07-28): strictly monotone J (+1.3 % ball / +19 % sector surrogate),
    violations peaking at 1.24 x band and riding at ~1.05 x band.

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
    if not 0.0 < restore_shrink < 1.0:
        raise ValueError("optimize_density: restore_shrink must be in (0, 1)")
    volumes = problem.element_volumes
    volume_max = float(volume_fraction * volumes.sum())
    if initial_density is None:
        rho = np.full(problem.n_el, float(volume_fraction))
    else:
        rho = np.asarray(initial_density, dtype=float).copy()
        if rho.shape != (problem.n_el,):
            raise ValueError("optimize_density: initial_density shape %r != "
                             "(%d,)" % (rho.shape, problem.n_el))
    if band_floor is None:
        band = float(band_relative) * np.maximum(np.abs(targets), 1e-300)
    else:
        band = np.broadcast_to(np.asarray(band_floor, dtype=float),
                               targets.shape).astype(float)
    if targets.size and not np.all(band > 0.0):
        raise ValueError("optimize_density: constraint bands must be positive")
    loads = [objective_load] + constraint_loads

    def evaluate(rho_vec, warm):
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
        lin = problem.linearize(density_to_s(rho_f, chi_iron), state_load,
                                loads, tol=tol, maxiter=cg_maxiter, warm=warm)

        def to_rho(g_s):
            g_rf = density_gradient_from_s_gradient(rho_f, g_s, chi_iron)
            if density_filter is None:
                return g_rf
            return density_filter.chain(np.where(unclipped, g_rf, 0.0))

        gJ = to_rho(lin.jacobians[0])
        gks = [to_rho(lin.jacobians[1 + k])
               for k in range(len(constraint_loads))]
        return lin, gJ, gks

    lin, gJ, gks = evaluate(rho, None)
    n_solves = 1
    history = []
    move = float(move_limit)
    converged = False
    for iteration in range(int(max_iterations)):
        t_iter = time.perf_counter()
        J = float(lin.values[0])
        viol = lin.values[1:] - targets
        gJ_lp = -gJ / max(np.abs(gJ).max(), 1e-300)

        def lp_rows(bands_eff):
            # normalized to O(1) for HiGHS's absolute tolerances
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

        accepted = False
        trials = 0
        band_mode = "band"
        while move >= move_min:
            trials += 1
            # absolute band; RESTORE (geometric shrink back) when a violation
            # sits outside it; HOLD (non-worsening, feasible at x = rho by
            # construction) when restoration is unreachable in this move box.
            bands_eff = np.where(
                np.abs(viol) > band,
                np.maximum(band, (1.0 - restore_shrink) * np.abs(viol)), band)
            band_mode = ("band" if not targets.size
                         or np.all(np.abs(viol) <= band) else "restore")
            A_ub, b_ub = lp_rows(bands_eff)
            try:
                update = solve_lp_update(rho, gJ_lp, volumes, volume_max,
                                         move_limit=move, A_ub=A_ub, b_ub=b_ub)
            except RuntimeError:
                bands_eff = np.maximum(band, np.abs(viol))
                band_mode = "hold"
                A_ub, b_ub = lp_rows(bands_eff)
                update = solve_lp_update(rho, gJ_lp, volumes, volume_max,
                                         move_limit=move, A_ub=A_ub, b_ub=b_ub)
            lin_new, gJ_new, gks_new = evaluate(update.density, warm=lin)
            n_solves += 1
            J_new = float(lin_new.values[0])
            viol_new = np.abs(lin_new.values[1:] - targets)
            ok_J = J_new >= J - objective_slack * abs(J)
            # absolute cap 1.25*band at all times OR strict geometric decrease
            # while outside the band -- no ratchet path exists.
            ok_g = np.all((viol_new <= 1.25 * band)
                          | (viol_new <= 0.97 * np.abs(viol)))
            if ok_J and ok_g:
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
        history.append(entry)
        if callback is not None:
            callback(entry)
    return DensityDesignResult(density=rho, history=tuple(history),
                               converged=converged, final_move=move,
                               solves=n_solves)
