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
    "density_to_s", "density_gradient_from_s_gradient",
    "gradient_pair_points", "dipole_array_field_cf",
    "field_functional_load", "uniform_field_load",
    "demag_field_from_solution",
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
def gradient_pair_points(points, weights, delta, axis=0):
    """Realize ``sum_i w_i d(field)/dx_axis (x_i)`` as +/- evaluation pairs.

    Returns ``(points_out, weights_out)`` with ``2 N`` rows: each input point
    splits into ``x_i +/- (delta/2) e_axis`` carrying weights ``+/- w_i/delta``
    (central finite difference of the point functional; ``delta`` is a
    physical stencil length on the orbit, not a numerical epsilon).
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    wts = np.asarray(weights, dtype=float).ravel()
    if len(pts) != len(wts):
        raise ValueError("gradient_pair_points: %d points but %d weights"
                         % (len(pts), len(wts)))
    if not delta > 0.0:
        raise ValueError("gradient_pair_points: delta must be positive")
    if axis not in (0, 1, 2):
        raise ValueError("gradient_pair_points: axis must be 0, 1, or 2")
    shift = np.zeros(3)
    shift[axis] = 0.5 * float(delta)
    return (np.concatenate([pts + shift, pts - shift]),
            np.concatenate([wts / delta, -wts / delta]))


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

    def objective_and_gradient(self, s, state_load, adjoint_load,
                               tol=1e-12, maxiter=5000, warm=None):
        """State + adjoint solve and the full per-element gradient.

        One weighted-mass assembly and one factorization serve both solves
        (the operator is self-adjoint).  ``warm`` is an optional
        :class:`AdjointGradientResult` from the previous design iterate; its
        fields seed both CG solves.  The objective is the design-dependent
        part ``adjoint_load^T m``; the gradient is
        ``dJ/ds_e = -Integrate(lambda . m, element_wise=True)``.
        """
        mass, A, pre = self._system(s)
        gfm = ng.GridFunction(self.fes)
        gfl = ng.GridFunction(self.fes)
        it_m = self._cg(A, pre, state_load.vec, gfm, tol, maxiter,
                        warm_vec=None if warm is None else warm.gfM.vec)
        it_l = self._cg(A, pre, adjoint_load.vec, gfl, tol, maxiter,
                        warm_vec=None if warm is None else warm.gfLambda.vec)
        objective = float(ng.InnerProduct(adjoint_load.vec, gfm.vec))
        gradient = -np.asarray(
            ng.Integrate(gfl * gfm, self.mesh, element_wise=True), dtype=float)
        return AdjointGradientResult(
            objective=objective, gradient=gradient, gfM=gfm, gfLambda=gfl,
            state_iterations=it_m, adjoint_iterations=it_l)
