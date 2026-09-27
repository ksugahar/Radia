"""Nonlinear transient heat conduction for induction-heating workpieces.

Enthalpy form of backward Euler, per time step ``t_n -> t_{n+1}``::

    int (H(T) - H(T_n)) / dt v w  +  int k(T) grad T . grad v w
      + int_conv h (T - T_ext) v w  +  int_rad eps sigma (T^4 - T_ext^4) v w
      = int_heat q(T) v w

with ``w = 1`` (3D) or ``w = 2 pi r`` (axisymmetric).  Writing the storage
term through the enthalpy ``H`` makes the discrete energy balance exact
for any ``cp(T)`` and latent heat.  Each step is solved by a Newton-type
iteration whose matrix is symmetric, so a Cholesky factorisation can be
used: the storage term uses the secant capacity
``(H(T_k) - H(T_n)) / (T_k - T_n)`` (the tangent ``dH/dT`` on the first
iterate), which converges monotonically for any non-decreasing ``H`` even
across the kinks of a latent band; radiation is linearised exactly;
``k(T)`` and ``q(T)`` are lagged at the current iterate.

The residual, the matrix and the energy audit all use the same explicit
integration rules.  NGSolve otherwise picks the quadrature order of each
form from its integrand, and a residual integrated more accurately than
its matrix (a nonlinear ``H(T)`` times the axisymmetric weight ``r``) made
the iteration diverge geometrically after the first few steps.  A step
that does not converge is split in halves, up to ``max_halvings`` times;
beyond that the run fails.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

SIGMA_SB = 5.670374419e-8
KELVIN = 273.15


@dataclass
class HeatBoundaryTerms:
    heat_flux: str = ""
    convection: str = ""
    h_conv: float = 0.0
    t_ext: float = 20.0
    radiation: str = ""
    emissivity: float = 0.0


@dataclass
class TransientAudit:
    steps: int = 0
    substeps: int = 0
    halvings: int = 0
    newton_iterations: list = field(default_factory=list)
    max_newton_iterations: int = 0
    energy_in_J: float = 0.0
    energy_loss_J: float = 0.0
    energy_stored_J: float = 0.0
    table_extrapolation_C: float = 0.0

    def as_dict(self):
        balance = self.energy_in_J - self.energy_loss_J
        rel = ((self.energy_stored_J - balance) / abs(balance)
               if balance else 0.0)
        return {
            "scheme": "enthalpy-backward-euler-newton",
            "steps": self.steps, "substeps": self.substeps,
            "halvings": self.halvings,
            "max_newton_iterations": self.max_newton_iterations,
            "mean_newton_iterations": (float(np.mean(self.newton_iterations))
                                       if self.newton_iterations else 0.0),
            "energy_in_J": self.energy_in_J,
            "energy_loss_J": self.energy_loss_J,
            "energy_stored_J": self.energy_stored_J,
            "energy_balance_relative_error": float(rel),
            "table_extrapolation_C": self.table_extrapolation_C,
        }


class NonlinearHeatStepper:
    """Advance ``gfT`` with temperature-dependent properties and source."""

    def __init__(self, gfT, material, boundaries: HeatBoundaryTerms, *,
                 weight=None, q_source=None,
                 q_update: Callable | None = None,
                 q_damping=None,
                 linear_solver: str = "sparsecholesky",
                 newton_tol_K: float = 1.0e-3, max_newton: int = 25,
                 max_halvings: int = 6):
        from ngsolve import CF, GridFunction

        self.gfT = gfT
        self.fes = gfT.space
        self.mesh = self.fes.mesh
        self.material = material
        self.bnd = boundaries
        self.w = CF(1.0) if weight is None else weight
        from ngsolve import Parameter
        # q_source is a coefficient function (a GridFunction when q depends
        # on T); q_update(gfT) refreshes it in place at every iterate.
        self.q_source = CF(0.0) if q_source is None else q_source
        self.q_update = q_update
        # optional non-negative -dq/dT on the heat-flux boundary (Newton
        # term of a temperature-dependent source)
        self.q_damping = q_damping
        self._rdt = Parameter(1.0)
        self._forms = None
        self.linear_solver = linear_solver
        self.tol = float(newton_tol_K)
        self.max_newton = int(max_newton)
        self.max_halvings = int(max_halvings)
        self.gf_old = GridFunction(self.fes)
        self.intorder = 2 * int(self.fes.globalorder) + 3
        self._dx, self._ds = self._measures()
        self.audit = TransientAudit()
        from ngsolve import NodeId, VERTEX
        self._vdofs = np.asarray([self.fes.GetDofNrs(NodeId(VERTEX, v.nr))[0]
                                  for v in self.mesh.vertices])
        from ngsolve import IfPos
        self.k_cf, self.H_cf, self.c_cf = material.coefficient_functions(gfT)
        self.k_old, self.H_old, _ = material.coefficient_functions(self.gf_old)
        dT = gfT - self.gf_old
        eps = 1.0e-6
        big = IfPos(dT * dT - eps * eps, 1.0, 0.0)
        dT_safe = IfPos(dT - eps, dT, IfPos(-dT - eps, dT, eps))
        self.c_secant = (big * (self.H_cf - self.H_old) / dT_safe
                         + (1.0 - big) * self.c_cf)

    # -- helpers ---------------------------------------------------------
    def _measures(self):
        from ngsolve import BND, VOL, IntegrationRule, dx, ds
        rules = {}
        for vb in (VOL, BND):
            types = {el.type for el in self.mesh.Elements(vb)}
            rules[vb] = {et: IntegrationRule(et, self.intorder)
                         for et in types}
        self._rules = rules
        return (lambda **kw: dx(intrules=rules[VOL], **kw),
                lambda region: ds(region, intrules=rules[BND]))

    def _nodal(self, gf=None):
        return np.asarray((gf or self.gfT).vec.FV().NumPy())[self._vdofs]

    def _check_range(self):
        T = self._nodal()
        lo, hi = self.material.T[0], self.material.T[-1]
        over = max(0.0, float(T.max() - hi), float(lo - T.min()))
        if over > 0.0:
            if not self.material.allow_extrapolation:
                raise ValueError(
                    f"temperature {T.min():.1f}..{T.max():.1f} C left the "
                    f"material table range {lo:.1f}..{hi:.1f} C; extend the "
                    "table or allow extrapolation explicitly")
            self.audit.table_extrapolation_C = max(
                self.audit.table_extrapolation_C, over)

    def energy(self, gf=None):
        from ngsolve import Integrate
        H = self.H_cf if gf is None else self.material.coefficient_functions(
            gf)[1]
        return float(Integrate(H * self.w, self.mesh,
                               order=self.intorder).real)

    def _flows(self, q_cf):
        """(heat input, losses) [W] at the current iterate."""
        from ngsolve import BND, Integrate
        b = self.bnd
        qin = 0.0
        if b.heat_flux:
            qin = float(Integrate(q_cf * self.w, self.mesh, BND,
                                  definedon=self.mesh.Boundaries(
                                      b.heat_flux),
                                  order=self.intorder).real)
        loss = 0.0
        if b.h_conv and b.convection:
            loss += float(Integrate(b.h_conv * (self.gfT - b.t_ext) * self.w,
                                    self.mesh, BND, definedon=self.mesh
                                    .Boundaries(b.convection),
                                    order=self.intorder).real)
        if b.emissivity and b.radiation:
            TK = self.gfT + KELVIN
            loss += float(Integrate(
                b.emissivity * SIGMA_SB * (TK ** 4 - (b.t_ext + KELVIN) ** 4)
                * self.w, self.mesh, BND,
                definedon=self.mesh.Boundaries(b.radiation),
                order=self.intorder).real)
        return qin, loss

    # -- one Newton-solved step -------------------------------------------
    def _build_forms(self):
        """Residual and matrix, built once; they re-read gfT on assembly."""
        from ngsolve import BilinearForm, InnerProduct, LinearForm, grad
        b = self.bnd
        u, v = self.fes.TnT()
        DX, DS = self._dx(), self._ds
        rdt = self._rdt
        R = LinearForm(self.fes)
        R += (self.H_cf - self.H_old) * rdt * v * self.w * DX
        R += self.k_cf * InnerProduct(grad(self.gfT), grad(v)) * self.w * DX
        J = BilinearForm(self.fes, symmetric=True)
        J += self.c_secant * rdt * u * v * self.w * DX
        J += self.k_cf * InnerProduct(grad(u), grad(v)) * self.w * DX
        if b.heat_flux:
            R += -self.q_source * v * self.w * DS(b.heat_flux)
            if self.q_damping is not None:
                J += self.q_damping * u * v * self.w * DS(b.heat_flux)
        if b.h_conv and b.convection:
            R += b.h_conv * (self.gfT - b.t_ext) * v * self.w \
                * DS(b.convection)
            J += b.h_conv * u * v * self.w * DS(b.convection)
        if b.emissivity and b.radiation:
            TK = self.gfT + KELVIN
            R += b.emissivity * SIGMA_SB * (TK ** 4 - (b.t_ext + KELVIN)
                                            ** 4) * v * self.w \
                * DS(b.radiation)
            J += 4.0 * b.emissivity * SIGMA_SB * TK ** 3 * u * v \
                * self.w * DS(b.radiation)
        return R, J

    def _solve_step(self, dt):
        from ngsolve import TaskManager
        if self._forms is None:
            self._forms = self._build_forms()
        R, J = self._forms
        self._rdt.Set(1.0 / dt)
        self.gf_old.vec.data = self.gfT.vec
        du = self.gfT.vec.CreateVector()
        for it in range(1, self.max_newton + 1):
            if self.q_update is not None:
                self.q_update(self.gfT)
            with TaskManager():
                R.Assemble()
                J.Assemble()
                inv = J.mat.Inverse(self.fes.FreeDofs(),
                                    inverse=self.linear_solver)
                du.data = inv * R.vec
                self.gfT.vec.data -= du
            step = float(np.max(np.abs(np.asarray(du.FV().NumPy())
                                       [self._vdofs])))
            if not math.isfinite(step):
                break
            if step < self.tol:
                if self.q_update is not None:
                    self.q_update(self.gfT)
                return it, self.q_source
        self.gfT.vec.data = self.gf_old.vec
        if self.q_update is not None:
            self.q_update(self.gfT)
        return None, None

    def advance(self, dt, depth=0):
        """Advance by ``dt``, halving on non-convergence."""
        if depth == 0:
            self.audit.steps += 1
        e0 = self.energy()
        its, q_cf = self._solve_step(dt)
        if its is None:
            if depth >= self.max_halvings:
                raise RuntimeError(
                    f"the nonlinear heat step did not converge within "
                    f"{self.max_newton} iterations even after {depth} "
                    f"halvings (dt={dt:.3e} s); reduce --dt or check the "
                    "material table")
            self.audit.halvings += 1
            self.advance(0.5 * dt, depth + 1)
            self.advance(0.5 * dt, depth + 1)
            return
        self._check_range()
        qin, loss = self._flows(q_cf)
        self.audit.substeps += 1
        self.audit.newton_iterations.append(its)
        self.audit.max_newton_iterations = max(
            self.audit.max_newton_iterations, its)
        self.audit.energy_in_J += qin * dt
        self.audit.energy_loss_J += loss * dt
        self.audit.energy_stored_J += self.energy() - e0
        self.last_q_cf = q_cf
        self.last_heat_input_W = qin
