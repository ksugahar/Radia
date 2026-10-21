"""Bounded depth refinement for the ESIM phasor cell (not a cycle solver)."""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded


def depth_mesh(length, delta, count, cylinder=False):
    """Fixed physical grading strength; increasing count never collapses cells."""
    if not np.isfinite(length) or length <= 0 or int(count) != count or count < 4:
        raise ValueError('ESIM depth requires positive length and >=4 nodes')
    strength = np.log1p(length / max(float(delta), length * 1e-12))
    t = np.linspace(0., 1., int(count))
    depth = length * np.expm1(strength * t) / np.expm1(strength) if strength else length*t
    return length - depth[::-1] if cylinder else depth


def linear_cell(solver, H0, mu):
    """Conservative nonuniform control volumes, row scaled before band solve.

    exp(+i omega t): -rho div(grad H) + i omega mu H = 0.
    Axis and insulated midplane have zero flux; halfspace end has H=0.
    """
    x = solver.mesh_points
    dx = np.diff(x)
    if not np.all(np.isfinite(dx)) or np.any(dx <= 64*np.finfo(float).eps*x[-1]):
        raise RuntimeError('ESIM depth mesh has unresolved floating-point intervals')
    n = len(x)
    cylinder = getattr(solver, 'geometry', '') == 'cylinder'
    finite = hasattr(solver, 'half_thickness')
    faces = (x[:-1] + x[1:])*.5
    left = np.r_[x[0], faces]
    right = np.r_[faces, x[-1]]
    volume = (right**2-left**2)*.5 if cylinder else right-left
    conductance = solver.rho * (faces if cylinder else 1.) / dx
    lower = np.zeros(n-1, complex)
    upper = np.zeros(n-1, complex)
    diagonal = 1j*solver.omega*mu*volume
    diagonal[:-1] += conductance
    diagonal[1:] += conductance
    lower[:] = -conductance
    upper[:] = -conductance
    rhs = np.zeros(n, complex)
    if cylinder:
        diagonal[-1], lower[-1], rhs[-1] = 1., 0., H0
    else:
        diagonal[0], upper[0], rhs[0] = 1., 0., H0
        if not finite:
            diagonal[-1], lower[-1], rhs[-1] = 1., 0., 0.
    scale = np.abs(diagonal)
    lower /= scale[1:]
    upper /= scale[:-1]
    diagonal /= scale
    rhs /= scale
    band = np.zeros((3, n), complex)
    band[0, 1:], band[1], band[2, :-1] = upper, diagonal, lower
    field = solve_banded((1, 1), band, rhs)
    residual = diagonal*field-rhs
    residual[1:] += lower*field[:-1]
    residual[:-1] += upper*field[1:]
    denominator = np.linalg.norm(rhs) + np.linalg.norm(diagonal*field)
    relative = np.linalg.norm(residual)/max(denominator, np.finfo(float).tiny)
    if not np.all(np.isfinite(field)) or relative > 1e-6:
        raise RuntimeError(f'ESIM algebraic residual failed: {relative:g}')
    solver._last_linear_residual = float(relative)
    return field


class AdaptiveDepth:
    def _configure_depth(self, mesh_tol, max_nodes, max_refinements, cell_max_iter):
        if not np.isfinite(self.sigma) or self.sigma <= 0 or not np.isfinite(self.frequency) or self.frequency <= 0:
            raise ValueError('ESIM diffusion requires positive finite conductivity and frequency')
        if not np.isfinite(mesh_tol) or not 0 < mesh_tol < 1:
            raise ValueError('mesh_tol must be in (0,1)')
        for name, value, minimum in [('max_nodes', max_nodes, 8),
                                     ('max_refinements', max_refinements, 1),
                                     ('cell_max_iter', cell_max_iter, 1)]:
            if isinstance(value, bool) or int(value) != value or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}')
        self.mesh_tol = float(mesh_tol)
        self.max_nodes = int(max_nodes)
        self.max_refinements = int(max_refinements)
        self.cell_max_iter = int(cell_max_iter)

    def solve(self, H0, tol=1e-6, max_iter=None, relaxation=.5, R_local=None):
        result = self._adaptive_solve(H0, tol, max_iter, relaxation, R_local)
        if hasattr(self, 'half_thickness'):
            return result
        # A half-space needs a second, independent domain-depth check.
        original, original_L = self.mesh_points, self.L
        depth_history = []
        try:
            for _ in range(4):
                self.L *= 2
                self.mesh_points = depth_mesh(self.L, self.delta_initial, self.n_nodes)
                extended = self._adaptive_solve(H0, tol, max_iter, relaxation)
                error = max(abs(extended[k]-result[k])/max(abs(extended[k]), np.finfo(float).tiny)
                            for k in ('Z', 'P_prime', 'Q_prime'))
                depth_history.append(dict(length_m=self.L, relative_error=float(error)))
                if error <= self.mesh_tol:
                    extended.update(depth_converged=True, depth_history=depth_history)
                    return extended
                result = extended
            raise RuntimeError('ESIM half-space depth budget exhausted')
        finally:
            self.mesh_points, self.L = original, original_L

    def _adaptive_solve(self, H0, tol=1e-6, max_iter=None, relaxation=.5, R_local=None):
        if not np.isfinite(H0) or H0 <= 0:
            raise ValueError('ESIM H0 must be finite and positive (peak A/m)')
        max_iter = self.cell_max_iter if max_iter is None else max_iter
        if int(max_iter) != max_iter or max_iter < 1 or not 0 < relaxation <= 1 or not 0 < tol < 1:
            raise ValueError('Invalid ESIM nonlinear iteration controls')
        if R_local is not None:
            if not hasattr(self, 'set_radius'):
                raise ValueError('R_local requires a finite cell')
            self.set_radius(float(R_local))
        original, original_n = self.mesh_points, self.n_nodes
        history = []
        try:
            coarse = self._solve_on_mesh(H0, tol, int(max_iter), relaxation)
            for level in range(self.max_refinements):
                if not coarse['converged']:
                    raise RuntimeError(f"ESIM nonlinear iteration limit: iterations={coarse['iterations']}, relative_change={coarse['relative_change']:g}")
                x = self.mesh_points.copy()
                if 2*len(x)-1 > self.max_nodes:
                    raise RuntimeError('ESIM depth node budget exhausted before mesh convergence')
                self.mesh_points = np.sort(np.r_[x, (x[:-1]+x[1:])*.5])
                self.n_nodes = len(self.mesh_points)
                fine = self._solve_on_mesh(H0, tol, int(max_iter), relaxation)
                if not fine['converged']:
                    raise RuntimeError(f"ESIM nonlinear iteration limit on refined mesh: iterations={fine['iterations']}")
                errors = [abs(fine[k]-coarse[k])/max(abs(fine[k]), np.finfo(float).tiny)
                          for k in ('Z', 'P_prime', 'Q_prime')]
                history.append(dict(nodes=self.n_nodes, Z_real=float(fine['Z'].real),
                    Z_imag=float(fine['Z'].imag), P_prime=float(fine['P_prime']),
                    relative_errors=errors, iterations=fine['iterations']))
                if max(errors) <= self.mesh_tol:
                    fine.update(mesh_converged=True, mesh_relative_error=max(errors),
                                mesh_history=history, linear_residual=self._last_linear_residual)
                    return fine
                # Compare local Joule integrals on nested pairs; bulk-mark the
                # dominant depth intervals instead of refining the inert core.
                h = coarse['H_solution']
                hf = fine['H_solution']
                q = abs(np.diff(h))**2/np.diff(x)
                qf = (abs(np.diff(hf))**2/np.diff(self.mesh_points)).reshape(-1, 2).sum(axis=1)
                weight = (x[1:]+x[:-1])/(2*x[-1]) if getattr(self, 'geometry', '') == 'cylinder' else 1.
                indicator = abs(qf-q)*weight
                order = np.argsort(indicator)[::-1]
                total = indicator.sum()
                marked = order[:np.searchsorted(np.cumsum(indicator[order]), .75*total)+1] if total else order
                self.mesh_points = np.sort(np.r_[x, (x[marked]+x[marked+1])*.5])
                self.n_nodes = len(self.mesh_points)
                coarse = self._solve_on_mesh(H0, tol, int(max_iter), relaxation)
            raise RuntimeError('ESIM depth refinement budget exhausted before mesh convergence')
        finally:
            # Numerical output owns its mesh; cache identity remains the
            # declared seed/accuracy controls, not last field's adaptive grid.
            self.mesh_points, self.n_nodes = original, original_n
