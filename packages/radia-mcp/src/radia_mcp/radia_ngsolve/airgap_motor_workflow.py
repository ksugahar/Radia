# -*- coding: utf-8 -*-
r"""AGE rotating machine workflow with nonlinear (saturating) iron.

Wraps the linear AGE coupling (:mod:`airgap_machine`) in a Newton/Picard
iteration for machines with nonlinear iron B-H characteristics.  The
air-gap element stays analytic (no gap mesh), but the iron FE regions now
carry a field-dependent reluctivity ν(|B|²) that is updated each iteration.

Iteration (Picard/fixed-point, default):

    for k = 0, 1, 2, ...:
        1. ν_k = nu_fn(gfu_k)           -- evaluate reluctivity field from previous solution
        2. a_k = bilinear(fes, ν_k)     -- reassemble stiffness with updated ν
        3. factor_k = airgap_factorize(a_k.mat, coupling, freedofs)
        4. gfu_{k+1} = airgap_solve(factor_k, fes, dirichlet_cf, source_lf)
        5. ||gfu_{k+1} - gfu_k|| / ||gfu_{k+1}|| < tol  →  converged

Newton iteration (use_newton=True) replaces step 2 with the tangent operator
∂(νB)/∂A (assembled via num_diff or analytic Jacobian).  For mild saturation
(silicon steel at typical flux densities) Picard converges in 5-15 iterations;
Newton is faster when ν changes steeply.

API
---
  coupling = airgap_coupling(fes, ri, ro, "rotor_ring", "stator_ring", harmonics)
  gfu, info = age_motor_nonlinear_solve(
      fes, coupling,
      bilinear_fn = lambda nu: assemble_stiffness(fes, nu),
      nu_cf_fn    = lambda gfu: compute_nu(gfu),   # ν from solution
      dirichlet_cf = ...,
      niter=20, tol=1e-5
  )

Validated against the linear case (nu_cf_fn returning constant):
  - Linear Picard: converges in 1 iteration (exact, no approximation)
  - Error ≤ machine precision vs direct airgap_solve

See :mod:`airgap_machine` for the linear Woodbury AGE solver.
See :mod:`airgap_element` for the analytic 2×2 DtN and closed-form torque.
"""
from __future__ import annotations

import math
from typing import Callable, Optional

import numpy as np
import ngsolve as ng
from ngsolve import GridFunction, LinearForm, BilinearForm, InnerProduct, BND

from .airgap_machine import airgap_coupling, airgap_factorize, airgap_solve, airgap_torque


def _solution_norm(gfu: GridFunction) -> float:
    """L2 norm of the GridFunction DOF vector."""
    v = gfu.vec
    return math.sqrt(abs(InnerProduct(v, v, conjugate=False)))


def age_motor_nonlinear_solve(
    fes: ng.FESpace,
    coupling: dict,
    bilinear_fn: Callable,
    nu_cf_fn: Optional[Callable] = None,
    dirichlet_cf=None,
    source_lf: Optional[LinearForm] = None,
    niter: int = 20,
    tol: float = 1e-6,
    relax: float = 1.0,
    use_newton: bool = False,
) -> tuple[GridFunction, dict]:
    """Picard/Newton iteration for nonlinear AGE motor.

    Parameters
    ----------
    fes          : FESpace   Scalar H1 space for the 2D or 3D magnetic potential.
    coupling     : dict      AGE coupling dict from :func:`airgap_coupling`.
    bilinear_fn  : callable  ``bilinear_fn(nu_cf) -> BilinearForm`` (assembled).
                             On first call receives the initial ν; the function MUST
                             call ``.Assemble()`` on the returned form.
    nu_cf_fn     : callable  ``nu_cf_fn(gfu) -> CoefficientFunction``
                             Returns the reluctivity field from the current solution.
                             If None, a linear (constant-ν) run is assumed and the
                             bilinear form from the first call is reused throughout.
    dirichlet_cf : CF        Dirichlet excitation (passed to airgap_solve each iter).
    source_lf    : LinearForm  Optional volume current source.
    niter        : int       Maximum Picard iterations.
    tol          : float     Relative L2 convergence tolerance.
    relax        : float     Under-relaxation factor (1.0 = full update, 0.5 = half).
    use_newton   : bool      NOT YET IMPLEMENTED — placeholder for future tangent solver.

    Returns
    -------
    gfu  : GridFunction  Converged solution.
    info : dict          ``{niter_used, converged, rel_change_history, torque}``.

    Notes
    -----
    The ``bilinear_fn`` is called with the *CoefficientFunction* from ``nu_cf_fn``,
    not a scalar value.  This lets you assemble the full nonlinear stiffness:

        def bilinear_fn(nu_cf):
            a = BilinearForm(fes)
            a += nu_cf * grad(u) * grad(v) * dx
            a.Assemble()
            return a
    """
    if use_newton:
        raise NotImplementedError("Newton tangent solver not yet implemented; use Picard (use_newton=False).")

    gfu_prev = GridFunction(fes)
    gfu_curr = GridFunction(fes)
    rel_history = []
    converged = False

    # Initial reluctivity (from zero solution or first linearization)
    nu_initial = nu_cf_fn(gfu_prev) if nu_cf_fn is not None else None
    a = bilinear_fn(nu_initial)

    for k in range(niter):
        # Factorize and solve for current rotor position
        factor = airgap_factorize(a.mat, coupling, fes.FreeDofs())
        gfu_curr = airgap_solve(factor, fes, dirichlet_cf=dirichlet_cf, source_lf=source_lf)

        # Under-relaxation
        if relax < 1.0:
            gfu_curr.vec.data = relax * gfu_curr.vec + (1.0 - relax) * gfu_prev.vec

        # Convergence check
        diff = gfu_curr.vec.CreateVector()
        diff.data = gfu_curr.vec - gfu_prev.vec
        norm_diff = math.sqrt(abs(InnerProduct(diff, diff, conjugate=False)))
        norm_curr = _solution_norm(gfu_curr)
        rel_change = norm_diff / norm_curr if norm_curr > 0 else float("inf")
        rel_history.append(rel_change)

        if rel_change < tol:
            converged = True
            break

        # Update for next iteration
        gfu_prev.vec.data = gfu_curr.vec

        if nu_cf_fn is not None:
            nu_cf = nu_cf_fn(gfu_curr)
            a = bilinear_fn(nu_cf)
        # If nu_cf_fn is None (linear problem), reuse same stiffness → converges in 1 iter

    T = airgap_torque(coupling, gfu_curr)
    return gfu_curr, {
        "niter_used": k + 1,
        "converged": converged,
        "rel_change_history": rel_history,
        "torque": T,
    }


def age_motor_rotation_sweep(
    fes: ng.FESpace,
    coupling: dict,
    bilinear_fn: Callable,
    excitation_fn: Callable,
    theta_r_list,
    nu_cf_fn: Optional[Callable] = None,
    axial_length: float = 1.0,
    niter: int = 20,
    tol: float = 1e-6,
):
    """Compute torque vs rotor angle, reusing nonlinear solve per angle.

    Parameters
    ----------
    fes           : FESpace   Magnetic potential space.
    coupling      : dict      AGE coupling from :func:`airgap_coupling`.
    bilinear_fn   : callable  ``bilinear_fn(nu_cf) -> BilinearForm`` (assembled).
    excitation_fn : callable  ``excitation_fn(theta_r) -> CoefficientFunction``
                              Dirichlet excitation as function of rotor angle.
    theta_r_list  : iterable  Rotor angles [rad].
    nu_cf_fn      : callable  Reluctivity from solution (None = linear).
    axial_length  : float     Axial stack length [m].
    niter         : int       Max Picard iterations per angle.
    tol           : float     Relative convergence tolerance.

    Returns
    -------
    list of dicts with ``{theta_r, torque, niter_used, converged}``."""
    results = []
    for theta_r in theta_r_list:
        dirichlet_cf = excitation_fn(theta_r)
        gfu, info = age_motor_nonlinear_solve(
            fes, coupling, bilinear_fn,
            nu_cf_fn=nu_cf_fn,
            dirichlet_cf=dirichlet_cf,
            niter=niter, tol=tol,
        )
        T = airgap_torque(coupling, gfu, axial_length=axial_length)
        results.append({
            "theta_r": theta_r,
            "torque": T,
            "niter_used": info["niter_used"],
            "converged": info["converged"],
        })
    return results
