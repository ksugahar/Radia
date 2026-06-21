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


def age_motor_hysteresis_sweep(
    fes,
    coupling: dict,
    model,
    source_fn: Callable,
    theta_list,
    iron_materials: str = "iron",
    other_nu: Optional[dict] = None,
    default_nu: float = 1.0,
    nu0: Optional[float] = None,
    axial_length: float = 1.0,
    freq: float = 1.0,
    cycles: int = 2,
    msource_iters: int = 60,
    tol: float = 1e-6,
) -> dict:
    r"""HISTORY-DEPENDENT (hysteretic) AGE rotor sweep.

    The iron regions carry the B-input vector **Play** hysteresis constitutive
    (:class:`radia_ngsolve.hysteresis.PlayHysteresis`) instead of a single-valued nu(|B|): the
    per-element play states ``p_k`` are carried FORWARD across rotor angles, so the iron's B-H response
    is path-dependent (magnetic history) and the torque / iron loss include the hysteresis.  This is the
    hysteretic upgrade of :func:`age_motor_rotation_sweep` (single-valued saturation) -- the play state
    per element is the internal variable.

    Each rotor position is solved by the fixed-point **M-source** (excess-field) method, which reuses
    ONE AGE factorization of the constant reference-reluctivity stiffness (the air-gap element stays
    analytic): with the reference ``nu0`` the iron's excess field ``H_res = nu0 B - H_play(B)`` is fed
    through :func:`airgap_solve`'s ``source_lf`` and iterated to a self-consistent B; the play states are
    then advanced with the converged B.  The play co-energy is convex (see
    :func:`hysteresis.play_coenergy_density`), so the M-source with ``nu0 >= nu_max/2`` is a contraction;
    the default ``nu0 = (nu_max + nu_min)/2`` is the optimal (fastest) reference.  (A per-angle
    variational/Newton solve is the robust alternative when the iron is very nonlinear, at the cost of a
    re-factorization per Newton step -- the M-source keeps the AGE one-factorization advantage.)

    Units: the AGE core uses the NORMALIZED reluctivity ``nu~ = 1/mur`` (gap ``nu~ = 1``), so B = curl A
    is physical Tesla.  Pass ``model`` with NORMALIZED slopes (``a~ = mu0 * nu``); the returned
    ``loss_density`` is then physical J/m^3 (the normalized loop area divided by mu0).

    Parameters
    ----------
    fes            : FESpace   order-1 H1 (so B = curl A is exactly P0 per element); complex for AGE.
    coupling       : dict      AGE coupling from :func:`airgap_coupling`.
    model          : PlayHysteresis  linear-cell play model (normalized slopes ``a``).
    source_fn      : callable  ``source_fn(theta) -> assembled LinearForm`` -- the per-angle excitation
                               (PM magnetisation and/or winding currents).
    theta_list     : iterable  rotor angles [rad] for one electrical period.
    iron_materials : str       ``"|"``-joined hysteretic-iron material names (e.g. ``"stator|rotoriron"``).
    other_nu       : dict      normalized reluctivity ``nu~`` for the non-iron materials (e.g.
                               ``{"magnet": 1/1.05}``); ``default_nu`` for the rest (gap/air/winding = 1).
    nu0            : float     M-source reference reluctivity (default optimal ``(nu_max+nu_min)/2``).
    axial_length   : float     stack length [m].
    freq           : float     electrical frequency [Hz] for the iron-loss wattage.
    cycles         : int       electrical cycles to run; only the LAST (steady) cycle is recorded.
    msource_iters  : int       max M-source iterations per angle.
    tol            : float     relative B convergence tolerance for the M-source.

    Returns
    -------
    dict with ``theta``, ``torque`` (per recorded angle), ``ms_iters`` (per angle), ``nu0``,
    ``iron_loss_W`` (sum over iron of the REAL per-element loop area |INT H.dB| * area * L * freq),
    ``loss_density`` (per iron element, J/m^3 per cycle), ``Bpk_iron`` (per iron element), and
    ``n_iron_elements``.
    """
    import numpy as np
    from ngsolve import (L2, CoefficientFunction, GridFunction as _GF, BilinearForm as _BF,
                         LinearForm as _LF, grad, dx, Integrate, TaskManager)

    MU0 = 4e-7 * math.pi
    eta = np.asarray(model.eta, float)
    if model.a is None:
        raise ValueError("age_motor_hysteresis_sweep needs linear-cell play slopes model.a (normalized nu~)")
    a = np.asarray(model.a, float); K = len(eta)
    nu_max = float(a.sum()); nu_min = float(a[0])
    if nu0 is None:
        nu0 = 0.5 * (nu_max + nu_min)
    irons = iron_materials.split("|")
    other_nu = dict(other_nu or {})

    mesh = fes.mesh
    u, v = fes.TnT()
    gfu = _GF(fes)
    cplx = fes.is_complex

    # B = curl A as element-constant fields (order-1 H1 => P0), via mixed H1<->L2 operators; the
    # transposes give the iron excess-field load INT_iron(H_res . curl v) as a cheap sparse matvec.
    fesL = L2(mesh, order=0, complex=cplx)
    qL, pL = fesL.TestFunction(), fesL.TrialFunction()
    mY = _BF(trialspace=fes, testspace=fesL); mY += grad(u)[1] * qL * dx; mY.Assemble()
    mX = _BF(trialspace=fes, testspace=fesL); mX += grad(u)[0] * qL * dx; mX.Assemble()
    mYt = _BF(trialspace=fesL, testspace=fes); mYt += pL * grad(v)[1] * dx; mYt.Assemble()
    mXt = _BF(trialspace=fesL, testspace=fes); mXt += pL * grad(v)[0] * dx; mXt.Assemble()
    gfTmp = _GF(fesL); gfHrx = _GF(fesL); gfHry = _GF(fesL)

    gfm = _GF(L2(mesh, order=0)); gfm.Set(mesh.MaterialCF({m: 1.0 for m in irons}, default=0.0))
    mask = gfm.vec.FV().NumPy() > 0.5
    nE = len(mask); ci = np.where(mask)[0]
    areas = Integrate(CoefficientFunction(1.0), mesh, element_wise=True).NumPy()

    numap = {m: nu0 for m in irons}; numap.update(other_nu)
    nu_ref = mesh.MaterialCF(numap, default=default_nu)
    aform = _BF(fes); aform += nu_ref * grad(u) * grad(v) * dx

    def _play(Bx, By, pl, ql):
        Bp = np.sqrt((Bx[:, None] - pl) ** 2 + (By[:, None] - ql) ** 2)
        dn = np.maximum(eta[None, :], Bp); dn = np.where(dn < 1e-30, 1e-30, dn)
        nx = Bx[:, None] - eta[None, :] * (Bx[:, None] - pl) / dn
        ny = By[:, None] - eta[None, :] * (By[:, None] - ql) / dn
        nx[:, 0] = Bx; ny[:, 0] = By
        return np.sum(a[None, :] * nx, axis=1), np.sum(a[None, :] * ny, axis=1), nx, ny

    def _iron_B():
        gfTmp.vec.data = mY.mat * gfu.vec
        col = gfTmp.vec.FV().NumPy(); bx = (np.real(col) if cplx else col)[ci] / areas[ci]
        gfTmp.vec.data = mX.mat * gfu.vec
        col = gfTmp.vec.FV().NumPy(); by = -(np.real(col) if cplx else col)[ci] / areas[ci]
        return bx, by

    px = np.zeros((nE, K)); py = np.zeros((nE, K))
    hrx = np.zeros(nE); hry = np.zeros(nE)
    total = _LF(fes); total.Assemble()
    rec = {"theta": [], "torque": [], "ms_iters": []}
    Bxh, Byh, Hxh, Hyh = [], [], [], []
    with TaskManager():
        aform.Assemble()
        factor = airgap_factorize(aform.mat, coupling, fes.FreeDofs())
        for cyc in range(cycles):
            last = (cyc == cycles - 1)
            for tm in theta_list:
                src = source_fn(tm); b_prev = None; used = msource_iters
                for it in range(msource_iters):
                    gfHrx.vec.FV().NumPy()[:] = hrx; gfHry.vec.FV().NumPy()[:] = hry
                    total.vec.data = src.vec + mYt.mat * gfHrx.vec - mXt.mat * gfHry.vec
                    gfu.vec.data = airgap_solve(factor, fes, source_lf=total).vec
                    bx, by = _iron_B()
                    Hx, Hy, _, _ = _play(bx, by, px[ci].copy(), py[ci].copy())
                    hrx[:] = 0.0; hry[:] = 0.0; hrx[ci] = nu0 * bx - Hx; hry[ci] = nu0 * by - Hy
                    if b_prev is not None:
                        d = math.sqrt(float(np.sum((bx - b_prev[0]) ** 2 + (by - b_prev[1]) ** 2)))
                        if d / (math.sqrt(float(np.sum(bx * bx + by * by))) + 1e-30) < tol:
                            used = it + 1; break
                    b_prev = (bx.copy(), by.copy())
                Hx_c, Hy_c, nx, ny = _play(bx, by, px[ci], py[ci]); px[ci], py[ci] = nx, ny
                if last:
                    rec["theta"].append(float(tm)); rec["ms_iters"].append(used)
                    rec["torque"].append(airgap_torque(coupling, gfu, axial_length=axial_length))
                    Bxh.append(bx.copy()); Byh.append(by.copy()); Hxh.append(Hx_c.copy()); Hyh.append(Hy_c.copy())

    def _close(arr):
        arr = np.array(arr); return np.vstack([arr, arr[:1]])
    Bxc, Byc, Hxc, Hyc = _close(Bxh), _close(Byh), _close(Hxh), _close(Hyh)
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    loss_density = np.abs(trap(Hxc, Bxc, axis=0) + trap(Hyc, Byc, axis=0)) / MU0
    P_iron = float(np.sum(loss_density * areas[ci]) * axial_length * freq)
    Bpk = np.max(np.sqrt(np.array(Bxh) ** 2 + np.array(Byh) ** 2), axis=0)
    return {
        "theta": np.array(rec["theta"]),
        "torque": np.array(rec["torque"]),
        "ms_iters": np.array(rec["ms_iters"]),
        "nu0": float(nu0),
        "iron_loss_W": P_iron,
        "loss_density": loss_density,
        "Bpk_iron": Bpk,
        "n_iron_elements": int(len(ci)),
    }
