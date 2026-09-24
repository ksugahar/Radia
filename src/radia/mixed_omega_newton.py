"""Quadrature-based Newton iteration for total/reduced scalar potentials.

The source split, interface constraint and finite-element spaces are fixed.
Only the constitutive tangent is reassembled; the PCHIP coenergy is
differentiated by NGSolve, without a projected permeability iteration.
"""
from __future__ import annotations

import math
import time
import numpy as np


class MixedOmegaNewtonNotConverged(RuntimeError):
    def __init__(self, message, state):
        super().__init__(message)
        self.state = state


def _memoize_linked_source(function):
    from .kelvin_solver import memoize_linked_source
    return memoize_linked_source(function)


@_memoize_linked_source
def solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin(
        mesh, H_s, source_potential, R_K, offset, *, bh_table,
        nonlinear_materials, reduced_materials, total_materials,
        interface_boundary, order=1, mu_r_initial=1000.0,
        tolerance=2e-5, residual_tolerance=1e-8, max_iterations=40,
        max_halvings=12, observation_points=None, progress_callback=None,
        bonus_intorder=4, inverse="pardiso", mu_r_by_material=None,
        condense_matching_trace=False, linear_solver="direct",
        **linear_options):
    """Newton with residual backtracking and the production PCHIP B(H) law.

    Nonlinear materials must lie in the physical total-potential region.
    Orders one and two use the same quadrature-evaluated constitutive law.
    Caller owns TaskManager, as for the linear mixed solver. Failure raises
    ``MixedOmegaNewtonNotConverged`` carrying iteration diagnostics.
    """
    import ngsolve as ng
    from .kelvin_solver import (
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
        solve_magnetostatic_matching_trace_total_reduced_omega,
        audit_mixed_omega_constitutive_field)
    from .scalar_potential_solver import (
        _build_bh_coefficient_function, _build_bh_coenergy_coefficient_function)

    mu0 = 4e-7 * math.pi
    table = np.asarray(bh_table, dtype=float)
    if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
            or not np.isfinite(table).all() or np.any(table[0] != 0)
            or np.any(np.diff(table[:, 0]) <= 0)
            or np.any(np.diff(table[:, 1]) <= 0)):
        raise ValueError("B(H) must start at [0,0] and have finite strictly increasing H and B")
    nonlinear = tuple(nonlinear_materials)
    reduced = tuple(reduced_materials)
    total = tuple(total_materials)
    if not nonlinear or not set(nonlinear) <= set(total):
        raise ValueError("Newton nonlinear_materials must be nonempty and in total_materials")
    if int(order) not in (1, 2):
        raise ValueError("Newton supports order 1 or 2")
    if (not math.isfinite(float(mu_r_initial)) or float(mu_r_initial) <= 0
            or not 0 < tolerance < 1 or not 0 < residual_tolerance < 1
            or int(max_iterations) < 1 or int(max_halvings) < 0):
        raise ValueError("invalid Newton initial permeability, tolerance or iteration limit")
    if progress_callback is not None and not callable(progress_callback):
        raise ValueError("progress_callback must be callable")
    if "return_system" in linear_options or "mu_cf" in linear_options:
        raise ValueError("Newton owns the assembled system and constitutive coefficient")
    if linear_solver not in ("direct", "cg") or (linear_solver == "cg" and not condense_matching_trace):
        raise ValueError("Newton CG requires matching-trace condensation")

    initial_mu = dict(mu_r_by_material or {})
    initial_mu.update({name: float(mu_r_initial) for name in nonlinear})
    started = time.perf_counter()
    if condense_matching_trace:
        allowed = {"kelvin_mats", "kelvin_match_exact", "dirichlet_bbbnd",
                   "total_source_h", "total_source_materials", "interface_constraint_scale"}
        active_options = {key for key, value in linear_options.items() if value is not None}
        if linear_options.get("kelvin_mats", ()) or active_options - allowed:
            raise ValueError("matching-trace Newton supports finite domains with a point gauge only")
        result = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, H_s, source_potential, reduced_materials=reduced, total_materials=total,
            interface_boundary=interface_boundary, order=order, mu_r_by_material=initial_mu,
            bonus_intorder=bonus_intorder, inverse=inverse, return_system=True,
            solver=linear_solver, cg_tolerance=1e-10,
            dirichlet_bbbnd=linear_options.get("dirichlet_bbbnd", "GND"),
            total_source_h=linear_options.get("total_source_h"),
            total_source_materials=linear_options.get("total_source_materials", ()))
        result.update(kelvin_materials=(), total_source_h=linear_options.get("total_source_h"),
                      total_source_materials=linear_options.get("total_source_materials", ()))
    else:
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, H_s, source_potential, R_K, offset,
            reduced_materials=reduced, total_materials=total,
            interface_boundary=interface_boundary, order=order,
            mu_r_by_material=initial_mu, bonus_intorder=bonus_intorder,
            inverse=inverse, return_system=True, **linear_options)
    initial_seconds = time.perf_counter() - started
    fes = result["fes"]
    solution = result["solution"]
    rhs = result["system"]["linear_form"].vec.CreateVector()
    rhs.data = result["system"]["linear_form"].vec
    free = np.fromiter(fes.FreeDofs(), dtype=bool, count=fes.ndof)
    selector = mesh.Materials("|".join(nonlinear))
    kelvin = tuple(result["kelvin_materials"])
    if set(nonlinear) & set(kelvin):
        raise ValueError("Kelvin exterior materials cannot be nonlinear")
    core = tuple(name for name in total if name not in set(kelvin))
    tangent = ng.BilinearForm(fes, symmetric=True)
    if condense_matching_trace:
        ut, vt = fes.TnT()
        tangent += result["mu_cf"] * ng.grad(ut) * ng.grad(vt) * ng.dx(bonus_intorder=bonus_intorder)
    else:
        (ur, ut, lm), (vr, vt, vm) = fes.TnT()
        tangent += result["mu_cf"] * ng.grad(ur) * ng.grad(vr) * ng.dx(
            definedon=mesh.Materials("|".join(reduced + kelvin)), bonus_intorder=bonus_intorder)
        tangent += result["mu_cf"] * ng.grad(ut) * ng.grad(vt) * ng.dx(
            definedon=mesh.Materials("|".join(core)), bonus_intorder=bonus_intorder)
        scale = linear_options.get("interface_constraint_scale")
        if scale is None:
            scale = mu0 / float(R_K)
        tangent += scale * (lm * (vt.Trace() - vr.Trace())
                            + vm * (ut.Trace() - ur.Trace())) * ng.ds(
            definedon=mesh.Boundaries(interface_boundary), bonus_intorder=bonus_intorder)
    harmonic = result["total_source_h"]
    harmonic_names = set(result["total_source_materials"])
    zero = ng.CF((0., 0., 0.))
    source = mesh.MaterialCF({name: harmonic if harmonic is not None
                            and name in harmonic_names else zero for name in core})
    trial_h = source - ng.grad(ut)
    # Smooth only the radial norm at machine-scale H, not the B(H) table.
    hmag = ng.sqrt(ng.InnerProduct(trial_h, trial_h) + 1e-24)
    correction = (_build_bh_coenergy_coefficient_function(hmag, table)
                  - 0.5 * mu0 * float(mu_r_initial) * ng.InnerProduct(trial_h, trial_h))
    tangent += ng.SymbolicEnergy(correction.Compile(), definedon=selector,
                                bonus_intorder=bonus_intorder)
    preconditioner = ng.Preconditioner(tangent, "local") if linear_solver == "cg" else None
    residual = rhs.CreateVector()
    trial = solution.vec.CreateVector()
    step = solution.vec.CreateVector()

    def evaluate(vector):
        tangent.Apply(vector, residual)
        residual.data -= rhs
        value = float(np.linalg.norm(residual.FV().NumPy()[free]))
        if not math.isfinite(value):
            raise RuntimeError("non-finite mixed Omega Newton residual")
        return value

    norm = evaluate(solution.vec)
    reference_norm = max(float(np.linalg.norm(rhs.FV().NumPy()[free])), norm, 1e-30)
    history = []
    converged = norm / reference_norm <= residual_tolerance
    relative_change = 0.0 if converged else None
    integration_order = max(4, 2 * int(order) + int(bonus_intorder))
    # H_cf refers to this persistent solution; updating it preserves all lifts.
    field_h = result["H_cf"]
    magnitude = ng.sqrt(ng.InnerProduct(field_h, field_h) + 1e-24)
    secant = (_build_bh_coefficient_function(magnitude, table) / magnitude).Compile()
    physical_mu = mesh.MaterialCF({name: secant for name in nonlinear}, default=result["mu_cf"])
    field_b = (physical_mu * field_h).Compile()
    previous = ng.GridFunction(fes)
    reason = "iteration limit"
    for iteration in range(1, int(max_iterations) + 1):
        if converged:
            break
        row = {"iteration": iteration, "residual_relative_before": norm / reference_norm}
        previous.vec.data = solution.vec
        before_h = source - ng.grad(previous if condense_matching_trace else previous.components[1])
        before_mag = ng.sqrt(ng.InnerProduct(before_h, before_h) + 1e-24)
        before_b = (_build_bh_coefficient_function(before_mag, table) / before_mag * before_h).Compile()
        t0 = time.perf_counter()
        tangent.AssembleLinearization(solution.vec)
        row["assembly_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        if linear_solver == "cg":
            from ngsolve.krylovspace import CGSolver
            preconditioner.Update()
            inv = CGSolver(mat=tangent.mat, pre=preconditioner.mat, tol=1e-10,
                           maxiter=3000, printrates=False)
        else:
            inv = tangent.mat.Inverse(fes.FreeDofs(), inverse=inverse)
        row["factorization_s"] = time.perf_counter() - t0
        linear_rhs = residual.CreateVector()
        linear_rhs.data = residual
        linear_rhs.FV().NumPy()[~free] = 0.0
        step.data = -inv * linear_rhs
        linear_error = residual.CreateVector()
        linear_error.data = tangent.mat * step + residual
        absolute_linear = float(np.linalg.norm(linear_error.FV().NumPy()[free]))
        relative_linear = float(absolute_linear
                                / max(np.linalg.norm(residual.FV().NumPy()[free]), 1e-30))
        row["linear_residual_relative"] = relative_linear
        row["linear_iterations"] = getattr(inv, "iterations", None)
        # Near nonlinear convergence, roundoff divided by the tiny Newton RHS
        # is not a meaningful failure. Also require the error to exceed the
        # original equation's absolute convergence budget.
        if (not math.isfinite(relative_linear) or
                (relative_linear > 1e-6 and absolute_linear > residual_tolerance * reference_norm)):
            raise MixedOmegaNewtonNotConverged("Newton linear solve did not converge",
                                              {"nonlinear_stats": {"converged": False,
                                               "history": history + [row]}})
        del inv
        alpha = 1.0
        accepted = False
        t0 = time.perf_counter()
        for halving in range(int(max_halvings) + 1):
            trial.data = solution.vec + alpha * step
            trial_norm = evaluate(trial)
            if trial_norm <= (1 - 1e-4 * alpha) * norm or trial_norm / reference_norm <= residual_tolerance:
                accepted = True
                break
            alpha *= 0.5
        row.update(line_search_s=time.perf_counter() - t0, step_length=alpha,
                   line_search_evaluations=halving + 1, line_search_accepted=accepted)
        if not accepted:
            history.append(row)
            reason = "residual backtracking failed"
            break
        solution.vec.data = trial
        norm = trial_norm
        delta = (field_b - before_b).Compile()
        numerator = float(ng.Integrate(ng.InnerProduct(delta, delta), mesh,
                                      definedon=selector, order=integration_order))
        denominator = float(ng.Integrate(ng.InnerProduct(field_b, field_b), mesh,
                                        definedon=selector, order=integration_order))
        relative_change = math.sqrt(max(numerator, 0.) / max(denominator, 1e-60))
        row.update(residual_relative=norm / reference_norm,
                   relative_B_change=relative_change)
        history.append(row)
        if progress_callback is not None:
            progress_callback(dict(row))
        converged = norm / reference_norm <= residual_tolerance and relative_change <= tolerance

    stats = dict(method="quadrature_pchip_newton", material_sampling="integration_point",
                 matching_trace_condensed=bool(condense_matching_trace), linear_solver=linear_solver,
                 bh_interpolation="pchip", converged=bool(converged), iterations=len(history),
                 residual_relative=norm / reference_norm, residual_tolerance=float(residual_tolerance),
                 relative_B_change=relative_change, tolerance=float(tolerance), history=history,
                 initial_linear_solve_s=initial_seconds, elapsed_s=time.perf_counter() - started)
    if not converged:
        raise MixedOmegaNewtonNotConverged(f"mixed Omega Newton: {reason}", {"nonlinear_stats": stats})
    result.update(mu_cf=physical_mu, B_cf=field_b, nonlinear_stats=stats, system=None,
                  assembled_energy=None, linear_residual=None)
    result["constitutive_field_audit"] = audit_mixed_omega_constitutive_field(
        mesh, field_h, field_b, table, nonlinear, integration_order=integration_order)
    if observation_points is not None:
        points = np.asarray(observation_points, dtype=float).reshape(-1, 3)
        stats["observation_field_T"] = [list(field_b(mesh(*p))) for p in points]
    return result
