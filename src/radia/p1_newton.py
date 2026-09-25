"""First-order (lowest-order Nedelec) magnetostatic Newton for reduced or total A.

``B = B_s + curl A`` (reduced, ``source_cf`` = B_s) or ``B = curl A`` with a
volume current (total, ``current_cf`` = J on ``current_materials``), homogeneous
tangential A on the Dirichlet boundary, a tabulated B-H law on the iron.

The discretisation is fixed by the element: curl of a lowest-order Nedelec
function is constant per element, so

* the flux density per iron element is ``curl A_e`` plus, for reduced A, the
  element mean of B_s projected once to discontinuous order ``p`` on the iron;
* the reluctivity is the exact inverse of the shared PCHIP law at that element
  flux (``_build_nu_of_b_interpolator``), sampled on a dense grid;
* one native pass over the elements (:class:`ElementCurl`) gives the element
  dofs, the element-constant curl ``C_e`` of the space's basis, the volumes,
  the element-graph matrix with the Jacobian's constant part (nu0 curl-curl off
  the iron) and every element's matrix positions; the residual and the element
  flux are gathers and scatters with ``C_e``, and the Newton Jacobian adds the
  closed-form iron element matrices ``vol_e C_e^T (nu_e I + q_e b_e b_e^T) C_e``
  at those positions (:class:`IronJacobian`).

Linear solves: compiled AMS + CG to a true relative residual (default,
``beta_zero`` on the ungauged compatible system, frozen AMG hierarchy across
Newton steps), compiled shifted IC(0) CG, or the SPD sparse direct solve.
Inner tolerances follow the nonlinear residual (inexact Newton) but never go
below ``0.1 * newton_tolerance * |R_0|``; Armijo backtracking on the residual.

TaskManager: call :func:`solve_p1_newton` **outside** ``ngsolve.TaskManager``.
The native AMS setup refuses to run inside one (it would mutate its hierarchy
in a parallel region), so this solver owns its parallel regions and raises if
it is entered inside a TaskManager.
"""
from __future__ import annotations

import math
import re
import time

import numpy as np
import ngsolve as ng
from ngsolve.krylovspace import CGSolver

MU0 = 4.0e-7 * math.pi
NU0 = 1.0 / MU0


def _sparsesolv():
    import radia.sparsesolv_ngsolve as ssn

    for name in ("LowestOrderGradient", "LowestOrderCurlSystem", "LowestOrderCurlJacobian",
                 "HypreBasedAMSPreconditioner",
                 "TaskManagerActive"):
        if not hasattr(ssn, name):
            raise RuntimeError(f"radia.sparsesolv_ngsolve lacks {name}; rebuild the native module")
    return ssn


def lowest_order_gradient(fes):
    """Discrete gradient H1(order 1) -> lowest-order HCurl, from the edge table.

    Equals ``fes.CreateGradient()[0]`` (bit for bit) and is built in parallel
    in C++; raises for a space that is not one dof per edge in edge order.
    """
    return _sparsesolv().LowestOrderGradient(fes)


class FluxSideLaw:
    """The shared monotone PCHIP B(H) law seen from the flux-density side.

    ``nu(|B|)`` is the exact inverse of the production PCHIP law sampled on a
    dense geometric grid and interpolated linearly; ``dH/dB`` is the numerical
    derivative of the sampled ``H = nu B``.
    """

    def __init__(self, bh_table, *, samples: int = 4001):
        from radia.vector_potential_solver import _build_nu_of_b_interpolator

        nu_of_b, b_saturation = _build_nu_of_b_interpolator(bh_table)
        self.b_saturation = float(b_saturation)
        self.nu_initial = float(nu_of_b(0.0))
        grid = np.concatenate(([0.0], np.geomspace(1.0e-7, 20.0 * self.b_saturation, int(samples))))
        nu = np.asarray([nu_of_b(value) for value in grid], dtype=float)
        dhdb = np.gradient(nu * grid, grid)
        dhdb[0] = nu[0]
        if not np.all(np.isfinite(nu)) or np.any(nu <= 0.0) or np.any(dhdb <= 0.0):
            raise RuntimeError("the inverted B-H law is not positive and monotone")
        self._grid, self._nu, self._dhdb = grid, nu, dhdb

    def reluctivity(self, magnitude: np.ndarray) -> np.ndarray:
        return np.interp(magnitude, self._grid, self._nu)

    def differential_reluctivity(self, magnitude: np.ndarray) -> np.ndarray:
        return np.interp(magnitude, self._grid, self._dhdb)


def constant_reluctivity(ne: int, iron: np.ndarray) -> np.ndarray:
    """Element reluctivity of the Newton Jacobian's constant part: nu0 off the iron, 0 on it."""
    values = np.full(int(ne), NU0)
    values[np.asarray(iron, dtype=np.int64)] = 0.0
    return values


class ElementCurl:
    """Element data of a lowest-order Nedelec space, from one native pass.

    ``LowestOrderCurlSystem`` evaluates, per straight tetrahedron, the six dofs
    (ascending), the element-constant curl ``C_e`` (3 x 6) of the space's own
    basis and the volume, and builds the element-graph matrix with each
    element's 36 value positions. With ``constant_reluctivity`` (one value per
    element) that matrix already holds ``sum_e c_e vol_e C_e^T C_e`` -- the
    Newton Jacobian's constant part, taken over by :class:`IronJacobian`.

    The element curl of ``A`` is ``C_e A[dofs_e]``; ``int_e f . curl v`` for an
    element-constant f is ``vol_e C_e^T f_e`` (``weight`` = ``volume``, the
    integral of the order-zero L2 function).
    """

    def __init__(self, fes, *, constant_reluctivity=None):
        system = _sparsesolv().LowestOrderCurlSystem(fes, coefficient=constant_reluctivity)
        self.ndof = int(fes.ndof)
        self.dofs = system["dofs"]
        self.coefficients = system["curl"]
        self.volume = system["volume"]
        self.weight = self.volume
        if np.any(self.volume <= 0.0):
            raise RuntimeError("non-positive element volumes")
        self.matrix = system["matrix"]
        self.positions = system["positions"]
        self.constant_reluctivity = (None if constant_reluctivity is None
                                     else np.asarray(constant_reluctivity, dtype=float).copy())

    def of(self, coefficients: np.ndarray) -> np.ndarray:
        """Element curl (ne, 3) of the coefficient vector."""
        return np.einsum("nkj,nj->nk", self.coefficients, coefficients[self.dofs])

    def pull_back(self, loads: np.ndarray) -> np.ndarray:
        """``sum_e w_e C_e^T loads_e`` (loads already carry vol / w factors)."""
        weights = np.einsum("nkj,nk->nj", self.coefficients, loads) * self.weight[:, None]
        return np.bincount(self.dofs.ravel(), weights=weights.ravel(), minlength=self.ndof)

    def element_coefficients(self, elements: np.ndarray):
        """Sorted edge dofs (n, 6) and curl coefficients (n, 3, 6) of the given tets."""
        return self.dofs[elements], self.coefficients[elements]


class IronJacobian:
    """Newton Jacobian: constant part once, iron element matrices in closed form.

    The matrix is the element-graph matrix of ``curl`` (built with
    :func:`constant_reluctivity`, so it holds nu0 curl-curl off the iron), plus
    the gauge mass when ``gauge_coefficient`` > 0. :meth:`refresh` adds
    ``vol_e C_e^T (nu_e I + q_e b_e b_e^T) C_e`` for every iron element at the
    element positions of the same pass.
    """

    def __init__(self, fes, curl: ElementCurl, iron: np.ndarray, gauge_coefficient: float):
        iron = np.asarray(iron, dtype=np.int64)
        expected = constant_reluctivity(len(curl.volume), iron)
        if curl.constant_reluctivity is None or not np.array_equal(curl.constant_reluctivity, expected):
            raise ValueError("IronJacobian needs ElementCurl(fes, constant_reluctivity="
                             "constant_reluctivity(mesh.ne, iron))")
        self.matrix = curl.matrix
        values = self.matrix.AsVector().FV().NumPy()
        if gauge_coefficient > 0.0:
            u, v = fes.TnT()
            mass = ng.BilinearForm(fes, symmetric=True)
            mass += gauge_coefficient * ng.InnerProduct(u, v) * ng.dx
            mass.Assemble()
            _, columns, offsets = self.matrix.CSR()
            mass_values, mass_columns, mass_offsets = mass.mat.CSR()
            if not (np.array_equal(np.asarray(columns), np.asarray(mass_columns))
                    and np.array_equal(np.asarray(offsets), np.asarray(mass_offsets))):
                raise RuntimeError("gauge mass pattern differs from the element graph")
            values += np.asarray(mass_values)
        # The native refresh saves the constant part of every row the iron
        # touches and rewrites those rows per Newton step, gathered row by row
        # (exactly symmetric, as the SPD direct factorisation requires).
        dofs, coefficients = curl.element_coefficients(iron)
        self._native = _sparsesolv().LowestOrderCurlJacobian(
            self.matrix, dofs, coefficients, curl.volume[iron], curl.positions[iron])

    def refresh(self, nu: np.ndarray, q: np.ndarray, b: np.ndarray) -> None:
        """Write the Jacobian for iron reluctivity ``nu``, rank-one ``q`` and flux ``b``."""
        self._native.Refresh(np.ascontiguousarray(nu, dtype=float), np.ascontiguousarray(q, dtype=float),
                             np.ascontiguousarray(b, dtype=float))


class _TrueResidualCG(CGSolver):
    """NGSolve CG recurrence stopped on the original free-DOF equation."""

    def __init__(self, *, rhs, solution, free, tolerance, **options):
        super().__init__(tol=tolerance, **options)
        self._rhs, self._solution, self._free = rhs, solution, free
        self._tolerance = tolerance
        self._rhs_norm = max(float(np.linalg.norm(rhs.FV().NumPy()[free])), 1e-300)
        self._true_residual = rhs.CreateVector()

    def CheckResidual(self, _preconditioned_residual):
        self.iterations += 1
        self._true_residual.data = self._rhs - self.mat * self._solution
        relative = float(np.linalg.norm(self._true_residual.FV().NumPy()[self._free]) / self._rhs_norm)
        if not math.isfinite(relative):
            raise RuntimeError("non-finite true CG residual")
        self.residuals.append(relative)
        return relative <= self._tolerance or self.iterations >= self.maxiter


def _material_numbers(mesh, names) -> np.ndarray:
    elements = mesh.ngmesh.Elements3D().NumPy()
    wanted = [index + 1 for index, name in enumerate(mesh.GetMaterials()) if name in set(names)]
    numbers = np.flatnonzero(np.isin(elements["index"], wanted)).astype(np.int64)
    if np.any(elements["np"][numbers] != 4):
        raise ValueError("the first-order Newton requires straight tetrahedra")
    return numbers


def solve_p1_newton(mesh, bh_table, *, iron=("iron",), source_cf=None, current_cf=None,
                    current_materials=None, dirichlet="outer", linear_solver="ams",
                    gauge_epsilon=None, beta_zero=None, reuse_hierarchy=True,
                    newton_tolerance=1e-3, field_tolerance=None, max_iterations=40,
                    max_halvings=6, cg_tolerance=1e-8, cg_max_iterations=2000, inexact=True,
                    linear_floor=True, ic_shift=1.05, source_projection_order=2,
                    observation_points=None, progress=None):
    """Nonlinear magnetostatics with lowest-order Nedelec A and Newton.

    Exactly one of ``source_cf`` (reduced A, vacuum source flux density B_s)
    and ``current_cf`` (total A, current density in A/m^2 on
    ``current_materials``) must be given. ``iron`` names the nonlinear
    materials; every other material has reluctivity nu0.

    ``linear_solver``: ``"ams"`` (default; with ``gauge_epsilon`` 0 it uses
    beta-zero AMS on the ungauged system, which needs a right-hand side
    orthogonal to the discrete gradients -- automatic for reduced A, checked
    and enforced for a current density), ``"iccg"`` (shifted IC(0) CG on the
    ungauged system) or ``"direct"`` (SPD sparse direct, gauged; a cross-check).
    ``gauge_epsilon`` (default 0 for ams/iccg, 1e-6 for direct) adds
    ``gauge_epsilon nu0 int A.v``.

    Convergence: relative residual <= ``newton_tolerance`` and, when
    ``field_tolerance`` is given, max change of |B| per iron element divided by
    the law's saturation flux <= ``field_tolerance``. Non-convergence raises.

    Must be called outside ``ngsolve.TaskManager`` (see module docstring).
    Returns ``A`` (GridFunction), ``B_cf``, ``stats`` (per-step history and
    phase timings) and, with ``observation_points``, ``observation_B_T``.
    """
    ssn = _sparsesolv()
    if ssn.TaskManagerActive():
        raise RuntimeError("solve_p1_newton must be called outside ngsolve.TaskManager: the native "
                           "AMS setup refuses to run inside one; the solver owns its parallel regions")
    if (source_cf is None) == (current_cf is None):
        raise ValueError("give exactly one of source_cf (reduced A) and current_cf (total A)")
    if linear_solver not in ("ams", "iccg", "direct"):
        raise ValueError("linear_solver must be 'ams', 'iccg' or 'direct'")
    if not 0.0 < float(newton_tolerance) < 1.0:
        raise ValueError("newton_tolerance must lie in (0, 1)")
    if field_tolerance is not None and not float(field_tolerance) > 0.0:
        raise ValueError("field_tolerance must be positive")
    if not (math.isfinite(cg_tolerance) and 0.0 < cg_tolerance < 1.0):
        raise ValueError("cg_tolerance must lie in (0, 1)")
    if gauge_epsilon is None:
        gauge_epsilon = 1.0e-6 if linear_solver == "direct" else 0.0
    gauge_epsilon = float(gauge_epsilon)
    if not math.isfinite(gauge_epsilon) or gauge_epsilon < 0.0:
        raise ValueError("gauge_epsilon must be finite and non-negative")
    if linear_solver == "direct" and gauge_epsilon == 0.0:
        raise ValueError("the direct solve needs gauge_epsilon > 0 (SPD system)")
    if beta_zero is None:
        beta_zero = linear_solver == "ams" and gauge_epsilon == 0.0
    if beta_zero and (linear_solver != "ams" or gauge_epsilon != 0.0):
        raise ValueError("beta_zero requires linear_solver='ams' and gauge_epsilon=0")
    if linear_solver == "ams" and gauge_epsilon == 0.0 and not beta_zero:
        raise ValueError("ungauged AMS needs beta_zero=True")
    if mesh.GetCurveOrder() != 1:
        raise ValueError("the first-order Newton requires a straight-sided mesh")
    iron_names = (iron,) if isinstance(iron, str) else tuple(iron)
    if not iron_names or not set(iron_names) <= set(mesh.GetMaterials()):
        raise ValueError(f"iron materials {iron_names} not all in the mesh")

    started = time.perf_counter()
    timing = {}
    law = FluxSideLaw(bh_table)
    iron_numbers = _material_numbers(mesh, iron_names)
    fes = ng.HCurl(mesh, order=1, dirichlet=dirichlet, nograds=True)
    free = np.fromiter(fes.FreeDofs(), dtype=bool, count=fes.ndof)
    with ng.TaskManager():
        t0 = time.perf_counter()
        curl = ElementCurl(fes, constant_reluctivity=constant_reluctivity(mesh.ne, iron_numbers))
        timing["element_system_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        gradient = lowest_order_gradient(fes)
        timing["gradient_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        load = np.zeros(fes.ndof)
        source_mean = None
        if source_cf is not None:
            projection = ng.GridFunction(ng.L2(mesh, order=int(source_projection_order), dim=3,
                                               definedon=mesh.Materials("|".join(map(re.escape, iron_names)))))
            projection.Set(source_cf, definedon=mesh.Materials("|".join(map(re.escape, iron_names))))
            means = np.stack([np.asarray(ng.Integrate(projection[k], mesh, ng.VOL, element_wise=True))
                              for k in range(3)], axis=1)
            source_mean = means[iron_numbers] / curl.volume[iron_numbers, None]
        else:
            names = current_materials
            if names is None:
                raise ValueError("current_materials names where current_cf lives")
            names = (names,) if isinstance(names, str) else tuple(names)
            if not set(names) <= set(mesh.GetMaterials()):
                raise ValueError(f"current materials {names} not all in the mesh")
            v = fes.TestFunction()
            current_form = ng.LinearForm(fes)
            current_form += ng.InnerProduct(current_cf, v) * ng.dx(
                definedon=mesh.Materials("|".join(map(re.escape, names))))
            current_form.Assemble()
            load = current_form.vec.FV().NumPy().copy()
            if gauge_epsilon == 0.0:
                # The ungauged curl-curl system is solvable only for a load
                # orthogonal to the discrete gradients (a weakly divergence-free J).
                from scipy.sparse import csr_matrix

                values, columns, offsets = gradient.CSR()
                g = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                               shape=(gradient.height, gradient.width))
                free_load = np.where(free, load, 0.0)
                defect = np.linalg.norm(g.T @ free_load)
                scale = np.linalg.norm(abs(g).T @ np.abs(free_load))
                if scale > 0.0 and defect > 1e-10 * scale:
                    raise ValueError(f"current load is not orthogonal to the discrete gradients "
                                     f"(relative defect {defect / scale:.2e}); use a weakly "
                                     "divergence-free current (radia.meshed_current) or gauge_epsilon > 0")
        timing["source_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        jacobian = IronJacobian(fes, curl, iron_numbers, gauge_epsilon * NU0)
        timing["jacobian_setup_s"] = time.perf_counter() - t0
        mass = None
        if gauge_epsilon > 0.0:
            from scipy.sparse import csr_matrix

            a, b = fes.TnT()
            mass_form = ng.BilinearForm(fes, symmetric=False)
            mass_form += gauge_epsilon * NU0 * ng.InnerProduct(a, b) * ng.dx
            mass_form.Assemble()
            values, columns, offsets = mass_form.mat.CSR()
            mass = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                              shape=(fes.ndof, fes.ndof))
    coordinates = np.asarray(mesh.ngmesh.Coordinates(), dtype=float)
    nu_all = np.full(mesh.ne, NU0)
    # curl.of already divides by w, so both terms pull back with vol / w.
    scale_all = curl.volume / curl.weight
    source_scale = scale_all[iron_numbers]

    def element_flux(coefficients):
        b = curl.of(coefficients)[iron_numbers]
        if source_mean is not None:
            b = b + source_mean
        return b, np.linalg.norm(b, axis=1)

    def material(b, magnitude):
        nu = law.reluctivity(magnitude)
        dhdb = law.differential_reluctivity(magnitude)
        safe = np.where(magnitude > 1.0e-9, magnitude, 1.0)
        q = np.where(magnitude > 1.0e-9, (dhdb - nu) / (safe * safe), 0.0)
        return nu, q

    def residual(coefficients, nu_iron):
        nu = nu_all.copy()
        nu[iron_numbers] = nu_iron
        curl_a = curl.of(coefficients)
        loads = (nu * scale_all)[:, None] * curl_a
        if source_mean is not None:
            loads[iron_numbers] += ((nu_iron - NU0) * source_scale)[:, None] * source_mean
        values = curl.pull_back(loads) - load
        if mass is not None:
            values += mass @ coefficients
        return values, float(np.linalg.norm(values[free]))

    solution = ng.GridFunction(fes, name="A_p1_newton")
    solution.vec[:] = 0.0
    x = solution.vec.FV().NumPy()
    b, magnitude = element_flux(x)
    nu, q = material(b, magnitude)
    r, r_norm = residual(x, nu)
    r0 = r_norm
    rhs = solution.vec.CreateVector()
    update = solution.vec.CreateVector()
    ams = None
    history, converged, change = [], r0 == 0.0, 0.0
    for iteration in range(1, int(max_iterations) + 1):
        if converged:
            break
        row = {"iteration": iteration, "residual_relative": r_norm / r0}
        t0 = time.perf_counter()
        with ng.TaskManager():
            jacobian.refresh(nu, q, b)
        row["assembly_s"] = time.perf_counter() - t0
        tolerance = float(cg_tolerance)
        if inexact:
            tolerance = max(tolerance, min(0.01, 0.1 * r_norm / r0))
        if linear_floor:
            tolerance = min(0.5, max(tolerance, 0.1 * float(newton_tolerance) * r0 / r_norm))
        row["linear_tolerance"] = tolerance
        rhs.FV().NumPy()[:] = np.where(free, -r, 0.0)
        update[:] = 0.0
        t0 = time.perf_counter()
        if linear_solver == "ams":
            if ams is None:
                ams = ssn.HypreBasedAMSPreconditioner(
                    mat=jacobian.matrix, grad_mat=gradient, freedofs=fes.FreeDofs(),
                    coord_x=coordinates[:, 0].tolist(), coord_y=coordinates[:, 1].tolist(),
                    coord_z=coordinates[:, 2].tolist(), cycle_type=1, print_level=0,
                    beta_zero=bool(beta_zero), reuse_hierarchy=bool(reuse_hierarchy))
            else:
                ams.Update(jacobian.matrix)
            row["preconditioner_s"] = time.perf_counter() - t0
            t0 = time.perf_counter()
            cg = _TrueResidualCG(mat=jacobian.matrix, pre=ams, maxiter=int(cg_max_iterations),
                                 rhs=rhs, solution=update, free=free, tolerance=tolerance,
                                 printrates=False)
            with ng.TaskManager():
                cg.Solve(rhs=rhs, sol=update, initialize=True)
            row["cg_iterations"] = int(cg.iterations)
        elif linear_solver == "iccg":
            solver = ssn.SparseSolvSolver(jacobian.matrix, method="ICCG", freedofs=fes.FreeDofs(),
                                          tol=tolerance, maxiter=int(cg_max_iterations),
                                          shift=float(ic_shift), save_best_result=True, use_abmc=True)
            row["preconditioner_s"] = time.perf_counter() - t0
            t0 = time.perf_counter()
            # The IC solver stops on its own residual measure; continue (warm)
            # with a tighter internal target until the true residual meets
            # the Newton inner tolerance.
            iterations, probe = 0, rhs.CreateVector()
            rhs_norm = max(float(np.linalg.norm(rhs.FV().NumPy()[free])), 1e-300)
            for _ in range(6):
                with ng.TaskManager():
                    result = solver.Solve(rhs, update)
                    probe.data = rhs - jacobian.matrix * update
                iterations += int(result.iterations)
                if float(np.linalg.norm(probe.FV().NumPy()[free])) / rhs_norm <= tolerance:
                    break
                solver.tol = solver.tol * 0.1
            row["cg_iterations"] = iterations
        else:
            copy = jacobian.matrix.CreateMatrix()
            copy.AsVector().data = jacobian.matrix.AsVector()
            row["preconditioner_s"] = 0.0
            t0 = time.perf_counter()
            with ng.TaskManager():
                from radia.vector_potential_solver import direct_inverse_type

                update.data = copy.Inverse(fes.FreeDofs(), inverse=direct_inverse_type()) * rhs
            row["cg_iterations"] = None
        row["solve_s"] = time.perf_counter() - t0
        check = rhs.CreateVector()
        with ng.TaskManager():
            check.data = rhs - jacobian.matrix * update
        rhs_norm = max(float(np.linalg.norm(rhs.FV().NumPy()[free])), 1e-300)
        linear = float(np.linalg.norm(check.FV().NumPy()[free]) / rhs_norm)
        row["linear_relative_residual"] = linear
        if not math.isfinite(linear) or linear > tolerance * (1.0 + 1e-9):
            raise RuntimeError(f"Newton linear solve ({linear_solver}) reached {linear:.2e}, "
                               f"above its tolerance {tolerance:.2e}")
        step = update.FV().NumPy().copy()
        t0 = time.perf_counter()
        alpha, accepted = 1.0, None
        for halving in range(int(max_halvings) + 1):
            trial = x + alpha * step
            b_trial, magnitude_trial = element_flux(trial)
            nu_trial, q_trial = material(b_trial, magnitude_trial)
            r_trial, n_trial = residual(trial, nu_trial)
            if n_trial <= (1.0 - 1.0e-4 * alpha) * r_norm:
                accepted = (trial, b_trial, magnitude_trial, nu_trial, q_trial, r_trial, n_trial)
                break
            alpha *= 0.5
        if accepted is None:
            raise RuntimeError(f"Newton line search failed at iteration {iteration} "
                               f"(relative residual {r_norm / r0:.2e})")
        trial, b_new, magnitude_new, nu, q, r, r_norm = accepted
        row["line_search_s"] = time.perf_counter() - t0
        row["step_length"] = alpha
        change = float(np.max(np.abs(magnitude_new - magnitude)) / law.b_saturation)
        row["relative_B_change"] = change
        row["residual_relative_after"] = r_norm / r0
        x[:] = trial
        b, magnitude = b_new, magnitude_new
        history.append(row)
        if progress is not None:
            progress(dict(row))
        converged = (r_norm <= float(newton_tolerance) * r0
                     and (field_tolerance is None or change <= float(field_tolerance)))
    stats = {"method": "first-order Newton (closed-form iron Jacobian)", "converged": bool(converged),
             "iterations": len(history), "final_residual_relative": (r_norm / r0) if r0 else 0.0,
             "final_relative_B_change": change, "newton_tolerance": float(newton_tolerance),
             "field_tolerance": field_tolerance, "linear_solver": linear_solver,
             "gauge_epsilon": gauge_epsilon, "beta_zero": bool(beta_zero),
             "reuse_hierarchy": bool(reuse_hierarchy), "inexact": bool(inexact),
             "linear_floor": bool(linear_floor), "ndof": int(fes.ndof),
             "iron_elements": int(len(iron_numbers)), "history": history, "setup_s": timing,
             "total_s": time.perf_counter() - started}
    if not converged:
        raise RuntimeError(f"first-order Newton did not converge in {max_iterations} iterations "
                           f"(relative residual {r_norm / r0:.2e}, B change {change:.2e})")
    B_cf = ng.curl(solution) if source_cf is None else source_cf + ng.curl(solution)
    nu_gf = ng.GridFunction(ng.L2(mesh, order=0), name="reluctivity")
    nu_values = np.full(mesh.ne, NU0)
    nu_values[iron_numbers] = nu
    nu_gf.vec.FV().NumPy()[:] = nu_values
    result = {"A": solution, "B_cf": B_cf, "nu": nu_gf, "H_cf": nu_gf * B_cf,
              "stats": stats, "fes": fes, "law": law}
    if observation_points is not None:
        points = np.asarray(observation_points, dtype=float).reshape(-1, 3)
        result["observation_B_T"] = np.asarray([[float(c) for c in B_cf(mesh(*p))] for p in points])
    return result
