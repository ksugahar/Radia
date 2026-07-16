"""Independent-space block coupling for HDiv permanent magnets and iron."""

from dataclasses import dataclass, field

import ngsolve as ng
import numpy as np

from ._field_batch import (
    field_coefficient_from_solution,
    field_from_solution,
)
from ._hysteresis import SolveHysteresis, _normalize_h_steps
from ._solve import hdiv_demag_solve


@dataclass
class CoupledBody:
    """One independently meshed body in an HDiv block solve.

    Give either ``bh_table`` for nonlinear iron, or ``mu_r`` for a linear body.
    Adding ``B_r`` selects the linear-recoil permanent-magnet law.  Independent
    spaces are mandatory so physical normal-magnetization jumps survive at
    touching iron/PM and segmented-PM interfaces.
    """

    mesh: object
    name: str
    mu_r: object = None
    B_r: object = None
    bh_table: object = None
    order: int = 1
    applied_field: object = None
    solve_options: dict = field(default_factory=dict)

    def __post_init__(self):
        if int(getattr(self.mesh, "dim", -1)) != 3:
            raise ValueError("vim.CoupledBody requires a 3D NGSolve mesh")
        if not self.name or not isinstance(self.name, str):
            raise ValueError("vim.CoupledBody requires a non-empty name")
        if (self.mu_r is None) == (self.bh_table is None):
            raise ValueError(
                "vim.CoupledBody requires exactly one of mu_r or bh_table")
        if self.B_r is not None and self.bh_table is not None:
            raise ValueError("vim.CoupledBody B_r cannot be combined with bh_table")
        self.order = int(self.order)
        if self.order not in (1, 2):
            raise ValueError("vim.CoupledBody order must be 1 or 2")
        self.solve_options = dict(self.solve_options)
        forbidden = {"mu_r", "B_r", "bh_table", "H_ext", "order", "_prepared_operator"}
        overlap = sorted(forbidden.intersection(self.solve_options))
        if overlap:
            raise ValueError(
                "vim.CoupledBody solve_options must not override %s" % overlap)

    @property
    def kind(self):
        if self.B_r is not None:
            return "linear-recoil-pm"
        if self.bh_table is not None:
            return "nonlinear-iron"
        return "linear-iron"


@dataclass
class CoupledHistoryBody:
    """One stateful B-input permanent magnet in a coupled history solve."""

    mesh: object
    name: str
    material: object
    order: int = 1
    applied_field: object = None
    initial_b_path: object = None
    initial_state: object = None
    solve_options: dict = field(default_factory=dict)

    def __post_init__(self):
        if int(getattr(self.mesh, "dim", -1)) != 3:
            raise ValueError("vim.CoupledHistoryBody requires a 3D NGSolve mesh")
        if not self.name or not isinstance(self.name, str):
            raise ValueError("vim.CoupledHistoryBody requires a non-empty name")
        required = ("state0", "forward", "commit", "nu_bound")
        missing = [name for name in required if not callable(
            getattr(self.material, name, None))]
        if missing:
            raise ValueError(
                "vim.CoupledHistoryBody material is missing %s" % missing)
        if self.initial_b_path is not None and self.initial_state is not None:
            raise ValueError(
                "vim.CoupledHistoryBody initial_b_path and initial_state are mutually exclusive")
        self.order = int(self.order)
        if self.order not in (1, 2):
            raise ValueError("vim.CoupledHistoryBody order must be 1 or 2")
        self.solve_options = dict(self.solve_options)
        forbidden = {
            "play", "material", "initial_b_path", "initial_state", "order",
            "_prepared_operator",
        }
        overlap = sorted(forbidden.intersection(self.solve_options))
        if overlap:
            raise ValueError(
                "vim.CoupledHistoryBody solve_options must not override %s" % overlap)


def _field3(value, label):
    if value is None:
        return ng.CF((0.0, 0.0, 0.0))
    if getattr(value, "dim", None) == 3:
        return value
    array = np.asarray(value, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError("%s must be a finite 3-vector or NGSolve CoefficientFunction" % label)
    return ng.CF(tuple(float(component) for component in array))


def solve_coupled(bodies, H_ext=None, *, tol=1.0e-6, maxit=50):
    """Solve mutually interacting independently meshed HDiv bodies.

    The block iteration is Gauss-Seidel.  Each geometry-only ChargeGram is
    prepared once and reused; only the body material solve and cross-body C++
    field evaluations repeat.  Non-convergence raises instead of returning a
    partially coupled state.  The caller owns ``with ngsolve.TaskManager():``.
    """
    bodies = tuple(bodies)
    if len(bodies) < 2 or not all(isinstance(body, CoupledBody) for body in bodies):
        raise ValueError("vim.SolveCoupled requires at least two CoupledBody objects")
    names = [body.name for body in bodies]
    if len(set(names)) != len(names):
        raise ValueError("vim.SolveCoupled body names must be unique")
    meshes = [body.mesh for body in bodies]
    if len({id(mesh) for mesh in meshes}) != len(meshes):
        raise ValueError(
            "vim.SolveCoupled requires a distinct mesh object/HDiv space per body")
    tol = float(tol)
    maxit = int(maxit)
    if not (tol > 0.0 and maxit > 0):
        raise ValueError("vim.SolveCoupled requires tol > 0 and maxit > 0")
    global_field = _field3(H_ext, "vim.SolveCoupled H_ext")

    results = [None] * len(bodies)
    relative_step = float("inf")
    history = []
    for iteration in range(1, maxit + 1):
        previous = [
            None if result is None else np.asarray(
                result["_m_coefficients"], dtype=float).copy()
            for result in results
        ]
        for index, body in enumerate(bodies):
            applied = global_field + _field3(
                body.applied_field, "vim.CoupledBody.applied_field")
            for other_index, other_result in enumerate(results):
                if other_index != index and other_result is not None:
                    applied = applied + field_coefficient_from_solution(
                        other_result, algorithm="direct")
            kwargs = dict(body.solve_options)
            kwargs.update(
                mu_r=body.mu_r, B_r=body.B_r, bh_table=body.bh_table,
                H_ext=applied, order=body.order,
                _prepared_operator=(None if results[index] is None
                                    else results[index]["_prepared_operator"]))
            results[index] = hdiv_demag_solve(body.mesh, **kwargs)

        if all(coefficients is not None for coefficients in previous):
            delta_squared = 0.0
            scale_squared = 0.0
            body_steps = {}
            for body, old, result in zip(bodies, previous, results):
                new = np.asarray(result["_m_coefficients"], dtype=float)
                delta = float(np.linalg.norm(new-old))
                scale = float(np.linalg.norm(new))
                body_steps[body.name] = delta/max(scale, 1.0e-300)
                delta_squared += delta*delta
                scale_squared += scale*scale
            relative_step = np.sqrt(delta_squared/max(scale_squared, 1.0e-300))
            history.append(dict(
                iteration=iteration, relative_step=float(relative_step),
                body_relative_steps=body_steps))
            if relative_step < tol:
                break
    else:
        raise RuntimeError(
            "vim.SolveCoupled did not converge in %d block iterations "
            "(relative step %.3e > %.3e)" % (maxit, relative_step, tol))

    return dict(
        bodies=tuple(results), body_names=tuple(names), iterations=int(iteration),
        relative_step=float(relative_step), converged=True,
        block_solver="gauss-seidel", convergence_history=history,
        permanent_magnet_body_count=sum(body.B_r is not None for body in bodies),
        nonlinear_iron_body_count=sum(body.bh_table is not None for body in bodies),
        body_kinds={body.name: body.kind for body in bodies},
        _body_specs=bodies,
    )


def _relative_coefficient_step(previous, current, names):
    if any(value is None for value in previous):
        return float("inf"), {}
    delta_squared = 0.0
    scale_squared = 0.0
    body_steps = {}
    for name, old, new in zip(names, previous, current):
        delta = float(np.linalg.norm(new-old))
        scale = float(np.linalg.norm(new))
        body_steps[name] = delta/max(scale, 1.0e-300)
        delta_squared += delta*delta
        scale_squared += scale*scale
    relative_step = np.sqrt(delta_squared/max(scale_squared, 1.0e-300))
    return float(relative_step), body_steps


def solve_coupled_hysteresis(history_body, bodies, h_steps, H_ext=None, *,
                             tol=1.0e-6, maxit=50):
    """Advance one stateful PM coupled to linear/nonlinear HDiv bodies.

    Every outer trial starts from the same committed history state.  The trial
    state is accepted only after the all-body coefficient fixed point converges,
    preventing one physical field step from committing the constitutive model
    multiple times.  The caller owns ``with ngsolve.TaskManager():``.
    """
    if not isinstance(history_body, CoupledHistoryBody):
        raise TypeError(
            "vim.SolveCoupledHysteresis requires a CoupledHistoryBody first")
    bodies = tuple(bodies)
    if not bodies or not all(isinstance(body, CoupledBody) for body in bodies):
        raise ValueError(
            "vim.SolveCoupledHysteresis requires at least one CoupledBody")
    names = [history_body.name] + [body.name for body in bodies]
    if len(set(names)) != len(names):
        raise ValueError("vim.SolveCoupledHysteresis body names must be unique")
    meshes = [history_body.mesh] + [body.mesh for body in bodies]
    if len({id(mesh) for mesh in meshes}) != len(meshes):
        raise ValueError(
            "vim.SolveCoupledHysteresis requires a distinct mesh object/HDiv space per body")
    tol = float(tol)
    maxit = int(maxit)
    if not (tol > 0.0 and maxit > 0):
        raise ValueError(
            "vim.SolveCoupledHysteresis requires tol > 0 and maxit > 0")

    global_field = _field3(H_ext, "vim.SolveCoupledHysteresis H_ext")
    step_records = _normalize_h_steps(h_steps)
    committed_state = history_body.initial_state
    history_result = None
    body_results = [None] * len(bodies)
    history_operator = None
    outputs = []

    for step_index, (step_field, step_uniform, _) in enumerate(step_records):
        applied_global = global_field + step_field
        convergence_history = []
        relative_step = float("inf")
        for iteration in range(1, maxit + 1):
            previous_results = [history_result] + list(body_results)
            previous = [
                None if result is None else np.asarray(
                    result["_m_coefficients"], dtype=float).copy()
                for result in previous_results
            ]

            pm_applied = applied_global + _field3(
                history_body.applied_field,
                "vim.CoupledHistoryBody.applied_field")
            for result in body_results:
                if result is not None:
                    pm_applied = pm_applied + field_coefficient_from_solution(
                        result, algorithm="direct")
            history_kwargs = dict(history_body.solve_options)
            history_kwargs.update(
                material=history_body.material, order=history_body.order,
                initial_state=committed_state,
                initial_b_path=(
                    history_body.initial_b_path
                    if committed_state is None else None),
                _prepared_operator=history_operator,
            )
            history_result = SolveHysteresis(
                history_body.mesh, [pm_applied], **history_kwargs)
            history_operator = history_result["_prepared_operator"]

            for index, body in enumerate(bodies):
                applied = applied_global + _field3(
                    body.applied_field, "vim.CoupledBody.applied_field")
                applied = applied + field_coefficient_from_solution(
                    history_result, algorithm="direct")
                for other_index, other_result in enumerate(body_results):
                    if other_index != index and other_result is not None:
                        applied = applied + field_coefficient_from_solution(
                            other_result, algorithm="direct")
                kwargs = dict(body.solve_options)
                kwargs.update(
                    mu_r=body.mu_r, B_r=body.B_r, bh_table=body.bh_table,
                    H_ext=applied, order=body.order,
                    _prepared_operator=(
                        None if body_results[index] is None
                        else body_results[index]["_prepared_operator"]),
                )
                body_results[index] = hdiv_demag_solve(body.mesh, **kwargs)

            current_results = [history_result] + list(body_results)
            current = [np.asarray(
                result["_m_coefficients"], dtype=float) for result in current_results]
            relative_step, body_steps = _relative_coefficient_step(
                previous, current, names)
            convergence_history.append(dict(
                iteration=iteration, relative_step=relative_step,
                body_relative_steps=body_steps))
            if relative_step < tol:
                break
        else:
            raise RuntimeError(
                "vim.SolveCoupledHysteresis step %d did not converge in %d block "
                "iterations (relative step %.3e > %.3e)"
                % (step_index, maxit, relative_step, tol))

        committed_state = history_result["state"]
        outputs.append(dict(
            step_index=int(step_index),
            h_applied=(None if step_uniform is None else step_uniform.copy()),
            history_body=history_result, bodies=tuple(body_results),
            iterations=int(iteration), relative_step=float(relative_step),
            convergence_history=convergence_history,
        ))

    return dict(
        steps=outputs, history_body=history_result, bodies=tuple(body_results),
        body_names=tuple(names), state=committed_state,
        converged=True, block_solver="gauss-seidel-history-trial",
        history_material_model=getattr(
            history_body.material, "permanent_magnet_model", "custom-b-input"),
        history_material_level=getattr(
            history_body.material, "permanent_magnet_level", None),
        nonlinear_iron_body_count=sum(body.bh_table is not None for body in bodies),
        _history_spec=history_body, _body_specs=bodies,
    )


def field_from_coupled_solution(result, points, algorithm="auto"):
    """Sum the magnetization fields of all bodies at ``points``."""
    if not isinstance(result, dict) or "bodies" not in result:
        raise TypeError("vim.FieldFromCoupledSolution requires SolveCoupled's result")
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    total = np.zeros((len(points), 3), dtype=float)
    for body_result in result["bodies"]:
        total += field_from_solution(body_result, points, algorithm=algorithm)
    return total


def field_from_coupled_hysteresis(result, points, algorithm="auto"):
    """Sum the final history-body and ordinary-body magnetization fields."""
    if not isinstance(result, dict) or "history_body" not in result:
        raise TypeError(
            "vim.FieldFromCoupledHysteresis requires SolveCoupledHysteresis's result")
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    total = field_from_solution(
        result["history_body"], points, algorithm=algorithm)
    for body_result in result["bodies"]:
        total += field_from_solution(body_result, points, algorithm=algorithm)
    return total
