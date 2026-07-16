"""Independent-space block coupling for HDiv permanent magnets and iron."""

from dataclasses import dataclass, field

import ngsolve as ng
import numpy as np

from ._field_batch import (
    field_coefficient_from_solution,
    field_from_solution,
)
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


def field_from_coupled_solution(result, points, algorithm="auto"):
    """Sum the magnetization fields of all bodies at ``points``."""
    if not isinstance(result, dict) or "bodies" not in result:
        raise TypeError("vim.FieldFromCoupledSolution requires SolveCoupled's result")
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    total = np.zeros((len(points), 3), dtype=float)
    for body_result in result["bodies"]:
        total += field_from_solution(body_result, points, algorithm=algorithm)
    return total

