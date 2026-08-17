"""Planar HDiv-VIM reduced motor model.

The rotor is represented in its own local frame.  Its symmetric charge Gram
is built once and reused while the applied stator MMF rotates through that
frame.  This is the production saliency/reluctance-motor path; the full
transient A-formulation remains a separate motor analysis.

Callers own the surrounding ``ngsolve.TaskManager`` region.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass

import ngsolve as ng
import numpy as np

from .force import force_torque_result
from .vim import PlanarDemagBody, maxwell_torque_circle

MU0 = 4.0e-7 * np.pi


def _rotation(angle: float) -> np.ndarray:
    c = math.cos(float(angle))
    s = math.sin(float(angle))
    return np.array(((c, -s), (s, c)), dtype=float)


def _vector2(value, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size != 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain two finite components")
    return result


@dataclass(slots=True)
class HDivMotorState:
    """Solved magnetization at one mechanical rotor angle."""

    rotor_angle_rad: float
    field_global_A_per_m: np.ndarray
    field_local_A_per_m: np.ndarray
    coefficients: np.ndarray
    magnetization_average_local_A_per_m: np.ndarray
    magnetization_average_global_A_per_m: np.ndarray
    coenergy_J: float
    torque_volume_Nm: float
    linear_iterations: int


class HDivReducedMotor:
    """Single-rotor planar HDiv-VIM motor with one cached demag operator.

    ``mesh`` is the rotor iron cross-section in its local reference frame.
    ``mu_r`` is linear in this first production path.  The applied stator MMF
    is supplied as a global two-component H field for each solve.
    """

    def __init__(self, mesh, mu_r: float, *, stack_length: float = 1.0,
                 center=(0.0, 0.0), eta: float = 2.0,
                 cg_tol: float = 1e-10, cg_maxit: int = 5000,
                 order: int = 1):
        if mesh.dim != 2:
            raise ValueError(f"HDivReducedMotor: mesh.dim must be 2 (got {mesh.dim})")
        if not float(mu_r) > 1.0:
            raise ValueError(f"HDivReducedMotor: mu_r must be > 1 (got {mu_r!r})")
        if not float(stack_length) > 0.0:
            raise ValueError("HDivReducedMotor: stack_length must be positive")
        self.mesh = mesh
        self.mu_r = float(mu_r)
        self.chi = self.mu_r - 1.0
        self.stack_length = float(stack_length)
        self.center = _vector2(center, "center")
        self.body = PlanarDemagBody(
            mesh, order=order, eta=eta, cg_tol=cg_tol, cg_maxit=cg_maxit)
        self.order = self.body.order
        self.area = float(np.sum(self.body.areas))
        self.gram_build_count = 1

    def solve_angle(self, rotor_angle: float, field_global) -> HDivMotorState:
        """Solve one angle while reusing the body's charge Gram."""
        angle = float(rotor_angle)
        if not math.isfinite(angle):
            raise ValueError("rotor_angle must be finite")
        field_global = _vector2(field_global, "field_global")
        rotation = _rotation(angle)
        field_local = rotation.T @ field_global
        projected = self.body.project(
            ng.CoefficientFunction((float(field_local[0]), float(field_local[1]))))
        coefficients = self.body.solve_linear(self.chi, projected)
        magnetization_local = np.asarray(self.body.M_avg(coefficients), dtype=float)
        magnetization_global = rotation @ magnetization_local
        coenergy = 0.5 * MU0 * float(
            coefficients @ (self.body.Mm @ projected)) * self.stack_length
        torque_volume = MU0 * self.area * float(
            magnetization_local[0] * field_local[1]
            - magnetization_local[1] * field_local[0]) * self.stack_length
        return HDivMotorState(
            rotor_angle_rad=angle,
            field_global_A_per_m=field_global.copy(),
            field_local_A_per_m=field_local,
            coefficients=coefficients,
            magnetization_average_local_A_per_m=magnetization_local,
            magnetization_average_global_A_per_m=magnetization_global,
            coenergy_J=coenergy,
            torque_volume_Nm=torque_volume,
            linear_iterations=self.body.last_linear_iterations,
        )

    def source_field_global(self, points, state: HDivMotorState) -> np.ndarray:
        """Evaluate the solved rotor field at global observation points."""
        points = np.ascontiguousarray(points, dtype=float).reshape(-1, 2)
        rotation = _rotation(state.rotor_angle_rad)
        local_points = (points-self.center) @ rotation + self.center
        local_field = self.body.H_at(local_points, state.coefficients)
        return local_field @ rotation.T

    def source_field_cf(self, state: HDivMotorState, *, target_angle: float = 0.0):
        """Native source H in another body's local frame for mutual coupling."""
        return self.body.field_cf(
            state.coefficients,
            source_angle=state.rotor_angle_rad,
            target_angle=float(target_angle),
            center=tuple(self.center),
        )

    def maxwell_torque(self, state: HDivMotorState, radius: float,
                       *, circle_points: int = 1440) -> float:
        """Air-gap Maxwell-stress torque for the solved state."""
        if not float(radius) > 0.0:
            raise ValueError("maxwell torque radius must be positive")
        if int(circle_points) < 8:
            raise ValueError("circle_points must be >= 8")

        def total_field(points):
            return (self.source_field_global(points, state)
                    + state.field_global_A_per_m)

        return self.stack_length * maxwell_torque_circle(
            total_field, float(radius), n=int(circle_points),
            center=tuple(self.center))

    def virtual_work_torque(self, rotor_angle: float, field_global,
                            *, delta_angle: float = math.radians(0.25)) -> float:
        """Central derivative of fixed-current magnetic coenergy."""
        delta = float(delta_angle)
        if not delta > 0.0:
            raise ValueError("virtual-work delta_angle must be positive")
        minus = self.solve_angle(float(rotor_angle)-delta, field_global)
        plus = self.solve_angle(float(rotor_angle)+delta, field_global)
        return (plus.coenergy_J-minus.coenergy_J)/(2.0*delta)

    def sweep(self, rotor_angles, field_global, *, maxwell_radius: float,
              circle_points: int = 1440,
              energy_delta_angle: float = math.radians(0.25)) -> dict:
        """Run a mechanical-angle sweep with three independent torque reads."""
        angles = np.asarray(rotor_angles, dtype=float).reshape(-1)
        if angles.size == 0 or not np.all(np.isfinite(angles)):
            raise ValueError("rotor_angles must contain finite values")
        field_global = _vector2(field_global, "field_global")
        started = time.perf_counter()
        rows = []
        for angle in angles:
            state = self.solve_angle(float(angle), field_global)
            torque_maxwell = self.maxwell_torque(
                state, maxwell_radius, circle_points=circle_points)
            torque_energy = self.virtual_work_torque(
                float(angle), field_global, delta_angle=energy_delta_angle)
            scale = max(abs(torque_maxwell), abs(torque_energy),
                        abs(state.torque_volume_Nm), 1e-30)
            pivot = [float(self.center[0]), float(self.center[1]), 0.0]
            force_torque_results = {
                "maxwell_surface": force_torque_result(
                    None,
                    [0.0, 0.0, torque_maxwell],
                    method="maxwell_surface_stress_air",
                    frame="global_cartesian",
                    pivot_m=pivot,
                    dimensionality="2d_planar",
                ),
                "magnetization_volume": force_torque_result(
                    None,
                    [0.0, 0.0, state.torque_volume_Nm],
                    method="magnetization_volume_moment",
                    frame="global_cartesian",
                    pivot_m=pivot,
                    dimensionality="2d_planar",
                ),
                "virtual_work": force_torque_result(
                    None,
                    [0.0, 0.0, torque_energy],
                    method="coenergy_virtual_work",
                    frame="global_cartesian",
                    pivot_m=pivot,
                    dimensionality="2d_planar",
                ),
            }
            rows.append({
                "rotor_angle_rad": float(angle),
                "rotor_angle_deg": math.degrees(float(angle)),
                "H_local_A_per_m": state.field_local_A_per_m.tolist(),
                "M_average_local_A_per_m": (
                    state.magnetization_average_local_A_per_m.tolist()),
                "coenergy_J": state.coenergy_J,
                "torque_maxwell_Nm": torque_maxwell,
                "torque_volume_Nm": state.torque_volume_Nm,
                "torque_virtual_work_Nm": torque_energy,
                "force_torque_results": force_torque_results,
                "torque_spread_relative": (
                    max(torque_maxwell, torque_energy, state.torque_volume_Nm)
                    - min(torque_maxwell, torque_energy, state.torque_volume_Nm)) / scale,
                "linear_iterations": state.linear_iterations,
            })
        return {
            "analysis": "hdiv_reduced_motor",
            "formulation": f"planar BDM{self.order} HDiv-VIM reluctance motor",
            "hdiv_order": self.order,
            "geometry_order": self.body.geometry_order,
            "mu_r": self.mu_r,
            "stack_length_m": self.stack_length,
            "center_m": self.center.tolist(),
            "field_global_A_per_m": field_global.tolist(),
            "ndof": self.body.ndof,
            "n_charge": self.body.n_charge,
            "n_elements": self.body.nel,
            "gram_build_count": self.gram_build_count,
            "gram_stats": dict(self.body.G.stats()),
            "maxwell_radius_m": float(maxwell_radius),
            "circle_points": int(circle_points),
            "energy_delta_angle_rad": float(energy_delta_angle),
            "force_result_schema": "radia.force-result/v1",
            "torque_axis": [0.0, 0.0, 1.0],
            "independent_torque_methods": [
                "maxwell_surface_stress_air",
                "magnetization_volume_moment",
                "coenergy_virtual_work",
            ],
            "angles": rows,
            "elapsed_s": time.perf_counter()-started,
        }


__all__ = ["HDivMotorState", "HDivReducedMotor"]
