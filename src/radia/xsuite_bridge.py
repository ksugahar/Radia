"""Thin Radia field-map bridge to CERN Xsuite particle tracking.

Radia owns the electromagnetic field solve.  Xsuite owns particle state and
the spatial Boris integrator.  Keeping that boundary explicit avoids growing a
second particle-tracking engine inside Radia while still allowing a solved
Radia object to be used directly by accelerator workflows.

All coordinates are metres and magnetic fields are tesla, matching the Radia
and Xsuite APIs used here.  The bridge currently covers magnetic tracking.  It
does not claim electrostatic acceleration or particle-matter interactions;
those require a separate Xsuite/Xfields or Xcoll lane.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


FieldEvaluator = Callable[[int, str, object], object]
FieldMap = Callable[[object, object, object], tuple[object, object, object]]


@dataclass(frozen=True)
class AxisAlignedBox:
    """Axis-aligned tracking boundary in metres."""

    minimum: tuple[float, float, float]
    maximum: tuple[float, float, float]

    def __post_init__(self) -> None:
        if len(self.minimum) != 3 or len(self.maximum) != 3:
            raise ValueError("minimum and maximum must contain three coordinates")
        if not all(np.isfinite((*self.minimum, *self.maximum))):
            raise ValueError("boundary coordinates must be finite")
        if any(lo >= hi for lo, hi in zip(self.minimum, self.maximum)):
            raise ValueError("each boundary minimum must be smaller than its maximum")


class RadiaBatchFieldMap:
    """Vectorized Xsuite field-map callable backed by ``radia.Fld``.

    ``evaluator`` is injectable so the adapter and its unit contract can be
    tested without constructing a native Radia object.  Production callers
    normally leave it unset, in which case ``radia.Fld`` is imported lazily.
    """

    def __init__(self, radia_object: int, evaluator: FieldEvaluator | None = None):
        self.radia_object = int(radia_object)
        self._evaluator = evaluator

    def _field_evaluator(self) -> FieldEvaluator:
        if self._evaluator is not None:
            return self._evaluator
        import radia

        return radia.Fld

    def __call__(self, x: object, y: object, z: object) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x_arr, y_arr, z_arr = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        points = np.column_stack((x_arr.ravel(), y_arr.ravel(), z_arr.ravel()))
        values = np.asarray(
            self._field_evaluator()(self.radia_object, "b", points),
            dtype=float,
        )
        if values.shape == (3,) and points.shape[0] == 1:
            values = values.reshape(1, 3)
        if values.shape != (points.shape[0], 3):
            raise ValueError(
                "Radia field evaluator must return one three-component B vector per point"
            )
        if not np.all(np.isfinite(values)):
            raise ValueError("Radia field evaluator returned a non-finite magnetic field")
        shape = x_arr.shape
        return tuple(values[:, component].reshape(shape) for component in range(3))  # type: ignore[return-value]


def build_spatial_boris_integrator(
    fieldmap_callable: FieldMap,
    *,
    s_start_m: float,
    s_end_m: float,
    n_steps: int,
    log_trajectories: bool = True,
) -> Any:
    """Build an Xsuite spatial Boris integrator for a magnetic field map."""
    if not callable(fieldmap_callable):
        raise TypeError("fieldmap_callable must be callable")
    if not np.isfinite((s_start_m, s_end_m)).all() or s_end_m <= s_start_m:
        raise ValueError("s_start_m and s_end_m must be finite and increasing")
    if isinstance(n_steps, bool) or int(n_steps) != n_steps or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    try:
        import xtrack as xt
    except ImportError as exc:
        raise ImportError(
            "Xsuite tracking is optional; install Radia with the 'beam' extra"
        ) from exc

    integrator = xt.BorisSpatialIntegrator(
        fieldmap_callable=fieldmap_callable,
        s_start=float(s_start_m),
        s_end=float(s_end_m),
        n_steps=int(n_steps),
    )
    integrator.log_trajectories = bool(log_trajectories)
    return integrator


def first_box_exit_events(
    x_log: object,
    y_log: object,
    z_log: object,
    boundary: AxisAlignedBox,
) -> list[dict[str, object]]:
    """Return the first logged AABB exit for every particle that leaves it."""
    x = np.asarray(x_log, dtype=float)
    y = np.asarray(y_log, dtype=float)
    z = np.asarray(z_log, dtype=float)
    if x.shape != y.shape or x.shape != z.shape or x.ndim != 2:
        raise ValueError("trajectory logs must share shape (step, particle)")

    lo = np.asarray(boundary.minimum, dtype=float)
    hi = np.asarray(boundary.maximum, dtype=float)
    positions = np.stack((x, y, z), axis=-1)
    outside = np.any((positions < lo) | (positions > hi), axis=-1)
    events: list[dict[str, object]] = []
    for particle_index in range(outside.shape[1]):
        steps = np.flatnonzero(outside[:, particle_index])
        if steps.size:
            step = int(steps[0])
            events.append(
                {
                    "particle_index": particle_index,
                    "step_index": step,
                    "position_m": positions[step, particle_index].tolist(),
                }
            )
    return events


def track_magnetic_fieldmap(
    particles: Any,
    fieldmap_callable: FieldMap,
    *,
    s_start_m: float,
    s_end_m: float,
    n_steps: int,
    boundary: AxisAlignedBox | None = None,
) -> dict[str, object]:
    """Track Xsuite particles and return a JSON-friendly verification record.

    The supplied ``particles`` object is advanced in place by Xsuite.  Logged
    trajectories and optional AABB exit events are returned for validation and
    MCP artifact construction.
    """
    integrator = build_spatial_boris_integrator(
        fieldmap_callable,
        s_start_m=s_start_m,
        s_end_m=s_end_m,
        n_steps=n_steps,
        log_trajectories=True,
    )
    delta_before = np.asarray(particles.delta, dtype=float).copy()
    integrator.track(particles)
    x_log = np.asarray(integrator.x_log, dtype=float)
    y_log = np.asarray(integrator.y_log, dtype=float)
    z_log = np.asarray(integrator.z_log, dtype=float)
    delta_after = np.asarray(particles.delta, dtype=float)
    events = first_box_exit_events(x_log, y_log, z_log, boundary) if boundary else []
    return {
        "schema": "radia-xsuite-magnetic-track/v1",
        "backend": "xtrack.BorisSpatialIntegrator",
        "units": {"position": "m", "magnetic_flux_density": "T"},
        "particle_count": int(x_log.shape[1]),
        "step_count": int(x_log.shape[0]),
        "s_start_m": float(s_start_m),
        "s_end_m": float(s_end_m),
        "trajectory": {
            "x_m": x_log.tolist(),
            "y_m": y_log.tolist(),
            "z_m": z_log.tolist(),
        },
        "relative_momentum_deviation_before": delta_before.tolist(),
        "relative_momentum_deviation_after": delta_after.tolist(),
        "boundary_exit_events": events,
        "limitations": [
            "magnetic_field_only",
            "particle_matter_interactions_require_xcoll",
        ],
    }


def radia_magnetic_fieldmap(
    radia_object: int,
    evaluator: FieldEvaluator | None = None,
) -> RadiaBatchFieldMap:
    """Construct a vectorized Xsuite field-map callable from a Radia object."""
    return RadiaBatchFieldMap(radia_object, evaluator=evaluator)

