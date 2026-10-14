"""Position-dependent force curves for moving conductive bodies."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def strictly_increasing(values: np.ndarray, name: str) -> None:
    """Raise unless ``values`` is a finite, strictly increasing 1-D sample."""

    if values.ndim != 1 or values.size < 2:
        raise ValueError(f"{name} must contain at least two values")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")


@dataclass(frozen=True)
class PositionForceCurve:
    """One-dimensional force curve with deterministic interpolation checks."""

    positions_m: np.ndarray
    force_N: np.ndarray
    name: str = "force"

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_m, dtype=float)
        forces = np.asarray(self.force_N, dtype=float)
        strictly_increasing(positions, "positions_m")
        if forces.shape != positions.shape:
            raise ValueError("force_N must match positions_m")
        if not np.all(np.isfinite(forces)):
            raise ValueError("force_N must contain only finite values")
        if not self.name:
            raise ValueError("name must not be empty")
        object.__setattr__(self, "positions_m", positions)
        object.__setattr__(self, "force_N", forces)

    def at(self, position_m):
        """Linearly interpolate force without extrapolation."""

        positions = np.asarray(position_m, dtype=float)
        if np.any(~np.isfinite(positions)):
            raise ValueError("position_m must be finite")
        if np.any(positions < self.positions_m[0]) or np.any(positions > self.positions_m[-1]):
            raise ValueError("position_m is outside the force-curve range")
        values = np.interp(positions, self.positions_m, self.force_N)
        return float(values) if values.ndim == 0 else values

    def force_result_at(
        self,
        position_m: float,
        *,
        lift_axis: int = 2,
        method: str = "interpolated_lorentz_force",
        frame: str = "global_cartesian",
    ) -> dict[str, object]:
        """Return one interpolated force as ``radia.force-result/v1``."""

        from radia.force import force_torque_result

        axis = int(lift_axis)
        if axis not in (0, 1, 2):
            raise ValueError("lift_axis must be 0, 1, or 2")
        position = float(position_m)
        force = [0.0, 0.0, 0.0]
        force[axis] = self.at(position)
        result = force_torque_result(
            force,
            None,
            method=method,
            frame=frame,
        )
        result.update(
            {
                "position_m": position,
                "force_curve": self.name,
                "lift_axis": axis,
            }
        )
        return result

    def crossings(self, target_force_N: float) -> np.ndarray:
        """Return all linearly interpolated positions where force reaches target."""

        target = float(target_force_N)
        if not np.isfinite(target):
            raise ValueError("target_force_N must be finite")
        roots: list[float] = []
        residual = self.force_N - target
        for i in range(self.positions_m.size - 1):
            left, right = residual[i], residual[i + 1]
            if left == 0.0:
                roots.append(float(self.positions_m[i]))
            if left * right < 0.0:
                fraction = -left / (right - left)
                roots.append(
                    float(
                        self.positions_m[i]
                        + fraction * (self.positions_m[i + 1] - self.positions_m[i])
                    )
                )
        if residual[-1] == 0.0:
            roots.append(float(self.positions_m[-1]))
        return np.asarray(roots, dtype=float)

    def compare(self, reference: PositionForceCurve) -> dict[str, float | int | str]:
        """Compare this curve to a reference on this curve's sample positions."""

        if not isinstance(reference, PositionForceCurve):
            raise TypeError("reference must be a PositionForceCurve")
        reference_force = np.asarray(reference.at(self.positions_m))
        error = self.force_N - reference_force
        scale = max(float(np.max(np.abs(reference_force))), np.finfo(float).tiny)
        return {
            "candidate": self.name,
            "reference": reference.name,
            "sample_count": int(self.positions_m.size),
            "max_abs_error_N": float(np.max(np.abs(error))),
            "rms_error_N": float(np.sqrt(np.mean(error**2))),
            "max_abs_error_normalized": float(np.max(np.abs(error)) / scale),
        }


__all__ = ["PositionForceCurve"]
