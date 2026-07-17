"""Position-dependent CLN coupling for moving conductive bodies.

HCurl Eddy Bubble reduces the spatial current space.  CLN/EVRS then carries
the reduced passive dynamics.  This module supplies the constant-basis motion
contract: models at neighboring positions may be interpolated only when their
state and port coordinates are identical.  Convex interpolation of the
Hermitian positive-semidefinite ``R``, ``L``, and surface blocks preserves
passivity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from radia.vim import HCurlEddyCLNModel


def _strictly_increasing(values: np.ndarray, name: str) -> None:
    if values.ndim != 1 or values.size < 2:
        raise ValueError(f"{name} must contain at least two values")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")


@dataclass(frozen=True)
class MovingHCurlCLNFamily:
    """Constant-basis HCurl Eddy Bubble/CLN models sampled over position."""

    positions_m: np.ndarray
    models: tuple[HCurlEddyCLNModel, ...]

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_m, dtype=float)
        _strictly_increasing(positions, "positions_m")
        models = tuple(self.models)
        if len(models) != positions.size:
            raise ValueError("models must match positions_m")
        if not all(isinstance(model, HCurlEddyCLNModel) for model in models):
            raise TypeError("models must contain HCurlEddyCLNModel objects")
        reference = models[0]
        for model in models:
            if not model.diagnostics()["passive"]:
                raise ValueError("constant-basis models must be passive")
        for model in models[1:]:
            if model.state_order != reference.state_order:
                raise ValueError("constant-basis models must have the same state order")
            if model.port_count != reference.port_count:
                raise ValueError("constant-basis models must have the same port count")
            if model.basis_names != reference.basis_names:
                raise ValueError("constant-basis models must have identical basis names")
            if model.blocks != reference.blocks:
                raise ValueError("constant-basis models must have identical block layout")
        object.__setattr__(self, "positions_m", positions)
        object.__setattr__(self, "models", models)

    @property
    def state_order(self) -> int:
        return self.models[0].state_order

    @property
    def port_count(self) -> int:
        return self.models[0].port_count

    def bracket(self, position_m: float) -> tuple[int, int, float]:
        """Return lower/upper sample indices and convex interpolation weight."""

        position = float(position_m)
        if not np.isfinite(position):
            raise ValueError("position_m must be finite")
        if position < self.positions_m[0] or position > self.positions_m[-1]:
            raise ValueError("position_m is outside the sampled CLN range")
        upper = int(np.searchsorted(self.positions_m, position, side="right"))
        if upper == 0:
            return 0, 0, 0.0
        if upper == self.positions_m.size:
            last = self.positions_m.size - 1
            return last, last, 0.0
        lower = upper - 1
        if position == self.positions_m[lower]:
            return lower, lower, 0.0
        span = self.positions_m[upper] - self.positions_m[lower]
        weight = (position - self.positions_m[lower]) / span
        return lower, upper, float(weight)

    def at(self, position_m: float) -> HCurlEddyCLNModel:
        """Interpolate a passive reduced model at one conductor position."""

        lower, upper, weight = self.bracket(position_m)
        if lower == upper:
            return self.models[lower]
        left = self.models[lower]
        right = self.models[upper]

        def blend(a: np.ndarray, b: np.ndarray) -> np.ndarray:
            return (1.0 - weight) * a + weight * b

        return HCurlEddyCLNModel(
            resistance=blend(left.resistance, right.resistance),
            inductance=blend(left.inductance, right.inductance),
            surface_mass=blend(left.surface_mass, right.surface_mass),
            port_rhs=blend(left.port_rhs, right.port_rhs),
            basis_names=left.basis_names,
            blocks=left.blocks,
        )

    def diagnostics(self) -> dict[str, object]:
        rows = [model.diagnostics() for model in self.models]
        return {
            "position_samples": int(self.positions_m.size),
            "position_min_m": float(self.positions_m[0]),
            "position_max_m": float(self.positions_m[-1]),
            "state_order": self.state_order,
            "port_count": self.port_count,
            "constant_basis": True,
            "all_samples_passive": all(row["passive"] for row in rows),
            "all_samples_finite_rl": all(row["finite_rl_state_space"] for row in rows),
        }


@dataclass(frozen=True)
class PositionForceCurve:
    """One-dimensional force curve with deterministic interpolation checks."""

    positions_m: np.ndarray
    force_N: np.ndarray
    name: str = "force"

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_m, dtype=float)
        forces = np.asarray(self.force_N, dtype=float)
        _strictly_increasing(positions, "positions_m")
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

    def compare(self, reference: "PositionForceCurve") -> dict[str, float | int | str]:
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


__all__ = ["MovingHCurlCLNFamily", "PositionForceCurve"]
