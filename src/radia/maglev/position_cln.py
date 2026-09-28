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

from .position_force import strictly_increasing as _strictly_increasing


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


__all__ = ["MovingHCurlCLNFamily"]
