"""Position-dependent Foster models for moving conductive bodies.

HCurl Eddy Bubble reduces the spatial current space; the Foster modal form
carries the reduced passive dynamics.  For a conductor meshed once and moved
relative to its source, the reduced ``R``/``L`` pair is shared by every
position and only the source port matrix changes.  This module enforces that
contract: one set of Foster modes serves the whole family, and the port
matrix interpolates linearly between sampled positions, which keeps every
interpolated model passive.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from radia.vim import HCurlEddyFosterModel

from .position_force import strictly_increasing


@dataclass(frozen=True)
class MovingHCurlFosterFamily:
    """Foster models sampled over position with one shared ``R``/``L`` pair."""

    positions_m: np.ndarray
    models: tuple[HCurlEddyFosterModel, ...]
    shared_operator_rtol: float = 1.0e-12

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_m, dtype=float)
        strictly_increasing(positions, "positions_m")
        models = tuple(self.models)
        if len(models) != positions.size:
            raise ValueError("models must match positions_m")
        if not all(isinstance(model, HCurlEddyFosterModel) for model in models):
            raise TypeError("models must contain HCurlEddyFosterModel objects")
        reference = models[0]
        for model in models:
            if not model.diagnostics()["passive"]:
                raise ValueError("family models must be passive")
        for model in models[1:]:
            if model.port_count != reference.port_count:
                raise ValueError("family models must have the same port count")
            if model.basis_names != reference.basis_names:
                raise ValueError("family models must have identical basis names")
            if model.blocks != reference.blocks:
                raise ValueError("family models must have identical block layout")
            for name in ("resistance", "inductance", "surface_mass"):
                ref = getattr(reference, name)
                cur = getattr(model, name)
                if cur.shape != ref.shape:
                    raise ValueError("family models must have the same state order")
                scale = max(float(np.linalg.norm(ref)), np.finfo(float).tiny)
                if float(np.linalg.norm(cur - ref)) > self.shared_operator_rtol * scale:
                    raise ValueError(
                        f"family models differ in {name}; a moving Foster family needs "
                        "one shared R/L pair (fixed conductor mesh)")
        object.__setattr__(self, "positions_m", positions)
        object.__setattr__(self, "models", models)

    @property
    def state_order(self) -> int:
        return self.models[0].state_order

    @property
    def port_count(self) -> int:
        return self.models[0].port_count

    @property
    def decay_rates(self) -> np.ndarray:
        return self.models[0].decay_rates

    def bracket(self, position_m: float) -> tuple[int, int, float]:
        """Return lower/upper sample indices and convex interpolation weight."""

        position = float(position_m)
        if not np.isfinite(position):
            raise ValueError("position_m must be finite")
        if position < self.positions_m[0] or position > self.positions_m[-1]:
            raise ValueError("position_m is outside the sampled family range")
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
        return lower, upper, float((position - self.positions_m[lower]) / span)

    def at(self, position_m: float) -> HCurlEddyFosterModel:
        """Return the Foster model at one conductor position (port matrix interpolated)."""

        lower, upper, weight = self.bracket(position_m)
        if lower == upper:
            return self.models[lower]
        left = self.models[lower]
        right = self.models[upper]
        return HCurlEddyFosterModel(
            resistance=left.resistance,
            inductance=left.inductance,
            surface_mass=left.surface_mass,
            port_rhs=(1.0 - weight) * left.port_rhs + weight * right.port_rhs,
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
            "shared_modes": True,
            "all_samples_passive": all(row["passive"] for row in rows),
            "all_samples_finite_rl": all(row["finite_rl_state_space"] for row in rows),
        }


__all__ = ["MovingHCurlFosterFamily"]
