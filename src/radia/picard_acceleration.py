"""Constrained Anderson acceleration for the real-valued material Picard loops.

The FEM engines of the static-electromagnet validations (mixed total/reduced Omega and
reduced-A) update a per-element permeability or reluctivity from the shared monotone B(H)
law with a fixed under-relaxation.  That damped Picard map contracts only slowly on a
saturating yoke: the slowest modes decay like ``1 - relaxation * (1 - L)`` per iteration,
so a step-size criterion of 2e-5 needs well over a hundred factorizations.

:class:`ConstrainedAndersonAccelerator` is the type-II Anderson mixing of Walker and Ni
(SIAM J. Numer. Anal. 49 (4), 2011) restricted to what a material coefficient can be:

* real arithmetic only -- the complex accelerator of :mod:`radia.esim_anderson` serves the
  surface-impedance Karl loop and would emit complex or negative permeabilities here;
* every accelerated iterate is projected onto the admissible interval of the B(H) law
  (``mu_r`` between one and the initial permeability, or the reluctivity range), so the
  next linear solve always sees a physical material;
* the extrapolation can be formed in ``log`` space (``transform="log"``): a shielded
  yoke spans three decades of secant permeability, and a multiplicative correction stays
  inside the interval where an additive one is clipped every step;
* an accelerated iterate is accepted only if the residual measured at it does not exceed
  ``accept_growth`` times the residual of the last accepted iterate.  A rejected iterate
  is replaced by the plain damped Picard step from the accepted one and the mixing
  history restarts, so a bad extrapolation costs one linear solve instead of an
  oscillation.  The history also restarts on a residual growth beyond ``restart_growth``,
  and the correction is clipped to ``step_clip`` damped-Picard steps;
* the residual history is kept so the caller can record it and estimate the contraction
  rate instead of guessing an iteration cap.

``depth=0`` is exactly the damped Picard update ``relaxation * g + (1 - relaxation) * x``
(projected), which keeps the existing engine contracts byte-identical when acceleration
is off.  The base Picard iterate is always formed in the original space; the transform only
shapes the Anderson correction.
"""
from __future__ import annotations

from collections import deque

import numpy as np

_TRANSFORMS = ("linear", "log")


class ConstrainedAndersonAccelerator:
    """Real-valued, projected, safeguarded Anderson type-II mixing for a map ``g(x)``."""

    def __init__(self, depth: int = 2, relaxation: float = 0.3, *, lower=None, upper=None,
                 step_clip: float = 100.0, restart_growth: float = 2.0,
                 transform: str = "linear", accept_growth: float = 1.0):
        depth = int(depth)
        if depth < 0:
            raise ValueError(f"depth must be >= 0, got {depth}")
        if not 0.0 < float(relaxation) <= 1.0:
            raise ValueError(f"relaxation must lie in (0, 1], got {relaxation}")
        if float(step_clip) <= 0.0:
            raise ValueError(f"step_clip must be positive, got {step_clip}")
        if float(restart_growth) <= 1.0:
            raise ValueError(f"restart_growth must exceed 1, got {restart_growth}")
        if float(accept_growth) <= 0.0:
            raise ValueError(f"accept_growth must be positive, got {accept_growth}")
        if transform not in _TRANSFORMS:
            raise ValueError(f"transform must be one of {_TRANSFORMS}, got {transform!r}")
        lower_value = None if lower is None else float(lower)
        upper_value = None if upper is None else float(upper)
        if lower_value is not None and upper_value is not None and lower_value > upper_value:
            raise ValueError(f"lower {lower_value} exceeds upper {upper_value}")
        if transform == "log" and (lower_value is None or lower_value <= 0.0):
            raise ValueError("transform='log' requires a positive lower bound")
        self.depth = depth
        self.relaxation = float(relaxation)
        self.lower = lower_value
        self.upper = upper_value
        self.step_clip = float(step_clip)
        self.restart_growth = float(restart_growth)
        self.transform = transform
        self.accept_growth = float(accept_growth)
        self._iterates: deque[np.ndarray] = deque(maxlen=depth + 1)
        self._residuals: deque[np.ndarray] = deque(maxlen=depth + 1)
        self._previous_residual_norm = None
        self._accepted = None          # (x, g, residual_norm) of the last accepted iterate
        self._candidate = None         # the accelerated iterate returned last, awaiting its residual
        self.residual_history: list[float] = []
        self.restarts = 0
        self.rejections = 0
        self.clips = 0
        self.projections = 0
        self.accelerated_steps = 0

    # ------------------------------------------------------------------ helpers
    def reset(self) -> None:
        """Forget the mixing history (the next step is a plain damped Picard step)."""
        self._iterates.clear()
        self._residuals.clear()
        self._previous_residual_norm = None
        self._candidate = None

    def project(self, values):
        """Clip ``values`` to the admissible interval, counting whether anything moved."""
        array = np.asarray(values, dtype=float)
        projected = np.clip(array, self.lower, self.upper)
        if np.any(projected != array):
            self.projections += 1
        return projected

    def _forward(self, values):
        if self.transform == "log":
            return np.log(np.maximum(values, self.lower))
        return np.asarray(values, dtype=float)

    def _inverse(self, values):
        if self.transform == "log":
            return np.exp(values)
        return np.asarray(values, dtype=float)

    def _picard(self, x, g):
        return self.relaxation * g + (1.0 - self.relaxation) * x

    # ------------------------------------------------------------------ the step
    def step(self, current, target):
        """Return the next iterate from ``x_k`` and ``g(x_k)``.

        ``target`` is the fixed-point map value (the material coefficient the linear
        solve just produced); ``current`` the coefficient that solve used.
        """
        x = np.asarray(current, dtype=float).reshape(-1)
        g = np.asarray(target, dtype=float).reshape(-1)
        if x.shape != g.shape:
            raise ValueError(f"current shape {x.shape} != target shape {g.shape}")
        if not (np.all(np.isfinite(x)) and np.all(np.isfinite(g))):
            raise ValueError("Anderson step received a non-finite iterate or map value")
        residual = g - x
        residual_norm = float(np.max(np.abs(residual)))
        self.residual_history.append(residual_norm)

        # A-posteriori safeguard: the residual just measured at the accelerated
        # candidate decides whether it is kept.  A worse residual reverts to the
        # damped Picard step from the last accepted iterate.
        if (self._candidate is not None and self._accepted is not None
                and np.array_equal(x, self._candidate)
                and residual_norm > self.accept_growth * self._accepted[2]):
            accepted_x, accepted_g, accepted_norm = self._accepted
            self.rejections += 1
            self.reset()
            self._previous_residual_norm = accepted_norm
            return self.project(self._picard(accepted_x, accepted_g))
        self._candidate = None
        self._accepted = (x.copy(), g.copy(), residual_norm)

        if (self._previous_residual_norm is not None and self._iterates
                and residual_norm > self.restart_growth * self._previous_residual_norm):
            self._iterates.clear()
            self._residuals.clear()
            self.restarts += 1
        self._previous_residual_norm = residual_norm

        u = self._forward(x)
        f = self._forward(g) - u
        self._iterates.append(u.copy())
        self._residuals.append(f.copy())

        # The damped Picard iterate is formed as the convex combination the engines
        # used before acceleration existed, so depth 0 reproduces them bit for bit.
        picard_next = self._picard(x, g)
        history = len(self._iterates) - 1
        if self.depth == 0 or history == 0:
            return self.project(picard_next)
        iterates = list(self._iterates)
        residuals = list(self._residuals)
        dU = np.column_stack([iterates[i + 1] - iterates[i] for i in range(history)])
        dF = np.column_stack([residuals[i + 1] - residuals[i] for i in range(history)])
        gamma, *_ = np.linalg.lstsq(dF, f, rcond=None)
        correction = -(dU + self.relaxation * dF) @ gamma
        limit = self.step_clip * max(float(np.linalg.norm(self.relaxation * f)), 1.0e-300)
        norm = float(np.linalg.norm(correction))
        if norm > limit:
            correction *= limit / norm
            self.clips += 1
        self.accelerated_steps += 1
        candidate = self.project(self._inverse(self._forward(picard_next) + correction))
        self._candidate = candidate.copy()
        return candidate

    def stats(self) -> dict[str, object]:
        """Counters and the residual history for the engine's ``nonlinear_stats``."""
        return {
            "anderson_depth": int(self.depth),
            "relaxation": float(self.relaxation),
            "transform": self.transform,
            "lower": self.lower,
            "upper": self.upper,
            "accept_growth": float(self.accept_growth),
            "accelerated_steps": int(self.accelerated_steps),
            "rejections": int(self.rejections),
            "restarts": int(self.restarts),
            "clips": int(self.clips),
            "projections": int(self.projections),
            "residual_history": [float(value) for value in self.residual_history],
        }


def estimate_contraction_rate(history, window: int = 6):
    """Geometric-mean ratio of the last ``window`` consecutive residuals, or None.

    The rate turns a step-size stopping criterion into an error estimate
    (``error ~ step * rate / (1 - rate)``) and predicts how many more iterations a
    tolerance needs; both belong in the engine's diagnostics rather than in a guess.
    """
    values = [float(v) for v in history if np.isfinite(v) and v > 0.0]
    if len(values) < 3:
        return None
    tail = values[-(int(window) + 1):]
    ratios = [tail[i + 1] / tail[i] for i in range(len(tail) - 1)]
    rate = float(np.exp(np.mean(np.log(ratios))))
    return rate if np.isfinite(rate) else None


__all__ = ["ConstrainedAndersonAccelerator", "estimate_contraction_rate"]
