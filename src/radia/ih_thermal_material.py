"""Temperature-dependent thermal material for induction-heating solves.

A :class:`ThermalMaterial` gives the conductivity ``k(T)`` and the
volumetric enthalpy ``H(T) = rho * h(T)`` of a workpiece.  ``h(T)`` is the
integral of the tabulated specific heat plus an optional latent heat
released uniformly over a temperature interval (the solid-state
transformation of steel, or melting), so energy is conserved exactly by an
enthalpy time integrator whatever the shape of ``cp(T)``.

Tables come from the user: this module ships no alloy data.  A CSV has a
header row and the columns ``T_C, k_W_mK, cp_J_kgK``; the density is a
single value because the thermal mesh does not move.  Temperatures outside
the table are an error unless extrapolation is explicitly allowed, in
which case the end values are held and the run records that it happened.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

import numpy as np

# Enthalpy is tabulated on this grid spacing [K] for the NGSolve splines.
_ENTHALPY_STEP_K = 1.0


@dataclass
class ThermalMaterial:
    rho: float                       # kg/m^3
    T: np.ndarray                    # degC, strictly increasing
    k: np.ndarray                    # W/(m K)
    cp: np.ndarray                   # J/(kg K)
    latent_heat: float = 0.0         # J/kg
    latent_range: tuple = (0.0, 0.0)  # degC interval of the latent release
    allow_extrapolation: bool = False
    source: str = "table"
    _grid: np.ndarray = field(default=None, repr=False)
    _H: np.ndarray = field(default=None, repr=False)

    # -- construction ---------------------------------------------------
    @classmethod
    def constant(cls, rho, cp, k):
        return cls(rho=float(rho), T=np.array([-273.15, 1.0e4]),
                   k=np.array([float(k)] * 2), cp=np.array([float(cp)] * 2),
                   allow_extrapolation=False, source="constant")

    @classmethod
    def from_csv(cls, path, *, rho, latent_heat=0.0, latent_range=None,
                 allow_extrapolation=False):
        with open(path, encoding="utf-8-sig", newline="") as handle:
            rows = [r for r in csv.reader(handle) if r and not
                    r[0].lstrip().startswith("#")]
        header = [h.strip().lower() for h in rows[0]]
        want = ("t_c", "k_w_mk", "cp_j_kgk")
        if tuple(header[:3]) != want:
            raise ValueError(f"{path}: header must start with "
                             f"{','.join(want)}, got {rows[0]}")
        data = np.asarray([[float(v) for v in r[:3]] for r in rows[1:]])
        return cls(rho=float(rho), T=data[:, 0], k=data[:, 1], cp=data[:, 2],
                   latent_heat=float(latent_heat),
                   latent_range=tuple(latent_range or (0.0, 0.0)),
                   allow_extrapolation=bool(allow_extrapolation),
                   source=str(path))

    def __post_init__(self):
        self.T = np.asarray(self.T, float)
        self.k = np.asarray(self.k, float)
        self.cp = np.asarray(self.cp, float)
        if not (len(self.T) == len(self.k) == len(self.cp) >= 2):
            raise ValueError("material table needs at least two rows")
        if not np.all(np.diff(self.T) > 0):
            raise ValueError("material table temperatures must increase")
        if np.any(self.k <= 0) or np.any(self.cp <= 0) or self.rho <= 0:
            raise ValueError("rho, k and cp must be positive")
        if self.latent_heat < 0:
            raise ValueError("latent heat must be non-negative")
        if self.latent_heat > 0:
            t0, t1 = map(float, self.latent_range)
            if not t1 > t0:
                raise ValueError("latent heat needs a range T_start < T_end")
            if t0 < self.T[0] or t1 > self.T[-1]:
                raise ValueError("the latent range must lie inside the table")
        self._build_enthalpy()

    # -- properties -----------------------------------------------------
    @property
    def is_constant(self) -> bool:
        return (np.all(self.k == self.k[0]) and np.all(self.cp == self.cp[0])
                and self.latent_heat == 0.0)

    def _build_enthalpy(self):
        lo, hi = self.T[0], self.T[-1]
        if self.source == "constant" and self.is_constant:
            grid = np.array([lo, hi])
        else:
            n = int(math.ceil((hi - lo) / _ENTHALPY_STEP_K)) + 1
            grid = np.unique(np.r_[np.linspace(lo, hi, n), self.T,
                                   list(self.latent_range)
                                   if self.latent_heat else []])
            grid = grid[(grid >= lo) & (grid <= hi)]
        cp = np.interp(grid, self.T, self.cp)
        # trapezoid integral of cp (exact for the piecewise-linear table)
        h = np.r_[0.0, np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(grid))]
        if self.latent_heat:
            t0, t1 = self.latent_range
            h = h + self.latent_heat * np.clip((grid - t0) / (t1 - t0), 0, 1)
        self._grid = grid
        self._H = self.rho * h

    def conductivity(self, T):
        return np.interp(self._check(T), self.T, self.k)

    def enthalpy(self, T):
        """Volumetric enthalpy [J/m^3] relative to the table start.

        With extrapolation allowed, ``cp`` is held at its end values, so the
        enthalpy continues linearly beyond the table.
        """
        T = np.asarray(T, float)
        Tc = self._check(T)
        H = np.interp(Tc, self._grid, self._H)
        lo, hi = self._grid[0], self._grid[-1]
        c_lo, c_hi = self.rho * self.cp[0], self.rho * self.cp[-1]
        return H + c_hi * np.maximum(T - hi, 0.0) + c_lo * np.minimum(T - lo,
                                                                      0.0)

    def capacity(self, T):
        """dH/dT [J/(m^3 K)] (piecewise constant between grid points)."""
        T = np.asarray(self._check(T), float)
        i = np.clip(np.searchsorted(self._grid, T, side="right") - 1, 0,
                    len(self._grid) - 2)
        return (self._H[i + 1] - self._H[i]) / (self._grid[i + 1]
                                                - self._grid[i])

    def _check(self, T):
        T = np.asarray(T, float)
        if not self.allow_extrapolation and T.size:
            lo, hi = float(np.min(T)), float(np.max(T))
            if lo < self.T[0] - 1e-9 or hi > self.T[-1] + 1e-9:
                raise ValueError(
                    f"temperature {lo:.1f}..{hi:.1f} C leaves the material "
                    f"table range {self.T[0]:.1f}..{self.T[-1]:.1f} C; "
                    "extend the table or allow extrapolation explicitly")
        return np.clip(T, self.T[0], self.T[-1])

    # -- NGSolve coefficient functions -------------------------------------
    def coefficient_functions(self, u):
        """(k(u), H(u), dH/dT(u)) as NGSolve coefficient functions of ``u``.

        ``u`` is clamped to the table range inside the CFs; whether that
        clamp was ever active is checked on the nodal values by the solver.
        """
        from ngsolve import BSpline, IfPos

        lo, hi = float(self.T[0]), float(self.T[-1])
        uc = IfPos(u - hi, hi, IfPos(lo - u, lo, u))
        k_spl = _linear_spline(BSpline, self.T, self.k)
        H_spl = _linear_spline(BSpline, self._grid, self._H)
        c_vals = np.diff(self._H) / np.diff(self._grid)
        # NGSolve B-splines are half-open on their last interval; append
        # one extra knot so the clamped end value evaluates correctly.
        g = [float(t) for t in self._grid] + [float(self._grid[-1]) + 1.0]
        c_spl = BSpline(1, g, [float(c) for c in c_vals] + [float(c_vals[-1])])
        c_lo, c_hi = self.rho * self.cp[0], self.rho * self.cp[-1]
        H_cf = (H_spl(uc) + IfPos(u - hi, c_hi * (u - hi), 0.0)
                + IfPos(lo - u, c_lo * (u - lo), 0.0))
        return k_spl(uc), H_cf, c_spl(uc)

    def audit(self) -> dict:
        return {"source": self.source, "rho_kg_m3": self.rho,
                "T_range_C": [float(self.T[0]), float(self.T[-1])],
                "k_range_W_mK": [float(self.k.min()), float(self.k.max())],
                "cp_range_J_kgK": [float(self.cp.min()),
                                   float(self.cp.max())],
                "latent_heat_J_kg": self.latent_heat,
                "latent_range_C": [float(v) for v in self.latent_range],
                "allow_extrapolation": self.allow_extrapolation,
                "temperature_dependent": not self.is_constant}


def _linear_spline(BSpline, x, y):
    """NGSolve piecewise-linear interpolant through (x_i, y_i)."""
    x = [float(v) for v in x]
    y = [float(v) for v in y]
    # An order-2 B-spline with knots [x0, x0, x1, ..., xn] interpolates the
    # control values at the knots, but its last interval is half-open, so
    # the end point gets one extra knot carrying the end value.
    x = x + [x[-1] + 1.0]
    y = y + [y[-1]]
    return BSpline(2, [x[0]] + x + [x[-1]], y)
