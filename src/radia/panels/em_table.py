"""
Lookup helpers for the 2D EM table produced by ``calc_em_table.py``.

Shared by the offline-table consistency test and the runtime heat
solver (``calc_heat_with_em_table.py``).  Keeps the interpolation
contract in one place so the test and the runtime cannot drift apart.

Usage::

    from em_table import load_em_table, interp_qsurf, interp_Zs

    tab = load_em_table("em_table_steel_50kHz.npz")
    q   = interp_qsurf(tab, H_t=5.0e4, T_celsius=350.0)
    Z   = interp_Zs(tab,    H_t=5.0e4, T_celsius=350.0)   # complex

Both H_t and T_celsius can be scalars OR NumPy arrays (broadcast-able);
the return mirrors the input shape.

Interpolation
-------------
Bilinear in (log H, T).  H is log-scaled because q_surf spans ~3-4
decades from low-current trim to full Curie-point operation; linear
H interp under-resolves the low-end and is wasteful at the high end.
Outside the grid the lookup CLAMPS (no extrapolation).  This matches
the sigma(T) clamping convention in calc_em_table and prevents the
heat solver from drifting into a physically-unsupported regime when
the user's I(t) trajectory accidentally pokes past the grid corner.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class EMTable:
    """In-memory representation of a calc_em_table .npz."""
    path: str
    H_grid: np.ndarray              # (n_H,), log-spaced, A/m
    T_grid: np.ndarray              # (n_T,), linear, degC
    Zs_re: np.ndarray               # (n_H, n_T), Ohm
    Zs_im: np.ndarray               # (n_H, n_T), Ohm
    q_surf: np.ndarray              # (n_H, n_T), W/m^2
    meta: dict

    @property
    def logH_grid(self) -> np.ndarray:
        return np.log(self.H_grid)

    @property
    def frequency(self) -> float:
        return float(self.meta.get("frequency", 0.0))

    @property
    def material(self) -> str:
        return str(self.meta.get("material", ""))


def load_em_table(path: str) -> EMTable:
    """Load a .npz produced by calc_em_table.py."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"EM table not found: {path}")
    d = np.load(path, allow_pickle=False)
    needed = ("H_grid", "T_grid", "Zs_re", "Zs_im", "q_surf", "meta")
    missing = [k for k in needed if k not in d.files]
    if missing:
        raise ValueError(
            f"EM table {path} is missing keys {missing}; "
            f"available: {d.files}")
    meta = json.loads(str(d["meta"]))
    H_grid = np.asarray(d["H_grid"], dtype=float)
    T_grid = np.asarray(d["T_grid"], dtype=float)
    if not np.all(np.diff(H_grid) > 0):
        raise ValueError(f"H_grid not strictly increasing in {path}")
    if not np.all(np.diff(T_grid) > 0):
        raise ValueError(f"T_grid not strictly increasing in {path}")
    return EMTable(
        path=os.path.abspath(path),
        H_grid=H_grid,
        T_grid=T_grid,
        Zs_re=np.asarray(d["Zs_re"], dtype=float),
        Zs_im=np.asarray(d["Zs_im"], dtype=float),
        q_surf=np.asarray(d["q_surf"], dtype=float),
        meta=meta,
    )


def _bilinear(tab_field: np.ndarray, logH_grid: np.ndarray,
              T_grid: np.ndarray, logH, T) -> np.ndarray:
    """Bilinear interp of ``tab_field`` on the (logH, T) grid.

    Clamps to grid edges outside the support.  Vectorised over logH/T.
    """
    logH = np.asarray(logH, dtype=float)
    T = np.asarray(T, dtype=float)
    logH_clip = np.clip(logH, logH_grid[0], logH_grid[-1])
    T_clip = np.clip(T, T_grid[0], T_grid[-1])

    # i,j index pairs bracketing the query
    i = np.searchsorted(logH_grid, logH_clip) - 1
    i = np.clip(i, 0, len(logH_grid) - 2)
    j = np.searchsorted(T_grid, T_clip) - 1
    j = np.clip(j, 0, len(T_grid) - 2)

    dH = logH_grid[i + 1] - logH_grid[i]
    dT = T_grid[j + 1] - T_grid[j]
    u = (logH_clip - logH_grid[i]) / dH
    v = (T_clip - T_grid[j]) / dT

    f00 = tab_field[i, j]
    f10 = tab_field[i + 1, j]
    f01 = tab_field[i, j + 1]
    f11 = tab_field[i + 1, j + 1]
    return ((1 - u) * (1 - v) * f00 + u * (1 - v) * f10
            + (1 - u) * v * f01 + u * v * f11)


def interp_qsurf(tab: EMTable, H_t, T_celsius) -> np.ndarray:
    """q_surf [W/m^2] at the given (|H_t|, T) points."""
    H, T = np.broadcast_arrays(
        np.abs(np.asarray(H_t, dtype=float)),
        np.asarray(T_celsius, dtype=float),
    )
    if not np.all(np.isfinite(H)) or not np.all(np.isfinite(T)):
        raise ValueError("H_t and T_celsius must be finite")
    zero_field = H == 0.0
    H_safe = np.where(zero_field, tab.H_grid[0], H)
    q_surf = _bilinear(
        tab.q_surf, tab.logH_grid, tab.T_grid, np.log(H_safe), T
    )
    return np.where(zero_field, 0.0, q_surf)


def interp_Zs(tab: EMTable, H_t, T_celsius) -> np.ndarray:
    """Complex Z_s [Ohm] at the given (|H_t|, T) points."""
    H, T = np.broadcast_arrays(
        np.abs(np.asarray(H_t, dtype=float)),
        np.asarray(T_celsius, dtype=float),
    )
    if not np.all(np.isfinite(H)) or not np.all(np.isfinite(T)):
        raise ValueError("H_t and T_celsius must be finite")
    H_safe = np.where(H > 0.0, H, tab.H_grid[0])
    logH = np.log(H_safe)
    Zr = _bilinear(tab.Zs_re, tab.logH_grid, tab.T_grid, logH, T)
    Zi = _bilinear(tab.Zs_im, tab.logH_grid, tab.T_grid, logH, T)
    return Zr + 1j * Zi

@dataclass(frozen=True)
class SurfaceImpedanceModel:
    """User-designed Zs(|Ht|, temperature) at one declared frequency.

    ``evaluate`` accepts broadcast arrays in A/m (peak) and degC and returns
    complex ohms with the same shape, or a scalar for a uniform model.
    Bounds are inclusive. No clipping, extrapolation or fitted correction is
    performed. This material interface does not freeze the electromagnetic
    field: a coupled solver must reevaluate the field after a material update.
    """

    evaluate: object
    frequency_hz: float
    h_bounds: tuple[float, float]
    temperature_bounds_c: tuple[float, float]
    name: str

    def __post_init__(self):
        if not callable(self.evaluate):
            raise TypeError("evaluate must be a callable Zs(H_peak_A_per_m, T_celsius)")
        if not np.isfinite(self.frequency_hz) or self.frequency_hz <= 0:
            raise ValueError("frequency_hz must be finite and positive")
        for label, bounds in (("H", self.h_bounds), ("temperature", self.temperature_bounds_c)):
            a = np.asarray(bounds, dtype=float)
            if a.shape != (2,) or not np.all(np.isfinite(a)) or a[0] >= a[1]:
                raise ValueError(f"{label} bounds must be two finite increasing values")
        if self.h_bounds[0] < 0 or self.temperature_bounds_c[0] < -273.15:
            raise ValueError("H must be non-negative and temperature above absolute zero")
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("An identifiable model name is required")

    def impedance(self, H_peak_A_per_m, T_celsius, *, frequency_hz):
        """Evaluate the model, rejecting invalid inputs and non-passive outputs."""
        if not np.isfinite(frequency_hz) or frequency_hz != self.frequency_hz:
            raise ValueError("Requested frequency does not match the impedance model")
        H, T = np.broadcast_arrays(np.asarray(H_peak_A_per_m, dtype=float),
                                   np.asarray(T_celsius, dtype=float))
        for label, value, bounds in (("H", H, self.h_bounds),
                                     ("temperature", T, self.temperature_bounds_c)):
            if (not np.all(np.isfinite(value)) or np.any(value < bounds[0])
                    or np.any(value > bounds[1])):
                raise ValueError(f"{label} outside the declared model domain {bounds}")
        value = np.asarray(self.evaluate(H.copy(), T.copy()), dtype=complex)
        if value.ndim == 0:
            value = np.full(H.shape, value, dtype=complex)
        if value.shape != H.shape:
            raise ValueError(f"Zs shape {value.shape} does not match input shape {H.shape}")
        if not np.all(np.isfinite(value)) or np.any(value.real < 0):
            raise ValueError("Zs must be finite and passive (real part >= 0)")
        return value.copy()

    def heat_flux(self, H_peak_A_per_m, T_celsius, *, frequency_hz):
        """Cycle-average local heat flux [W/m^2], for peak phasors."""
        z = self.impedance(H_peak_A_per_m, T_celsius, frequency_hz=frequency_hz)
        with np.errstate(over="ignore", invalid="ignore"):
            q = 0.5 * z.real * np.square(np.asarray(H_peak_A_per_m, dtype=float))
        if not np.all(np.isfinite(q)):
            raise ValueError("Surface heat flux overflowed")
        return q

    @classmethod
    def from_table(cls, table: EMTable, *, name: str):
        """Use existing EMTable Zs interpolation with strict range validation.

        q_surf is deliberately not interpolated independently: heat is derived
        from the evaluated complex Zs and the actual field, preserving power.
        The arrays are copied so later edits to the source table cannot alter
        a running model. A log-H table supports only strictly positive H.
        """
        h, t = np.array(table.H_grid, dtype=float), np.array(table.T_grid, dtype=float)
        for label, a in (("H", h), ("temperature", t)):
            if (a.ndim != 1 or len(a) < 2 or not np.all(np.isfinite(a))
                    or np.any(np.diff(a) <= 0)):
                raise ValueError(f"Invalid {label} table grid")
        if h[0] <= 0:
            raise ValueError("Log-H table requires positive field nodes")
        z = np.array(table.Zs_re, dtype=float) + 1j * np.array(table.Zs_im, dtype=float)
        if (np.shape(table.Zs_re) != (len(h), len(t))
                or np.shape(table.Zs_im) != (len(h), len(t))
                or not np.all(np.isfinite(z)) or np.any(z.real < 0)):
            raise ValueError("Invalid impedance table shape, finite values or passivity")
        copied = EMTable("memory", h, t, z.real.copy(), z.imag.copy(),
                         np.zeros_like(z.real), dict(table.meta))
        return cls(lambda H, T: interp_Zs(copied, H, T), table.frequency,
                   (float(h[0]), float(h[-1])), (float(t[0]), float(t[-1])), name)
