"""The production soft-iron B(H) interpolant shared by every Radia route.

The table is interpolated through its magnetization M = B/mu0 - H with
monotone PCHIP, and B(H) = mu0 (H + M(H)) is returned as a piecewise cubic on
the table's own knots.  When the tabulated M is non-decreasing this keeps
dB/dH >= mu0 everywhere: PCHIP applied to B(H) only guarantees dB/dH >= 0 and
drops below mu0 between the rows of real saturating tables (0.69 mu0 on ESRF
example 3, 0 on example 7), i.e. its M falls where the data does not.  A table
whose own M falls (reported by the soft-iron contract) keeps PCHIP on B(H), the
only one of the two that is then still monotone in B.
"""
from __future__ import annotations

import numpy as np

MU_0 = 4e-7 * np.pi

MAGNETIZATION_PCHIP = "magnetization_pchip"
FLUX_DENSITY_PCHIP = "flux_density_pchip"


def _magnetization_rises(H, B):
    M = B / MU_0 - H
    scale = max(1.0, float(np.max(np.abs(M))))
    return bool(np.all(np.diff(M) >= -1e-9 * scale))


def bh_interpolation_kind(H, B):
    """Which monotone PCHIP the table takes (see the module docstring)."""
    H = np.asarray(H, dtype=float)
    B = np.asarray(B, dtype=float)
    return MAGNETIZATION_PCHIP if _magnetization_rises(H, B) else FLUX_DENSITY_PCHIP


def monotone_bh_pchip(H, B, extrapolate=None):
    """B(H) through the table as a ``scipy.interpolate.PPoly`` on the table knots.

    Drop-in for ``PchipInterpolator(H, B, extrapolate=...)``: it has the same
    breakpoints, coefficient layout ``c[k, i]`` (power of ``H - H[i]``), and
    ``derivative``/``antiderivative``.  The caller owns the extension beyond
    the table.
    """
    from scipy.interpolate import PchipInterpolator, PPoly

    H = np.asarray(H, dtype=float)
    B = np.asarray(B, dtype=float)
    if bh_interpolation_kind(H, B) == FLUX_DENSITY_PCHIP:
        return PchipInterpolator(H, B, extrapolate=extrapolate)
    magnetization = PchipInterpolator(H, B / MU_0 - H, extrapolate=extrapolate)
    coefficients = MU_0 * np.array(magnetization.c, dtype=float)
    coefficients[-2] += MU_0                  # d/dH of mu0 H on every piece
    coefficients[-1] += MU_0 * H[:-1]         # mu0 H at each piece's left knot
    return PPoly(coefficients, H.copy(), extrapolate=extrapolate)
