"""Major / minor axis of the time-locus of an AC vector quantity.

Part 5, eq 29-37 of the Wakao-Igarashi-Fujiwara-Kameari series.

In linear steady-state eddy current analysis a Cartesian component
``b_i(t)`` of a vector quantity (B, J, ...) is a phase-shifted sinusoid

    b_i(t) = Re(B_i exp(j omega t)),    i = x, y, z

with complex phasors ``B_i = B_{i,m} exp(j theta_i)``. The tip of
``b(t)`` traces a planar ellipse (3D locus of a single Fourier mode).
This module computes the major axis ``B_max`` and minor axis ``B_min``
of that ellipse.

Closed form (eq 37)
-------------------
Let

    s2  = |B_x|^2 + |B_y|^2 + |B_z|^2          (real, time-average * 2)
    z2  = B_x^2 + B_y^2 + B_z^2                (complex; complex dot of B with itself)

Then

    B_max = sqrt((s2 + |z2|) / 2)
    B_min = sqrt((s2 - |z2|) / 2)

Limit cases
-----------
* ``B_min = 0`` : linearly polarised (z2 has the same magnitude as s2;
  the locus collapses to a line).
* ``B_max = B_min`` : circularly polarised in some plane (``z2 = 0``).

Why this matters
----------------
The peak field over a period is what drives core loss when the magnetic
circuit is dominated by a near-sinusoidal flux waveform; ``B_max`` is
the right value for hysteresis-loss correlations. In induction heating
the peak ``J_max`` is what bounds local heating-density uniformisation
constraints.
"""

from __future__ import annotations

import numpy as np


def ac_locus_axes(B: complex | np.ndarray) -> tuple[float, float]:
    """Major and minor axis of a single 3-component AC phasor.

    Parameters
    ----------
    B : array_like, shape (3,), complex
        Complex phasor of a vector quantity.

    Returns
    -------
    B_max, B_min : float
        Major and minor axis of the time-locus ellipse, in the same
        units as the input.
    """
    B = np.asarray(B, dtype=complex)
    if B.shape != (3,):
        raise ValueError(f"B must have shape (3,), got {B.shape}")
    s2 = float(np.vdot(B, B).real)         # |Bx|^2 + |By|^2 + |Bz|^2
    z2 = complex(np.dot(B, B))              # Bx^2 + By^2 + Bz^2 (complex)
    abs_z2 = abs(z2)
    # Numerical safety: half-difference can be slightly negative on
    # circularly-polarised input due to rounding.
    arg_max = max(0.0, 0.5 * (s2 + abs_z2))
    arg_min = max(0.0, 0.5 * (s2 - abs_z2))
    return np.sqrt(arg_max), np.sqrt(arg_min)


def ac_locus_axes_batch(B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised major / minor axis over a batch of phasors.

    Parameters
    ----------
    B : ndarray, shape (..., 3), complex
        Last axis is the Cartesian component. Leading axes are batch
        dimensions (e.g. mesh nodes, time samples).

    Returns
    -------
    B_max, B_min : ndarray, shape (...,)
        Real-valued major / minor axes.
    """
    B = np.asarray(B, dtype=complex)
    if B.shape[-1] != 3:
        raise ValueError(f"last axis of B must be 3, got shape {B.shape}")
    s2 = np.real(np.einsum("...i,...i->...", np.conj(B), B))
    z2 = np.einsum("...i,...i->...", B, B)
    abs_z2 = np.abs(z2)
    arg_max = np.maximum(0.0, 0.5 * (s2 + abs_z2))
    arg_min = np.maximum(0.0, 0.5 * (s2 - abs_z2))
    return np.sqrt(arg_max), np.sqrt(arg_min)
