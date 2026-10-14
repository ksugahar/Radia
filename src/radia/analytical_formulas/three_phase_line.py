"""Three-phase line magnetic field: direct evaluation + far-field asymptotes.

Part 4 §5 of the Wakao-Igarashi-Fujiwara-Kameari series (straight
arrangements: triangle, planar, hexagon) and Part 5 §3 (helical).

Conventions
-----------
Phasor amplitudes use the peak convention (eq 33 of Part 4):

    I_k(t) = Re(I_peak,k exp(j omega t))

If the user prefers an RMS-based input ``I_rms``, scale by ``sqrt(2)``
on the way in:

    I_peak = sqrt(2) * I_rms

For balanced three-phase currents,

    I_peak,k = I_peak * exp(- j 2 k pi / 3),    k = 0, 1, 2.

For a circularly-polarised far field the **time-constant** magnitude
``|B(t)|`` equals the RMS of ``|B(t)|`` (since ``|B(t)|`` is constant
in time); both are the value reported by the closed-form formulas
below, called consistently the "amplitude".

Functions
---------
``vector_potential_z(positions, currents, observation)``
    Exact A_z of an arbitrary set of straight currents (z-directed).

``field_xy(positions, currents, observation)``
    Exact ``(B_x, B_y)`` from the analytic gradient of A_z.

``triangle_far_field_amplitude(r, I_peak, a)``,
``planar_far_field_amplitude(r, I_peak, d)``
    Time-constant ``|B|`` at far observation ``r`` (eq 38-43 of Part 4).
    Both decay as ``1 / r**2``.

``helical_far_field_amplitude(r, I_peak, a, p)``
    Time-constant ``|B|`` at far ``r`` for a triangle-arrangement
    *helical* three-phase line of pitch ``p`` (Part 5, eq 28). Decay
    is exponential in ``r / p``.

Hexagonal arrangement
---------------------
A naive "two anti-phase three-phase sets at opposite vertices"
arrangement (built by :func:`hexagon_positions` plus
:func:`balanced_six_phase_currents`) gives ``1 / r**2`` decay -- the
opposite-position anti-phase wires reinforce, not cancel, the 2D
dipole moment. The 1/r**3 hexagonal asymptote in PDF eq 44 corresponds
to a specific phase-on-vertex assignment described in Fig 5(b) of the
PDF that needs separate construction; see the open issue note in the
test ``test_hexagon_naive_arrangement_is_dipolar`` for details.

References
----------
Part 4 (SA-04-1 / RM-04-1, 2004), eq 30-44.
Part 5 (SA-04-44 / RM-04-68, 2004), eq 11-28.
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4.0e-7 * math.pi
SQRT2 = math.sqrt(2.0)


# ---------------------------------------------------------------------------
# Direct evaluation: arbitrary set of straight currents
# ---------------------------------------------------------------------------


def vector_potential_z(
    positions,
    currents,
    observation,
):
    """Vector-potential A_z of arbitrary z-directed straight currents.

    Parameters
    ----------
    positions : array_like, shape (n_lines, 2)
        Cartesian (x, y) positions of the lines.
    currents : array_like of shape (n_lines,)
        Phasor (peak) currents (Amperes); accepts real or complex.
    observation : array_like, shape (..., 2)
        Observation point(s) in the (x, y) plane.

    Returns
    -------
    A_z : ndarray
        Same leading shape as ``observation``; phasor T m, with a
        position-independent constant dropped (irrelevant to ``B``).
    """
    positions = np.asarray(positions, dtype=float)
    currents = np.asarray(currents)
    obs = np.asarray(observation, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions must have shape (n_lines, 2)")
    if obs.shape[-1] != 2:
        raise ValueError("observation must have shape (..., 2)")
    dx = obs[..., None, 0] - positions[None, :, 0]
    dy = obs[..., None, 1] - positions[None, :, 1]
    r2 = dx * dx + dy * dy
    log_r2 = np.where(r2 > 0.0, np.log(np.where(r2 > 0.0, r2, 1.0)), 0.0)
    return -MU_0 / (4.0 * math.pi) * np.einsum(
        "...n,n->...", log_r2, currents
    )


def field_xy(
    positions,
    currents,
    observation,
) -> np.ndarray:
    """Magnetic flux density ``(B_x, B_y)`` from straight z-currents.

    Differentiates ``A_z`` analytically; output dtype is real for real
    currents and complex for complex currents.
    """
    positions = np.asarray(positions, dtype=float)
    currents = np.asarray(currents)
    obs = np.asarray(observation, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions must have shape (n_lines, 2)")
    if obs.shape[-1] != 2:
        raise ValueError("observation must have shape (..., 2)")

    dx = obs[..., None, 0] - positions[None, :, 0]
    dy = obs[..., None, 1] - positions[None, :, 1]
    r2 = dx * dx + dy * dy
    inv_r2 = np.where(r2 > 0.0, 1.0 / np.where(r2 > 0.0, r2, 1.0), 0.0)
    pref = MU_0 / (2.0 * math.pi)
    Bx = -pref * np.einsum("...n,...n,n->...", dy, inv_r2, currents)
    By = +pref * np.einsum("...n,...n,n->...", dx, inv_r2, currents)
    out = np.zeros(obs.shape, dtype=Bx.dtype)
    out[..., 0] = Bx
    out[..., 1] = By
    return out


# ---------------------------------------------------------------------------
# Far-field amplitudes: arrangement-specific closed forms
# ---------------------------------------------------------------------------


def triangle_far_field_amplitude(r: float, I_peak: float, a: float) -> float:
    """Time-constant ``|B|`` far from a triangle-arrangement three-phase
    line (Part 4 eq 38-39).

    Three lines on a circle of radius ``a`` at angles
    ``2 k pi / 3``, currents ``I_peak * exp(- j 2 k pi / 3)``.

        |B|(t) ~ 3 mu_0 a I_peak / (4 pi r**2)

    The far-field is circularly polarised, so ``|B|`` is constant in
    time and equals the RMS of ``|B|``.
    """
    return 3.0 * MU_0 * a * I_peak / (4.0 * math.pi * r * r)


def planar_far_field_amplitude(r: float, I_peak: float, d: float) -> float:
    """Time-locus axis ``|B|`` far from a coplanar three-phase line
    (Part 4 eq 41-43).

    Three lines on a row spaced ``d`` apart at ``-d, 0, +d``, currents
    ``I_peak * exp(- j 2 k pi / 3)``. Unlike the triangle case the
    far-field is linearly polarised; the major-axis amplitude is

        |B|_max ~ sqrt(3) mu_0 d I_peak / (2 pi r**2)

    (sqrt(3) coefficient comes from the off-diagonal phasor
    combination ``q* - q = j sqrt(3)`` of the three-phase set).
    """
    return math.sqrt(3.0) * MU_0 * d * I_peak / (2.0 * math.pi * r * r)


def helical_near_field_amplitude(
    r: float, I_peak: float, a: float, p: float
) -> float:
    """Time-constant ``|B|`` close to (but outside) a helical three-phase
    line, in the regime ``a << r << p`` (Part 5 eq 25-26).

    Approximated by the straight-triangle formula:

        |B| ~ 3 mu_0 a I_peak / (4 pi r**2).
    """
    if p <= 0.0:
        raise ValueError(f"pitch p must be > 0, got {p}")
    return triangle_far_field_amplitude(r, I_peak, a)


def helical_far_field_amplitude(
    r: float, I_peak: float, a: float, p: float
) -> float:
    """Time-constant ``|B|`` far from a helical three-phase line in the
    regime ``a << p << r`` (Part 5 eq 28). Exponential decay:

        |B| ~ (3 pi mu_0 a I_peak / (4 sqrt(p r))) * exp(-2 pi r / p).
    """
    if p <= 0.0:
        raise ValueError(f"pitch p must be > 0, got {p}")
    return (
        3.0 * math.pi * MU_0 * a * I_peak
        / (4.0 * math.sqrt(p * r))
        * math.exp(-2.0 * math.pi * r / p)
    )


# ---------------------------------------------------------------------------
# Builders for arrangement positions and balanced phasor currents
# ---------------------------------------------------------------------------


def triangle_positions(a: float) -> np.ndarray:
    """Cartesian positions of three conductors at radius ``a`` and
    angles 0, 2 pi / 3, 4 pi / 3.
    """
    angles = np.array([0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0])
    return np.stack([a * np.cos(angles), a * np.sin(angles)], axis=-1)


def planar_positions(d: float) -> np.ndarray:
    """Three coplanar conductors at ``(-d, 0)``, ``(0, 0)``, ``(d, 0)``."""
    return np.array([[-d, 0.0], [0.0, 0.0], [d, 0.0]])


def hexagon_positions(a: float) -> np.ndarray:
    """Six conductors at the vertices of a regular hexagon of radius ``a``,
    listed in order ``angle = 0, pi/3, 2 pi/3, pi, 4 pi/3, 5 pi/3``.

    Note
    ----
    The 1/r**3 hexagon asymptote in PDF eq 44 requires a specific
    phase-on-vertex assignment that is not encoded by
    :func:`balanced_six_phase_currents` below. Combining these two
    helpers gives a 1/r**2 (dipole-dominant) far field. See the
    module docstring for context.
    """
    angles = np.arange(6) * (math.pi / 3.0)
    return np.stack([a * np.cos(angles), a * np.sin(angles)], axis=-1)


def balanced_three_phase_currents(I_peak: float) -> np.ndarray:
    """Three balanced phasor currents of peak amplitude ``I_peak``.

    For an RMS current ``I_rms`` pass ``I_peak = sqrt(2) * I_rms``.
    """
    q = np.exp(-1j * 2.0 * math.pi / 3.0)
    return np.array([I_peak, I_peak * q, I_peak * q ** 2])


def balanced_six_phase_currents(I_peak: float) -> np.ndarray:
    """Six phasor currents corresponding to a balanced three-phase set
    plus its anti-phase replica:

        [I_peak, I_peak q, I_peak q**2, -I_peak, -I_peak q, -I_peak q**2]

    with ``q = exp(-j 2 pi / 3)``. Useful as a building block; see
    the note in :func:`hexagon_positions` about quadrupole arrangements.
    """
    three = balanced_three_phase_currents(I_peak)
    return np.concatenate([three, -three])
