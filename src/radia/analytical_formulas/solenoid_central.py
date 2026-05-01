"""Central- and axial-field of a rectangular cross-section solenoid.

Part 4, eq 25-29 of the Wakao-Igarashi-Fujiwara-Kameari series.
Original reference: D.B. Montgomery, "Solenoid Magnet Design",
John Wiley & Sons, 1969 (Fabri's form factor).

Geometry
--------
Air-core solenoid with rectangular cross-section, centred on the
origin with the symmetry axis along ``z``:

    inner radius : a_1
    outer radius : a_2 (>= a_1)
    half-length  : b   (so coil extends from z = -b to z = +b)
    current density : J (A / m**2), azimuthal

Define the dimensionless aspect ratios

    alpha = a_2 / a_1,    beta = b / a_1

Central field (eq 26-27)
------------------------
The on-axis field at the origin is

    B_0 = mu_0 * a_1 * J * F(alpha, beta)

with Fabri's form factor

    F(alpha, beta) = beta * ln[ (alpha + sqrt(alpha**2 + beta**2))
                              / (1   + sqrt(1       + beta**2))  ]

Limit checks:

* infinitely long thin coil (``beta -> infinity``) :
  ``F -> alpha - 1``, recovering the textbook  B = mu_0 (a_2 - a_1) J.
* very short pancake (``beta -> 0``) :
  ``F ~ beta * ln(alpha)``, i.e. the on-axis field of a thin radial
  ring set.

Axial field (closed form)
-------------------------
Integrating Biot-Savart for a single circular hoop along ``z`` and then
over the rectangular cross-section gives

    B_z(0, z) = (mu_0 J / 2) * [
          (b - z) * ln( (a_2 + sqrt(a_2**2 + (b - z)**2))
                       /(a_1 + sqrt(a_1**2 + (b - z)**2)) )
        + (b + z) * ln( (a_2 + sqrt(a_2**2 + (b + z)**2))
                       /(a_1 + sqrt(a_1**2 + (b + z)**2)) )
    ]

At ``z = 0`` this collapses to ``mu_0 a_1 J * F(alpha, beta)``.

Notes
-----
The off-axis expansion via Legendre polynomials (PDF eq 25, 28-29 and
the Fabri-Taylor coefficient table) is not implemented here; the
analytical-magnet machinery elsewhere in Radia covers full 3D field
evaluation when needed.
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4.0e-7 * math.pi


def fabri_F(alpha: float, beta: float) -> float:
    """Fabri form factor ``F(alpha, beta)`` (eq 27).

    Parameters
    ----------
    alpha : float
        Outer-to-inner radius ratio ``a_2 / a_1`` (>= 1).
    beta : float
        Half-length over inner radius ``b / a_1`` (> 0).
    """
    if alpha < 1.0:
        raise ValueError(f"alpha = a_2/a_1 must be >= 1, got {alpha}")
    if beta <= 0.0:
        raise ValueError(f"beta = b/a_1 must be > 0, got {beta}")
    return beta * math.log(
        (alpha + math.sqrt(alpha ** 2 + beta ** 2))
        / (1.0 + math.sqrt(1.0 + beta ** 2))
    )


def central_field(a1: float, a2: float, b: float, J: float) -> float:
    """``B_z`` at the centre of a rectangular-section solenoid (eq 26).

    Parameters
    ----------
    a1, a2 : float
        Inner / outer radii [m] (``a_2 >= a_1 > 0``).
    b : float
        Half-length [m] (``b > 0``).
    J : float
        Current density [A / m**2], azimuthal.

    Returns
    -------
    B_0 : float
        On-axis ``B_z`` at the origin [Tesla].
    """
    if a1 <= 0 or a2 < a1 or b <= 0:
        raise ValueError(
            f"require 0 < a_1 <= a_2 and b > 0, got a_1={a1}, a_2={a2}, b={b}"
        )
    return MU_0 * a1 * J * fabri_F(a2 / a1, b / a1)


def axial_field(z, a1: float, a2: float, b: float, J: float) -> np.ndarray:
    """``B_z`` on the symmetry axis at distance(s) ``z`` from centre.

    Closed-form integral of Biot-Savart over the rectangular cross-
    section.

    Parameters
    ----------
    z : float or array_like
        Axial position [m]; ``z = 0`` is the geometric centre.
    a1, a2, b, J
        Geometry and current density as in :func:`central_field`.

    Returns
    -------
    B_z : ndarray (or scalar) of float
        Axial field [Tesla].
    """
    if a1 <= 0 or a2 < a1 or b <= 0:
        raise ValueError(
            f"require 0 < a_1 <= a_2 and b > 0, got a_1={a1}, a_2={a2}, b={b}"
        )
    z_arr = np.asarray(z, dtype=float)

    def _ratio_log(dz):
        return np.log(
            (a2 + np.sqrt(a2 ** 2 + dz ** 2))
            / (a1 + np.sqrt(a1 ** 2 + dz ** 2))
        )

    p = b - z_arr           # distance to top end
    q = b + z_arr           # distance to bottom end
    return 0.5 * MU_0 * J * (p * _ratio_log(p) + q * _ratio_log(q))
