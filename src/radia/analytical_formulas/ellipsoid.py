"""Rotational ellipsoid: demagnetization factor, internal field, torque.

Part 5, eq 38-44 of the Wakao-Igarashi-Fujiwara-Kameari series.

Geometry
--------
A rotational ellipsoid has semi-axes (a, a, c) with c the polar (symmetry)
axis. Following the PDF, c is denoted ``b`` (the polar semi-axis); ``a`` is
the equatorial semi-axis.

    prolate (cigar):  c > a   --  axis longer than equator
    oblate  (disk):   c < a   --  axis shorter than equator
    sphere:           c = a

Demagnetization factor
----------------------
The demagnetization tensor of an ellipsoid in a uniform external field is
diagonal in the principal-axis frame. By azimuthal symmetry the two
equatorial factors equal ``(1 - N_z) / 2`` and only the polar factor
``N_z`` is independent. The sum of the three factors equals 1.

Prolate, ``e = sqrt(1 - (a/c)**2)`` (eq 40):

    N_z = (1 - e**2) / e**3 * (atanh(e) - e)

Oblate, ``e = sqrt((a/c)**2 - 1)`` (eq 39):

    N_z = (1 + e**2) / e**3 * (e - atan(e))

Both branches share the limit ``N_z -> 1/3`` as ``c/a -> 1`` (sphere); for
aspect ratios near unity the implementation switches to a Taylor series
to avoid 0/0.

Internal field (eq 38)
----------------------
For a linear magnetic material with relative susceptibility ``chi_r``
(``= mu_r - 1``) placed in a uniform external field ``H_0`` aligned with
a principal axis ``i``,

    H_i = H_0 / (1 + chi_r * N_i)

Torque (eq 44)
--------------
With the symmetry axis tilted by angle ``alpha`` from the uniform field,

    T_z = - mu_0 * chi_r**2 * (1 - 3 N_z)
          / (2 (1 + chi_r N_z) (2 + chi_r (1 - N_z)))
          * sin(2 alpha) * H_0**2 * V

with ``V = (4/3) pi a**2 c``. ``T_z < 0`` (restoring) for prolate
(``N_z < 1/3``); the axis is pulled toward the field.

References
----------
Part 5 (SA-04-44 / RM-04-68, 2004), eq 38-44.
Osborn, J.A., Phys. Rev. 67 (1945) 351-357 -- general (a, b, c) form.
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4.0e-7 * math.pi


def _demag_series_near_sphere(eps: float) -> float:
    """Polar demag factor near sphere, ``eps = (c/a)**2 - 1``.

    Series: ``N_z = 1/3 - 2 eps / 15 + 8 eps**2 / 105 - ...``

    Valid for both prolate (``eps > 0``) and oblate (``eps < 0``);
    derived by expanding ``e**2 = eps / (1 + eps)`` and substituting in
    ``N_z = 1/3 - 2 e**2 / 15 - 2 e**4 / 35 - ...`` (the ``e``-series
    obtained from atanh / atan around 0).
    """
    return 1.0 / 3.0 - 2.0 * eps / 15.0 + 8.0 * eps * eps / 105.0


def demag_factor_prolate(c: float, a: float) -> float:
    """Polar demag factor of a prolate rotational ellipsoid.

    Parameters
    ----------
    c : float
        Polar semi-axis (``c > a``).
    a : float
        Equatorial semi-axis.

    Returns
    -------
    N_z : float
        Demagnetization factor along the polar axis. Equatorial factors
        are ``(1 - N_z) / 2``.
    """
    if c <= 0.0 or a <= 0.0:
        raise ValueError("c, a must be positive")
    if c < a:
        raise ValueError(
            f"prolate requires c >= a, got c={c}, a={a}; "
            f"call demag_factor_oblate instead"
        )
    eps = (c / a) ** 2 - 1.0
    if eps < 1.0e-4:
        return _demag_series_near_sphere(eps)
    e = math.sqrt(1.0 - (a / c) ** 2)
    return (1.0 - e * e) / (e ** 3) * (math.atanh(e) - e)


def demag_factor_oblate(c: float, a: float) -> float:
    """Polar demag factor of an oblate rotational ellipsoid.

    Parameters
    ----------
    c : float
        Polar semi-axis (``c < a``).
    a : float
        Equatorial semi-axis.

    Returns
    -------
    N_z : float
        Demagnetization factor along the polar axis.
    """
    if c <= 0.0 or a <= 0.0:
        raise ValueError("c, a must be positive")
    if c > a:
        raise ValueError(
            f"oblate requires c <= a, got c={c}, a={a}; "
            f"call demag_factor_prolate instead"
        )
    eps = (c / a) ** 2 - 1.0  # negative
    if abs(eps) < 1.0e-4:
        return _demag_series_near_sphere(eps)
    e = math.sqrt((a / c) ** 2 - 1.0)
    return (1.0 + e * e) / (e ** 3) * (e - math.atan(e))


def demag_factor_rotational(c: float, a: float) -> tuple[float, float, float]:
    """Demagnetization factors ``(N_x, N_y, N_z)`` of a rotational ellipsoid.

    ``c`` is the polar (z) semi-axis, ``a`` the equatorial (x = y).
    The sphere case is selected when ``|c/a - 1| < 1e-4``.

    Returns
    -------
    (N_x, N_y, N_z) : tuple of float
        With ``N_x = N_y = (1 - N_z) / 2`` and ``N_x + N_y + N_z = 1``.
    """
    if c >= a:
        Nz = demag_factor_prolate(c, a)
    else:
        Nz = demag_factor_oblate(c, a)
    Nxy = 0.5 * (1.0 - Nz)
    return (Nxy, Nxy, Nz)


def ellipsoid_internal_field(
    H0: float,
    chi_r: float,
    N: float,
) -> float:
    """Internal H along a principal axis (eq 38).

    ``H_i = H_0 / (1 + chi_r * N_i)``

    Parameters
    ----------
    H0 : float
        External field along the axis [A/m].
    chi_r : float
        Relative susceptibility, ``chi_r = mu_r - 1``.
    N : float
        Demag factor along that axis.
    """
    return H0 / (1.0 + chi_r * N)


def ellipsoid_torque(
    H0: float,
    alpha: float,
    chi_r: float,
    a: float,
    c: float,
) -> float:
    """Torque on a tilted rotational ellipsoid in a uniform field (eq 44).

    The symmetry axis is tilted by ``alpha`` (rad) from the uniform
    external field of magnitude ``H0``. ``T_z`` is the torque component
    perpendicular to both the symmetry axis and the field; the sign is
    such that ``T_z < 0`` is a restoring torque pulling ``alpha -> 0``
    (axis aligns with field), which is the prolate / soft-iron case.

    Parameters
    ----------
    H0 : float
        External field magnitude [A/m].
    alpha : float
        Tilt angle [rad] between symmetry axis and external field.
    chi_r : float
        Relative susceptibility, ``chi_r = mu_r - 1``.
    a : float
        Equatorial semi-axis [m].
    c : float
        Polar semi-axis [m].

    Returns
    -------
    T_z : float
        Torque [N m].
    """
    Nz = (
        demag_factor_prolate(c, a) if c >= a else demag_factor_oblate(c, a)
    )
    V = (4.0 / 3.0) * math.pi * a * a * c
    denom = 2.0 * (1.0 + chi_r * Nz) * (2.0 + chi_r * (1.0 - Nz))
    return (
        -MU_0 * chi_r * chi_r * (1.0 - 3.0 * Nz) / denom
        * math.sin(2.0 * alpha) * H0 * H0 * V
    )
