"""2D uniformly magnetized rectangular bar.

Part 2, eq 2-3 of the Wakao-Igarashi-Fujiwara-Kameari series.

Geometry
--------
A bar of cross-section ``2 a x 2 b`` centred on the origin, infinitely
long along z, with uniform 2D magnetization ``M = (M_x, M_y)`` in the
x-y plane. ``M`` is in A / m (Radia convention); the polarisation is
``J = mu_0 M`` in Tesla.

Vector potential A_z (eq 2)
---------------------------
Define the corner-relative coordinates

    u_i = x_p - x_i,   x_i in {-a, +a},   i = 1, 2
    v_j = y_p - y_j,   y_j in {-b, +b},   j = 1, 2

and the auxiliary functions

    G(u, v)      = u ln(u**2 + v**2) + 2 v atan(u / v)
    G_swap(u, v) = v ln(u**2 + v**2) + 2 u atan(v / u)

The vector potential is

    A_z = (mu_0 / 4 pi) * sum_{i,j=1,2} (-1)**(i+j)
          * [M_x * G(u_i, v_j) - M_y * G_swap(u_i, v_j)]

The ``M_x`` term is exactly Part 2 eq (2). The ``M_y`` term comes by
the same derivation applied to the surface magnetization currents on
the left / right faces (see derivation in the doc string at the bottom
of the file).

Field B (eq 3)
--------------
By ``B = curl A`` with ``A = A_z e_z``:

    B_x =  partial A_z / partial y
        =  (mu_0 / 2 pi) * sum (-1)**(i+j) * [
                M_x * atan(u_i / v_j) - M_y * (1/2) ln(u_i**2 + v_j**2)
            ]
    B_y = -partial A_z / partial x
        = -(mu_0 / 4 pi) * sum (-1)**(i+j) * [
                M_x * ln(u_i**2 + v_j**2) - 2 M_y * atan(v_j / u_i)
            ]

The constants in ``ln`` and ``atan`` derivatives drop under the
alternating ``(-1)**(i+j)`` sum.

Branch choice
-------------
The arctangent uses the principal branch ``[-pi/2, pi/2]`` (i.e.
``np.arctan``, not ``np.arctan2``). Across the line ``v = 0`` the
``atan(u/v)`` term has a jump of pi but is multiplied by ``v``, giving
a continuous ``v atan(u/v)``; numerically ``np.arctan(np.inf) == pi/2``
and the product ``2 v * pi/2 -> 0`` as ``v -> 0`` is recovered to
machine precision when ``v == 0`` is masked. The formulae are exact
everywhere outside the rectangle and on its boundary except the four
corners where ``u_i == 0`` and ``v_j == 0`` simultaneously.

Notes
-----
* M_x and M_y are independent linear contributions; the formulae below
  evaluate both in a single pass.
* Long limit: as ``b -> infinity`` the bar becomes a 1D slab and ``B``
  approaches the textbook step ``B_x = mu_0 M_x`` inside.

References
----------
Part 2 (SA-03-08 / RM-03-08, 2003), eq 2-3.
Stratton J.A., "Electromagnetic Theory", McGraw-Hill (1941).
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4.0e-7 * math.pi


def _safe_atan(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """Principal-branch arctangent of ``num / den`` with den == 0 mapped
    to ``+- pi/2``. Used so that the multiplicative factor ``den`` makes
    the product ``den * arctan(num/den)`` continuous at ``den == 0``.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = num / den
        out = np.arctan(ratio)
    # arctan(nan) is nan; that only arises at corner (num == 0 == den)
    return np.where(np.isfinite(out), out, 0.0)


def _safe_log_sq(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """``ln(u**2 + v**2)`` with the corner ``u == v == 0`` clamped to
    zero. Caller guarantees the corner only contributes a 0 * log(0)
    type term that drops out of the alternating sum.
    """
    r2 = u * u + v * v
    return np.where(r2 > 0.0, np.log(np.where(r2 > 0.0, r2, 1.0)), 0.0)


def _corner_sum(values: np.ndarray) -> np.ndarray:
    """Apply the ``(-1)**(i+j)`` alternating sign over a 2 x 2 corner
    array stacked on the LAST two axes.
    """
    # values[..., i, j] with i, j in {0, 1} corresponding to (1, 2) of PDF.
    # (-1)**(i+j): + - / - +
    sign = np.array([[+1.0, -1.0], [-1.0, +1.0]])
    return np.einsum("...ij,ij->...", values, sign)


def _gather_corner_uvs(
    point: np.ndarray,
    a: float,
    b: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(u, v)`` arrays of shape ``(..., 2, 2)`` for the four
    corners of the ``2 a x 2 b`` bar.
    """
    point = np.asarray(point, dtype=float)
    if point.shape[-1] != 2:
        raise ValueError(
            f"point must have shape (..., 2), got {point.shape}"
        )
    x = point[..., 0]
    y = point[..., 1]
    x_corners = np.array([-a, +a])
    y_corners = np.array([-b, +b])
    u = x[..., None, None] - x_corners[None, :, None]   # (..., 2, 1)
    v = y[..., None, None] - y_corners[None, None, :]   # (..., 1, 2)
    u, v = np.broadcast_arrays(u, v)
    return u, v


def rect_magnet_2d_A(
    point,
    a: float,
    b: float,
    M: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Vector potential ``A_z`` of a 2D uniformly magnetised bar.

    Parameters
    ----------
    point : array_like, shape (..., 2)
        Observation point(s) ``(x, y)`` in meters.
    a : float
        Half-extent along x.
    b : float
        Half-extent along y.
    M : (M_x, M_y)
        Magnetisation in A / m.

    Returns
    -------
    A_z : ndarray, shape (...)
        Vector potential in T m.
    """
    Mx, My = M
    u, v = _gather_corner_uvs(point, a, b)
    r2 = u * u + v * v
    log_r2 = np.where(r2 > 0.0, np.log(np.where(r2 > 0.0, r2, 1.0)), 0.0)
    G_x = u * log_r2 + 2.0 * v * _safe_atan(u, v)
    G_y = v * log_r2 + 2.0 * u * _safe_atan(v, u)
    K = MU_0 / (4.0 * math.pi)
    return K * _corner_sum(Mx * G_x - My * G_y)


def rect_magnet_2d_B(
    point,
    a: float,
    b: float,
    M: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Field ``B = (B_x, B_y)`` of a 2D uniformly magnetised bar.

    Parameters
    ----------
    point : array_like, shape (..., 2)
        Observation point(s) ``(x, y)`` in meters.
    a : float
        Half-extent along x.
    b : float
        Half-extent along y.
    M : (M_x, M_y)
        Magnetisation in A / m.

    Returns
    -------
    B : ndarray, shape (..., 2)
        ``(B_x, B_y)`` in Tesla.
    """
    Mx, My = M
    u, v = _gather_corner_uvs(point, a, b)
    r2 = u * u + v * v
    log_r2 = np.where(r2 > 0.0, np.log(np.where(r2 > 0.0, r2, 1.0)), 0.0)
    atan_uv = _safe_atan(u, v)
    atan_vu = _safe_atan(v, u)

    # B_x = +dA_z/dy
    B_x = (
        MU_0 / (2.0 * math.pi)
    ) * _corner_sum(Mx * atan_uv - 0.5 * My * log_r2)
    # B_y = -dA_z/dx
    B_y = -(
        MU_0 / (4.0 * math.pi)
    ) * _corner_sum(Mx * log_r2 - 2.0 * My * atan_vu)

    point = np.asarray(point, dtype=float)
    out = np.zeros_like(point)
    out[..., 0] = B_x
    out[..., 1] = B_y
    return out
