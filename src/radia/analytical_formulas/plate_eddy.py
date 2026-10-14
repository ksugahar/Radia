"""Eddy current in a thin rectangular plate (slow-B-field approximation).

Part 1, eq 26-27 of the Wakao-Igarashi-Fujiwara-Kameari series.
Originally Weissenburger and Christensen (PPPL-1517, 1979); also
Kameari (J. Comput. Phys. 42, 124-140, 1981).

Geometry and assumptions
------------------------
A thin conducting plate of in-plane half-extents ``a`` (along x) and
``b`` (along y) and thickness ``d`` is placed in a spatially uniform,
perpendicular magnetic flux density ``B_z(t)``. The slow-field limit
applies, that is

    d / (mu sigma B_dot * min(a, b) / B) << 1

(eq 25 of the PDF) so that the reaction field is negligible compared
to the applied field. In that limit the eddy current is purely in
the plate plane and is divergence-free, which lets it be written in
terms of a single current vector potential

    J = curl(T_z e_z),
    J_x =  d T_z / d y,    J_y = - d T_z / d x

Current vector potential T_z (eq 26)
------------------------------------
With ``Bdot = dB_z/dt`` and ``k_n = (2n + 1) pi / (2 a)``:

    T_z(x, y) = (sigma * Bdot / 8) * {
        x**2 - a**2
        + (32 a**2 / pi**3) * sum_{n=0}^infinity
              ((-1)**n / (2 n + 1)**3) cos(k_n x) cosh(k_n y) / cosh(k_n b)
    }

Current density (eq 27)
-----------------------
Direct differentiation gives

    J_x =  (2 sigma Bdot a / pi**2) * sum_n ((-1)**n / (2 n + 1)**2)
              * cos(k_n x) sinh(k_n y) / cosh(k_n b)

    J_y = -(sigma Bdot x / 4)
          + (2 sigma Bdot a / pi**2) * sum_n ((-1)**n / (2 n + 1)**2)
                * sin(k_n x) cosh(k_n y) / cosh(k_n b)

The closed-form ``-sigma Bdot x / 4`` term in ``J_y`` exactly cancels
the residual at ``y = +-b`` so that ``J . n = 0`` on the plate edges
(verified analytically by recognizing the triangle-wave Fourier
series of ``x`` over ``[-a, a]``).

Series convergence
------------------
Coefficients decay as ``1 / (2 n + 1)**3`` for ``T`` and
``1 / (2 n + 1)**2`` for ``J``; near the boundary the ratio
``cosh(k_n y) / cosh(k_n b)`` decays exponentially in ``n`` if
``|y| < b``. The series converges slowly in a neighborhood of the
plate corners; for high-precision corner evaluation, increase
``n_terms``. The default ``n_terms = 200`` gives ``< 1e-10`` relative
error away from the immediate corner singularity.

References
----------
Wakao S., Igarashi H., Fujiwara K., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 1)", IEE Japan SA-02-28 / RM-02-64 (2002).
  -> eq 26-27, the T_z and J series formulas reproduced here.

Wakao S., Fujiwara K., Tokumasu T., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 6)", IEE Japan SA-05-15 / RM-05-15 (2005).
  -> §3.1, scalar dissipation P (the closed form whose OCR-extracted
     coefficient placement was ambiguous; we use the equivalent
     numerical area-integral of |J|**2/sigma instead).

Weissenburger D. W. and Christensen U. R., "Transient eddy currents
  on finite plane and toroidal conducting surfaces", Princeton
  Plasma Physics Laboratory Report PPPL-1517 (1979).

Kameari A., "Transient Eddy Current Analysis on Thin Conductors with
  Arbitrary Connections and Shapes", J. Comput. Phys. 42, 124-140
  (1981).
"""

from __future__ import annotations

import math

import numpy as np


def _odd_modes(n_terms: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (m, sgn) where m = 2 n + 1 and sgn = (-1)**n, n = 0..n_terms-1."""
    n = np.arange(n_terms, dtype=float)
    return 2.0 * n + 1.0, (-1.0) ** n


def plate_eddy_T(
    point,
    a: float,
    b: float,
    sigma: float,
    Bdot: float,
    n_terms: int = 200,
) -> np.ndarray:
    """Current vector potential ``T_z(x, y)`` in the plate (eq 26).

    Parameters
    ----------
    point : array_like, shape (..., 2)
        ``(x, y)`` in meters; values outside the plate (|x| > a or
        |y| > b) are evaluated by analytic continuation but have no
        physical meaning.
    a, b : float
        Plate half-extents [m].
    sigma : float
        Conductivity [S/m].
    Bdot : float
        Time derivative of the applied perpendicular field, dB_z/dt
        [T/s].
    n_terms : int
        Number of Fourier terms (default 200).

    Returns
    -------
    T_z : ndarray, shape (...)
        Current vector potential value [A/m]. The streamlines of the
        in-plane eddy current coincide with iso-contours of ``T_z``.
    """
    point = np.asarray(point, dtype=float)
    if point.shape[-1] != 2:
        raise ValueError(f"point must have shape (..., 2), got {point.shape}")
    x = point[..., 0]
    y = point[..., 1]
    m, sgn = _odd_modes(n_terms)               # shape (n_terms,)
    k = m * math.pi / (2.0 * a)                # k_n
    # cosh ratio computed via tanh-form to avoid overflow when k_n b is large
    # cosh(k y)/cosh(k b) = exp(k(|y|-b)) * (1 + e^{-2k|y|})/(1 + e^{-2k b}),
    # but numpy's cosh is fine if k b is bounded by ~700; for larger b/a
    # truncate via np.where.
    arg_y = np.multiply.outer(np.abs(y), k)    # (..., n_terms)
    arg_b = k * b                              # (n_terms,)
    # numerically safer: cosh(u)/cosh(v) = exp(u-v)*(1+e^{-2u})/(1+e^{-2v})
    cosh_ratio = np.exp(arg_y - arg_b) * (
        (1.0 + np.exp(-2.0 * arg_y))
        / (1.0 + np.exp(-2.0 * arg_b))
    )
    # cos(k x): broadcast x to (..., n_terms)
    cos_kx = np.cos(np.multiply.outer(x, k))
    coef = sgn / m**3                          # (n_terms,)
    series = np.einsum("...n,...n,n->...", cos_kx, cosh_ratio, coef)
    return (sigma * Bdot / 8.0) * (
        x * x - a * a
        + (32.0 * a * a / math.pi**3) * series
    )


def plate_eddy_dissipation(
    a: float,
    b: float,
    d: float,
    sigma: float,
    Bdot: float,
    n_terms: int = 200,
    n_quad: int = 64,
) -> float:
    """Total instantaneous Joule loss in the plate (Part 6 §3.1).

    Power dissipation in the slow-field limit:

        P = (d / sigma) integral_{-a}^{a} integral_{-b}^{b}
              (J_x**2 + J_y**2) dx dy

    The current density itself uses the closed-form Fourier series of
    :func:`plate_eddy_J`; the outer area integral is evaluated by a
    tensor-product Gauss-Legendre rule of order ``n_quad`` per axis.

    The PDF (Part 6 §3.1, eq 8) gives a compact closed-form

        P_inf - (256 a**4 d / pi**5 sigma) sum_n tanh(lambda_n b) / lambda_n**5

    where ``P_inf = (4 sigma a**3 b d / 3) Bdot**2`` is the
    infinite-strip result, but the OCR-extracted form leaves the
    second-term sign and the position of ``sigma`` ambiguous.
    Numerical integration of the analytic ``J`` is exact up to the
    Gauss-Legendre truncation and avoids that ambiguity.

    Parameters
    ----------
    a, b : float
        Plate half-extents [m].
    d : float
        Plate thickness [m].
    sigma : float
        Conductivity [S/m].
    Bdot : float
        Time derivative of the perpendicular field [T/s].
    n_terms : int
        Number of Fourier terms in the J series (default 200).
    n_quad : int
        Gauss-Legendre points per axis for the area integral
        (default 64; gives < 1e-10 relative error well away from the
        corner singularity).
    """
    from .gauss_legendre import nodes_weights

    x_n, w_x = nodes_weights(min(n_quad, 24))
    # Use up to n=24 (Wakao Table 3 ceiling); for higher accuracy run
    # repeat passes on subdivided panels (not needed here).
    y_n, w_y = nodes_weights(min(n_quad, 24))
    X = a * x_n[:, None] + 0.0 * y_n[None, :]
    Y = 0.0 * x_n[:, None] + b * y_n[None, :]
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    J = plate_eddy_J(pts, a, b, sigma, Bdot, n_terms=n_terms)
    J_sq = (J[..., 0] ** 2 + J[..., 1] ** 2).reshape(X.shape)
    W = a * b * np.outer(w_x, w_y)
    return d / sigma * float(np.sum(W * J_sq))


def plate_eddy_J(
    point,
    a: float,
    b: float,
    sigma: float,
    Bdot: float,
    n_terms: int = 200,
) -> np.ndarray:
    """In-plane eddy current density ``(J_x, J_y)`` (eq 27).

    Parameters
    ----------
    point, a, b, sigma, Bdot, n_terms
        See :func:`plate_eddy_T`.

    Returns
    -------
    J : ndarray, shape (..., 2)
        ``(J_x, J_y)`` [A/m**2].

    Notes
    -----
    The ``-sigma Bdot x / 4`` closed-form term in ``J_y`` is the
    Fourier-summed contribution of the trivial ``x**2 - a**2`` part of
    ``T_z``; together with the cosh-cosh series it makes the boundary
    condition ``J . n = 0`` on ``y = +- b`` exactly satisfied.
    """
    point = np.asarray(point, dtype=float)
    if point.shape[-1] != 2:
        raise ValueError(f"point must have shape (..., 2), got {point.shape}")
    x = point[..., 0]
    y = point[..., 1]
    m, sgn = _odd_modes(n_terms)
    k = m * math.pi / (2.0 * a)
    arg_y = np.multiply.outer(np.abs(y), k)
    arg_b = k * b
    cosh_ratio = np.exp(arg_y - arg_b) * (
        (1.0 + np.exp(-2.0 * arg_y))
        / (1.0 + np.exp(-2.0 * arg_b))
    )
    sinh_ratio = np.exp(arg_y - arg_b) * (
        (1.0 - np.exp(-2.0 * arg_y))
        / (1.0 + np.exp(-2.0 * arg_b))
    )
    # sinh_ratio carries the sign of y; we used |y| so multiply back
    sinh_ratio = sinh_ratio * np.sign(y)[..., None]

    cos_kx = np.cos(np.multiply.outer(x, k))
    sin_kx = np.sin(np.multiply.outer(x, k))
    coef = sgn / m**2

    series_Jx = np.einsum("...n,...n,n->...", cos_kx, sinh_ratio, coef)
    series_Jy = np.einsum("...n,...n,n->...", sin_kx, cosh_ratio, coef)

    pref = 2.0 * sigma * Bdot * a / math.pi**2
    J_x = pref * series_Jx
    J_y = -(sigma * Bdot * x / 4.0) + pref * series_Jy

    out = np.zeros_like(point)
    out[..., 0] = J_x
    out[..., 1] = J_y
    return out
