"""Gauss-Legendre quadrature nodes and weights for n = 1..24.

Part 3, eq 14-21 + Table 3 of the Wakao-Igarashi-Fujiwara-Kameari
series.

The ``n``-point Gauss-Legendre rule on the canonical interval
``[-1, +1]`` uses the n roots ``x_i`` of the Legendre polynomial
``P_n(x)`` as quadrature points, with weights

    w_i = 2 (1 - x_i**2) / (n * P_{n-1}(x_i))**2

A polynomial of degree ``2 n - 1`` is integrated exactly. To use the
rule on a general interval ``[a, b]`` apply the affine map
``t = (b - a)/2 * x + (b + a)/2`` and rescale weights by ``(b - a)/2``;
helper :func:`integrate` does both.

This module pre-computes the table that the original Wakao paper
prints (Tables 3-a and 3-b) at runtime via NumPy's
``numpy.polynomial.legendre.leggauss``, which returns the same nodes
and weights to machine precision.  Holding the values at runtime
rather than literal-encoding the printed table preserves the full
machine-precision accuracy ``leggauss`` provides while keeping the
``n_max = 24`` ceiling that the PDF documents.
"""

from __future__ import annotations

import numpy as np


N_MAX = 24


def nodes_weights(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(x_i, w_i)`` for n-point Gauss-Legendre on [-1, +1].

    Parameters
    ----------
    n : int
        Number of quadrature points, ``1 <= n <= 24``.

    Returns
    -------
    nodes, weights : ndarray of shape (n,)
        Nodes are sorted ascending.
    """
    if not (1 <= n <= N_MAX):
        raise ValueError(
            f"n must satisfy 1 <= n <= {N_MAX} (the PDF Table 3 ceiling), "
            f"got {n}"
        )
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def integrate(f, a: float, b: float, n: int = 16) -> float:
    """Approximate ``int_a^b f(x) dx`` by an n-point Gauss-Legendre rule.

    Parameters
    ----------
    f : callable
        Integrand; called with a 1D ndarray of evaluation points and
        expected to return a 1D ndarray of values.
    a, b : float
        Integration limits.
    n : int
        Number of quadrature points. The default ``n = 16`` is exact
        for polynomials of degree up to 31.

    Returns
    -------
    integral : float
    """
    x, w = nodes_weights(n)
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    pts = half * x + mid
    return half * float(np.sum(w * f(pts)))


def integrate_2d(f, a: float, b: float, c: float, d: float,
                 nx: int = 16, ny: int = 16) -> float:
    """Approximate ``int_c^d int_a^b f(x, y) dx dy`` by a tensor-product
    Gauss-Legendre rule.

    Parameters
    ----------
    f : callable
        Integrand; called with two 2D ndarrays ``(X, Y)`` of the same
        shape and expected to return an ndarray of the same shape.
    a, b, c, d : float
        Integration limits.
    nx, ny : int
        Number of points in each direction.
    """
    x_, wx = nodes_weights(nx)
    y_, wy = nodes_weights(ny)
    halfx, midx = 0.5 * (b - a), 0.5 * (a + b)
    halfy, midy = 0.5 * (d - c), 0.5 * (c + d)
    px = halfx * x_ + midx
    py = halfy * y_ + midy
    X, Y = np.meshgrid(px, py, indexing="xy")
    W = np.outer(wy, wx)
    return halfx * halfy * float(np.sum(W * f(X, Y)))
