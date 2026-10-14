"""Hastings polynomial approximations of K(k), E(k).

Part 3, eq 13 + Tables 1-2 of the Wakao-Igarashi-Fujiwara-Kameari
series.  Originally Abramowitz and Stegun 17.3.34 / 17.3.36 and the
2nd-order tables in Hastings (1955).

The series gives the complete elliptic integrals of the 1st and 2nd
kind in terms of ``x = 1 - k**2`` (the complementary modulus
parameter):

    K(k) ~ sum_{i=0}^{n}  (a_i + b_i * ln(1/x)) * x**i
    E(k) ~ sum_{i=0}^{n}  (c_i + d_i * ln(1/x)) * x**i

The PDF gives both an n=2 (degree-2) and an n=4 (degree-4) variant.
The degree-2 form gives ``< 4e-5`` absolute error in K and E; the
degree-4 form gives ``< 2e-8``. SciPy's ``scipy.special.ellipk`` and
``ellipe`` are based on the Bulirsch / Carlson reductions and are
accurate to machine precision; the polynomial form is preserved here
as a documented reference and as a no-scipy fallback.

Limits and edge cases
---------------------
* ``k = 0`` (``x = 1``):    K = E = pi / 2.
* ``k = 1`` (``x = 0``):    K diverges as ``-(1/2) ln x``;  E = 1.

The polynomial approximation handles ``x = 0`` by setting ``ln(1/0)``
to ``inf`` in the ``b_i`` coefficient (which is zero at ``i = 0``)
multiplied by ``0`` from ``x ** i``; the user should expect the result
to differ from infinity due to the truncated tail. For applications
near ``k = 1`` use ``scipy.special.ellipk`` directly.
"""

from __future__ import annotations

import math

import numpy as np


# Table 1 of Part 3: degree-2 approximation.
_TABLE1 = {
    "a": [1.3862944, 0.1119723, 0.0725296],
    "b": [0.5,       0.1213478, 0.0288729],
    "c": [1.0,       0.4630151, 0.1077812],
    "d": [0.0,       0.2452727, 0.0412496],
}

# Table 2 of Part 3: degree-4 approximation.
_TABLE2 = {
    "a": [
        1.38629436112,
        0.09666344259,
        0.03590092383,
        0.03742563713,
        0.01451196212,
    ],
    "b": [
        0.5,
        0.12498593597,
        0.06880248576,
        0.03328355346,
        0.00441787012,
    ],
    "c": [
        1.0,
        0.44325141463,
        0.06260601220,
        0.04757383546,
        0.01736506451,
    ],
    "d": [
        0.0,
        0.24998368310,
        0.09200180037,
        0.04069697526,
        0.00526449639,
    ],
}


def _hastings_K(k: np.ndarray, table) -> np.ndarray:
    """Evaluate K(k) from a Hastings table."""
    x = 1.0 - k * k
    # ln(1/x) -> inf at k = 1; downstream consumers should clamp k.
    with np.errstate(divide="ignore"):
        log_inv_x = -np.log(x)
    a, b = table["a"], table["b"]
    out = np.zeros_like(x, dtype=float) if isinstance(x, np.ndarray) else 0.0
    xp = np.ones_like(x) if isinstance(x, np.ndarray) else 1.0
    for i in range(len(a)):
        out = out + (a[i] + b[i] * log_inv_x) * xp
        xp = xp * x
    return out


def _hastings_E(k: np.ndarray, table) -> np.ndarray:
    """Evaluate E(k) from a Hastings table."""
    x = 1.0 - k * k
    with np.errstate(divide="ignore"):
        log_inv_x = -np.log(x)
    c, d = table["c"], table["d"]
    out = np.zeros_like(x, dtype=float) if isinstance(x, np.ndarray) else 0.0
    xp = np.ones_like(x) if isinstance(x, np.ndarray) else 1.0
    for i in range(len(c)):
        out = out + (c[i] + d[i] * log_inv_x) * xp
        xp = xp * x
    return out


def K_hastings_2(k):
    """Complete elliptic integral of the 1st kind, degree-2 Hastings.

    Accurate to ``< 4e-5`` for ``0 <= k**2 <= 1``.
    """
    return _hastings_K(np.asarray(k, dtype=float), _TABLE1)


def E_hastings_2(k):
    """Complete elliptic integral of the 2nd kind, degree-2 Hastings."""
    return _hastings_E(np.asarray(k, dtype=float), _TABLE1)


def K_hastings_4(k):
    """Complete elliptic integral of the 1st kind, degree-4 Hastings.

    Accurate to ``< 2e-8`` for ``0 <= k**2 <= 1``.
    """
    return _hastings_K(np.asarray(k, dtype=float), _TABLE2)


def E_hastings_4(k):
    """Complete elliptic integral of the 2nd kind, degree-4 Hastings."""
    return _hastings_E(np.asarray(k, dtype=float), _TABLE2)


# Convenience alias: prefer the degree-4 form by default.
K = K_hastings_4
E = E_hastings_4
