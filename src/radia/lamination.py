"""Material homogenization helpers for laminated magnetic steel."""

from __future__ import annotations

import cmath
import math


MU0 = 4.0e-7 * math.pi


def laminated_mu_eff(
    mu_r: float,
    sigma: float,
    omega: float,
    d_lam: float,
    fill: float = 1.0,
) -> complex:
    """Return the in-plane complex permeability of a laminated stack.

    ``d_lam`` is the conducting sheet thickness and ``fill`` is the steel
    volume fraction. The convention is ``exp(+j omega t)``; eddy-current
    shielding therefore gives a negative imaginary permeability.
    """

    if omega == 0 or sigma == 0:
        return complex(MU0 * (fill * mu_r + (1.0 - fill)))
    b = (d_lam / 2.0) * cmath.sqrt(1j * omega * MU0 * mu_r * sigma)
    factor = cmath.tanh(b) / b
    return MU0 * (fill * mu_r * factor + (1.0 - fill))
