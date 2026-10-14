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


def laminated_mu_perpendicular(mu_r: float, fill: float = 1.0) -> float:
    """Return the static normal permeability of a steel/vacuum layered stack.

    Normal B is continuous across the layers. Averaging H gives
    H_avg = B * (fill / (MU0 * mu_r) + (1 - fill) / MU0), hence the
    harmonic mixture below. ``fill`` is the steel volume fraction.
    This is a linear static material law; it predicts no transverse eddy loss.
    The same module is available through MATLAB's named Python batch adapter.
    """
    mu_r = float(mu_r)
    fill = float(fill)
    if not math.isfinite(mu_r) or mu_r <= 0.0:
        raise ValueError("mu_r must be finite and positive")
    if not math.isfinite(fill) or not 0.0 <= fill <= 1.0:
        raise ValueError("fill must be finite and between zero and one")
    if fill == 0.0:
        return MU0
    if fill == 1.0:
        return MU0 * mu_r
    return MU0 / (fill / mu_r + (1.0 - fill))
