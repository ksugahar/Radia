"""AC impedance of conductors via skin-effect 1D / Bessel solutions.

Part 6 §4 (planar skin depth) and §5 (cylindrical conductor with
Bessel solution) of the Wakao-Igarashi-Fujiwara-Kameari series.

Planar surface impedance (eq 17-22)
-----------------------------------
A semi-infinite conductor of conductivity ``sigma`` in a tangential
sinusoidal magnetic field of angular frequency ``omega`` has the
1D diffusion solution

    E_x(z) = E_x0 exp(-k z),   k = (1 - j) / d_skin,
    d_skin = sqrt(2 / (omega mu_0 sigma)).

The surface impedance ratio E_t / H_t at z = 0 is

    Z_s = (1 + j) / (sigma d_skin)
        = (1 + j) sqrt(omega mu_0 / (2 sigma)).

``|Z_s| = sqrt(omega mu_0 / sigma)``;  ``arg(Z_s) = pi / 4``.

Cylindrical conductor (eq 28-34)
--------------------------------
An infinite z-directed current ``I`` in a solid cylinder of radius
``a`` and conductivity ``sigma``: the 1D diffusion equation in
cylindrical coordinates has Bessel solutions

    E_z(r) = (I k / (2 pi a sigma)) J_0(k r) / J_1(k a),
    k**2 = - j omega mu_0 sigma  (so k = sqrt(-j omega mu_0 sigma)),

giving impedance per unit length

    Z = E_z(a) / I = k J_0(ka) / (2 pi a sigma J_1(ka)).

Limits:

* ``ka << 1`` (low frequency / thick skin): ``Z -> R_dc + j omega L_int``
  with DC resistance ``R_dc = 1 / (pi a**2 sigma)`` and internal
  inductance ``L_int = mu_0 / (8 pi)`` per unit length.
* ``ka >> 1`` (high frequency / thin skin): ``Z -> Z_s / (2 pi a)``
  with the planar surface impedance.

References
----------
Wakao S., Fujiwara K., Tokumasu T., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 6)", IEE Japan SA-05-15 / RM-05-15 (2005).
  -> §4 planar surface impedance (eq 17-22), §5 cylindrical conductor
     AC impedance (eq 28-34).

Stratton J. A., "Electromagnetic Theory", McGraw-Hill (1941),
  §5.16-5.19 (skin effect in straight wires; the full Bessel solution
  for solid cylinders).

Wheeler H. A., "Formulas for the Skin Effect", Proc. IRE 30 (1942),
  pp. 412-424.

Schelkunoff S. A., "The Electromagnetic Theory of Coaxial Transmission
  Lines and Cylindrical Shields", Bell Syst. Tech. J. 13 (1934),
  pp. 532-579.

Montgomery D. B., "Solenoid Magnet Design", Wiley (1969).

Dowell P. L., "Effects of eddy currents in transformer windings",
  Proc. IEE 113(8), pp. 1387-1394 (1966).
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

MU_0 = 4.0e-7 * math.pi


def skin_depth(sigma: float, omega: float, mu_r: float = 1.0) -> float:
    """Skin depth ``d = sqrt(2 / (omega mu_0 mu_r sigma))`` in meters."""
    if sigma <= 0 or omega < 0 or mu_r < 1:
        raise ValueError(
            f"sigma > 0, omega >= 0, mu_r >= 1 required;"
            f" got sigma={sigma}, omega={omega}, mu_r={mu_r}"
        )
    if omega == 0:
        return float("inf")
    return math.sqrt(2.0 / (omega * MU_0 * mu_r * sigma))


def planar_surface_impedance(
    sigma: float, omega: float, mu_r: float = 1.0
) -> complex:
    """Planar surface impedance ``Z_s = (1 + j) / (sigma d_skin)``
    (Part 6 §4, eq 22).

    Returns
    -------
    Z_s : complex
        Surface impedance (Ohm). ``|Z_s| = sqrt(omega mu_0 mu_r / sigma)``;
        ``arg(Z_s) = pi / 4``.
    """
    d = skin_depth(sigma, omega, mu_r)
    if math.isinf(d):
        return 0.0 + 0.0j
    return (1.0 + 1.0j) / (sigma * d)


def dowell_rectangular_ac_impedance(
    width: float,
    thickness: float,
    sigma: float,
    omega: float,
    mu_r: float = 1.0,
) -> complex:
    """Single-layer Dowell impedance of a rectangular conductor [Ohm/m].

    The conductor carries current along its length.  ``width`` and
    ``thickness`` are normalized internally so the smaller dimension is
    the diffusion thickness.  The returned impedance is

    ``Z = R_dc q coth(q)``, ``q = (1+j) thickness / (2 delta)``.

    This is the one-dimensional broad-face solution for a rectangular
    strip.  It retains the DC limit and internal reactance and approaches
    two independent planar surfaces at high frequency.  It does not invent
    a corner radius or an equivalent circular section; edge crowding and
    proximity redistribution remain the responsibility of the PEEC
    perimeter solve.
    """
    if width <= 0 or thickness <= 0 or sigma <= 0:
        raise ValueError(
            "width, thickness, sigma must be > 0; got "
            f"width={width}, thickness={thickness}, sigma={sigma}"
        )
    if omega < 0 or mu_r < 1:
        raise ValueError(
            f"omega >= 0, mu_r >= 1 required; got omega={omega}, mu_r={mu_r}"
        )

    broad_width = max(float(width), float(thickness))
    diffusion_thickness = min(float(width), float(thickness))
    r_dc = 1.0 / (sigma * broad_width * diffusion_thickness)
    if omega == 0.0:
        return complex(r_dc, 0.0)

    delta = skin_depth(sigma, omega, mu_r)
    q = (1.0 + 1.0j) * diffusion_thickness / (2.0 * delta)
    abs_q = abs(q)
    if abs_q < 1.0e-3:
        # q*coth(q) = 1 + q^2/3 - q^4/45 + 2q^6/945 + O(q^8).
        q2 = q * q
        factor = 1.0 + q2 / 3.0 - q2 * q2 / 45.0 \
            + 2.0 * q2 * q2 * q2 / 945.0
    elif q.real > 30.0:
        # coth(q) is one to machine precision and direct sinh/cosh can
        # overflow in the strong-skin regime.
        factor = q
    else:
        factor = q / np.tanh(q)
    return complex(r_dc * factor)


def cylinder_ac_impedance(
    a: float,
    sigma: float,
    omega: float,
    mu_r: float = 1.0,
    bessel_func: Callable | None = None,
) -> complex:
    """AC impedance per unit length of a solid cylindrical conductor
    (Part 6 §5, eq 34).

        Z = k J_0(k a) / (2 pi a sigma J_1(k a)),    k = sqrt(-j omega mu sigma).

    Parameters
    ----------
    a : float
        Conductor radius [m].
    sigma : float
        Conductivity [S/m].
    omega : float
        Angular frequency [rad/s]; pass 0 for DC.
    mu_r : float
        Relative permeability of the conductor (default 1).
    bessel_func : callable, optional
        Override for ``scipy.special.jv``-style Bessel evaluator
        ``f(order, complex_arg) -> complex``. Defaults to
        ``scipy.special.jv``.

    Returns
    -------
    Z : complex
        Impedance per unit length [Ohm/m].
    """
    if a <= 0 or sigma <= 0:
        raise ValueError(
            f"a, sigma must be > 0; got a={a}, sigma={sigma}"
        )
    if omega == 0.0:
        # DC limit:  R = 1 / (pi a**2 sigma); zero (or only DC) inductance.
        return 1.0 / (math.pi * a * a * sigma) + 0.0j
    # k**2 = -j omega mu sigma. Pick the principal branch sqrt
    # giving Re(k) > 0 so the field decays into the conductor.
    k_sq = -1.0j * omega * MU_0 * mu_r * sigma
    k = np.sqrt(k_sq)
    ka = k * a
    if bessel_func is None:
        from scipy.special import jv as bessel_func  # noqa: WPS433
    J0 = bessel_func(0, ka)
    J1 = bessel_func(1, ka)
    return complex(ka * J0 / (2.0 * math.pi * a * a * sigma * J1))


def cylinder_dc_resistance(a: float, sigma: float) -> float:
    """DC limit of :func:`cylinder_ac_impedance`."""
    return 1.0 / (math.pi * a * a * sigma)


def cylinder_internal_inductance() -> float:
    """DC (low-frequency) internal inductance per unit length of a solid
    cylindrical conductor: ``mu_0 / (8 pi)``.
    """
    return MU_0 / (8.0 * math.pi)
