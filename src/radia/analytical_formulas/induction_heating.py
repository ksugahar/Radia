"""Canonical induction-heating Joule-loss formulas (textbook references).

Phase γ deliverable (v4.23.0). After Phase α/β resolved the cuboid-
average-B closed form by in-house sympy derivation, an attempt to
absorb Stafl 1967 Chapter 4 §4.1-§4.4 closed-forms via PDF OCR was
made; sanity-test against canonical limits showed the OCR'd Stafl
formulas were systematically wrong by orders of magnitude (sphere
~1e10×, cylinder ~800× at low frequency). This module replaces that
work with **closed forms derived from first principles** and validated
against canonical textbook references.

Cylinder of conductivity ``sigma`` and permeability ``mu`` placed
inside a long axial-current solenoid producing surface current
``H_0 = N I`` (turns-per-unit-length × peak current) on the cylinder
surface. We solve the 1-D Helmholtz equation for ``H_z(r)`` inside the
cylinder, evaluate the Poynting inflow at ``r = a``, and obtain the
time-averaged Joule loss per unit length.

Derivation (recorded for reproducibility)
------------------------------------------

In sinusoidal steady state ``e^{j omega t}``, inside the cylinder
``∇²H_z - jωσμH_z = 0`` (axisymmetric, axial). Define the Stafl wave
number ``p = sqrt(omega sigma mu)``. The regular solution is

    H_z(r) = H_0 J_0(j^{3/2} p r) / J_0(j^{3/2} p a),

with Kelvin functions ``J_0(j^{3/2} x) = ber(x) + j bei(x)``. From
``∇×H = J = σ E`` (quasi-static) with ``H = H_z ẑ``:

    E_phi(r) = -(1/sigma) dH_z/dr.

Time-averaged Poynting power inflow at ``r = a`` per unit length:

    P = pi a Re[ -E_phi(a) H_z*(a) ]
      = pi a H_0^2 / sigma · Re[ -p ber'(pa)+ j bei'(pa)
                                  / (ber(pa) + j bei(pa)) ]
      = pi a H_0^2 p / sigma · (ber·ber' + bei·bei') / (ber^2 + bei^2),

i.e.

    P = (π / σ) H_0² (k a)
        · (ber(ka) ber'(ka) + bei(ka) bei'(ka))
        / (ber(ka)² + bei(ka)²),  k = sqrt(omega sigma mu).

(Note: this differs from Stafl 1967 eq 4-157 which has
``ber·bei' - bei·ber'`` in the numerator; that combination corresponds
to the *imaginary* part of the Poynting inflow, i.e. the reactive
power, NOT the real Joule loss. The PDF-OCR'd Stafl 4-157 was used in
an aborted attempt; this module's formula is the textbook-canonical
real-power form.)

Limits (verified by tests)
--------------------------

* **Small ka** (thick skin, ``a << δ``): Taylor expansion gives
  ``P → (π σ ω² μ² H_0² a^4) / 16``. This matches the elementary
  eddy-current derivation by Faraday's law plus Joule integral over a
  cylinder with uniform internal B (``B = μ_0 H_0``):

      E_phi(r) = (r / 2) ω B_0 sin(ω t),
      J_phi   = σ E_phi,
      <P/V>   = (1/2) σ E_phi^2_peak,
      P/L     = ∫_0^a (σ ω² B_0² r² / 8) 2π r dr
              = π σ ω² B_0² a^4 / 16.

* **Large ka** (thin skin, ``a >> δ``): asymptotically the Kelvin-
  function ratio approaches ``-1/√2`` and ``P → π a H_0² Re(Z_s)`` with
  the planar surface resistance ``Re(Z_s) = sqrt(omega mu / (2 sigma))``.

References
----------
Smythe W. R., "Static and Dynamic Electricity", McGraw-Hill (3rd ed.
1968), §11.07 -- cylinder in uniform AC magnetic field via Bessel
functions.

Landau L. D., Lifshitz E. M., "Electrodynamics of Continuous Media"
(Pergamon, 2nd ed. 1984), §59 -- skin effect in cylindrical / spherical
conductors.

Jackson J. D., "Classical Electrodynamics" (Wiley, 3rd ed. 1999),
§5.18 -- diffusion of fields in conductors.
"""

from __future__ import annotations

import math

from scipy import special


MU_0 = 4.0e-7 * math.pi


def _ka(a: float, omega: float, sigma: float, mu_r: float) -> float:
    """Stafl wave-number times radius, ``ka = a sqrt(omega sigma mu)``."""
    if a <= 0 or omega <= 0 or sigma <= 0 or mu_r < 1:
        raise ValueError(
            f"a>0, omega>0, sigma>0, mu_r>=1 required; got a={a}, "
            f"omega={omega}, sigma={sigma}, mu_r={mu_r}"
        )
    return a * math.sqrt(omega * sigma * mu_r * MU_0)


def cylinder_axial_eddy_loss(
    H_0: float,
    a: float,
    omega: float,
    sigma: float,
    mu_r: float = 1.0,
) -> float:
    """Time-averaged Joule loss per unit length of an infinite conductive
    cylinder of radius ``a`` placed inside a long solenoid that produces
    axial surface field amplitude ``H_0`` on the cylinder boundary.

    Closed-form Bessel solution (see module docstring for the derivation
    and the textbook references). For ``mu_r = 1`` this is the
    canonical induction-furnace work-piece formula.

    .. math::

        P = \\frac{\\pi}{\\kappa} H_0^2 (ka)
            \\frac{\\mathrm{ber}(ka) \\mathrm{ber}'(ka)
                  + \\mathrm{bei}(ka) \\mathrm{bei}'(ka)}
                 {\\mathrm{ber}(ka)^2 + \\mathrm{bei}(ka)^2},

    with ``k = sqrt(omega sigma mu)``. See ``cylinder_axial_eddy_loss_small_ka``
    and ``cylinder_axial_eddy_loss_thin_skin`` for the asymptotic limits.

    Parameters
    ----------
    H_0 : float
        Peak amplitude of the axial field at the cylinder surface
        [A/m]. Equal to ``NI = N × I`` for a long solenoid of
        ``N`` turns per unit length carrying peak current ``I``.
    a : float
        Cylinder radius [m].
    omega : float
        Angular frequency [rad/s].
    sigma : float
        Cylinder conductivity [S/m].
    mu_r : float
        Cylinder relative permeability (default 1).

    Returns
    -------
    P : float
        Time-averaged Joule loss per unit length [W / m].
    """
    if H_0 < 0:
        raise ValueError(f"H_0>=0 required; got {H_0}")
    ka = _ka(a, omega, sigma, mu_r)
    ber = special.ber(ka)
    bei = special.bei(ka)
    berp = special.berp(ka)
    beip = special.beip(ka)
    denom = ber * ber + bei * bei
    if denom == 0.0:
        return 0.0
    num = ber * berp + bei * beip
    return math.pi / sigma * H_0 ** 2 * ka * num / denom


def cylinder_axial_eddy_loss_small_ka(
    H_0: float,
    a: float,
    omega: float,
    sigma: float,
    mu_r: float = 1.0,
) -> float:
    """Small-``ka`` (thick-skin / low-frequency) asymptote of
    ``cylinder_axial_eddy_loss``.

    .. math::

        P_{small} = \\frac{\\pi}{16} \\sigma \\omega^2 \\mu^2 H_0^2 a^4.

    Derivation: in the thick-skin limit ``B`` is uniform inside the
    cylinder (no shielding), the induced ``E_phi(r) = (r/2) omega B``
    by Faraday, the Joule integrand ``sigma E^2 / 2`` integrated over
    the cross-section gives the formula above. This matches
    ``cylinder_axial_eddy_loss`` for ``ka << 1``.
    """
    if H_0 < 0 or a <= 0 or omega < 0 or sigma <= 0 or mu_r < 1:
        raise ValueError(
            f"H_0>=0, a>0, omega>=0, sigma>0, mu_r>=1 required; got "
            f"H_0={H_0}, a={a}, omega={omega}, sigma={sigma}, mu_r={mu_r}"
        )
    mu = mu_r * MU_0
    return math.pi / 16.0 * sigma * omega ** 2 * mu ** 2 * H_0 ** 2 * a ** 4


def cylinder_axial_eddy_loss_thin_skin(
    H_0: float,
    a: float,
    omega: float,
    sigma: float,
    mu_r: float = 1.0,
) -> float:
    """Large-``ka`` (thin-skin / high-frequency) asymptote of
    ``cylinder_axial_eddy_loss``.

    .. math::

        P_{thin} = \\pi a H_0^2 \\, \\mathrm{Re}(Z_s),
        \\quad \\mathrm{Re}(Z_s) = \\sqrt{\\omega \\mu / (2 \\sigma)}.

    Equivalent to wrapping the planar surface impedance
    :math:`Z_s = (1 + j) / (\\sigma \\delta)` around the cylinder
    perimeter ``2 pi a`` with surface current amplitude ``H_0`` per
    unit length and computing the time-averaged surface dissipation.
    Matches ``cylinder_axial_eddy_loss`` for ``ka >> 1``.
    """
    if H_0 < 0 or a <= 0 or omega < 0 or sigma <= 0 or mu_r < 1:
        raise ValueError(
            f"H_0>=0, a>0, omega>=0, sigma>0, mu_r>=1 required; got "
            f"H_0={H_0}, a={a}, omega={omega}, sigma={sigma}, mu_r={mu_r}"
        )
    mu = mu_r * MU_0
    Re_Zs = math.sqrt(omega * mu / (2.0 * sigma))
    return math.pi * a * H_0 ** 2 * Re_Zs
