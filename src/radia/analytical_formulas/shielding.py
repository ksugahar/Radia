"""Static and AC-thin-shell shielding factor of a magnetic / conducting shell.

Part 1, eq 23-24 (static, magnetic shell), Part 6 §2 eq 5 (interior
field of magnetic spherical shell) and Part 8 §2 eq 4-5 (AC, thin
conducting shell) of the Wakao-Igarashi-Fujiwara-Kameari series.

References
----------
Wakao S., Igarashi H., Fujiwara K., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 1)", IEE Japan SA-02-28 / RM-02-64 (2002).
  -> static magnetic-shell formulas (eq 23-24).

Wakao S., Fujiwara K., Tokumasu T., Kameari A.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 6)", IEE Japan SA-05-15 / RM-05-15 (2005).
  -> internal field of a magnetic spherical shell (eq 5).

Matsuo T.,
  "Useful Formulas of Analytical Integration in Electromagnetic Field
  Computations (Part 8)", IEE Japan SA-06-76 / RM-06-78 (2006).
  -> thin-shell AC shielding formulas (eq 4-5).

Rikitake T., "Magnetic and Electromagnetic Shielding",
  Terra Scientific Publishing Co. (1987). -- classical reference for
  finite cylinders, multilayer shells, and high-frequency limits.

Hoburg J. F., "Principles of Quasistatic Magnetic Shielding with
  Cylindrical and Spherical Shields", IEEE Trans. EMC 37 (1995),
  pp. 574-579.

Definitions
-----------
The shielding factor ``S`` is the ratio of the magnitude of the field
inside the shell to the magnitude of the applied uniform external field:

    S = |H_inside| / |H_applied|

A perfect shield gives ``S = 0``; no shield gives ``S = 1``. The smaller
``S`` the better the shielding.

Cylindrical shell (infinite length, transverse field, eq 23)
------------------------------------------------------------
Inner radius ``a``, outer radius ``b``, relative permeability ``mu_r``,
applied field perpendicular to the cylinder axis:

    S_cyl = 4 mu_r / [(mu_r + 1)**2 - (a/b)**2 (mu_r - 1)**2]

Spherical shell (eq 24)
-----------------------
Inner radius ``a``, outer radius ``b``, relative permeability ``mu_r``:

    S_sph = 9 mu_r
            / [(mu_r + 2)(2 mu_r + 1) - 2 (a/b)**3 (mu_r - 1)**2]

Both formulas approach ``1`` as ``mu_r -> 1`` (no shielding) and
``-> 0`` as ``mu_r -> infinity`` with finite wall thickness, as
expected. For thin walls ``b - a = t << a`` the high-mu_r asymptote is

    S_cyl ~ 1 / (mu_r * t / a),   S_sph ~ 3 / (2 mu_r * t / a).

Notes
-----
The cylindrical formula assumes infinite length; finite cylinders and
caps are discussed in Rikitake [13]; for a quick estimate the infinite
formula is conservative (under-predicts shielding compared to a sealed
cylinder).
"""

from __future__ import annotations

import math

MU_0 = 4.0e-7 * math.pi


def shielding_factor_cylinder(a: float, b: float, mu_r: float) -> float:
    """Shielding factor of an infinite cylindrical magnetic shell (eq 23).

    Parameters
    ----------
    a : float
        Inner radius.
    b : float
        Outer radius (``b > a``).
    mu_r : float
        Relative permeability of the shell material (``>= 1``).

    Returns
    -------
    S : float
        Shielding factor ``|H_inside| / |H_applied|``.
    """
    if a <= 0 or b <= a:
        raise ValueError(f"require 0 < a < b, got a={a}, b={b}")
    if mu_r < 1.0:
        raise ValueError(f"mu_r must be >= 1, got {mu_r}")
    return 4.0 * mu_r / (
        (mu_r + 1.0) ** 2 - (a / b) ** 2 * (mu_r - 1.0) ** 2
    )


def shielding_factor_sphere(a: float, b: float, mu_r: float) -> float:
    """Shielding factor of a spherical magnetic shell (eq 24).

    Parameters
    ----------
    a : float
        Inner radius.
    b : float
        Outer radius (``b > a``).
    mu_r : float
        Relative permeability.

    Returns
    -------
    S : float
        Shielding factor ``|H_inside| / |H_applied|``.
    """
    if a <= 0 or b <= a:
        raise ValueError(f"require 0 < a < b, got a={a}, b={b}")
    if mu_r < 1.0:
        raise ValueError(f"mu_r must be >= 1, got {mu_r}")
    return 9.0 * mu_r / (
        (mu_r + 2.0) * (2.0 * mu_r + 1.0)
        - 2.0 * (a / b) ** 3 * (mu_r - 1.0) ** 2
    )


# ---------------------------------------------------------------------------
# Part 8 §2 -- AC thin-shell conductors
# ---------------------------------------------------------------------------


def shielding_factor_sphere_thin_ac(
    a: float,
    delta: float,
    sigma: float,
    omega: float,
) -> complex:
    """AC shielding factor of a thin spherical conducting shell (Part 8 eq 4).

        S = 1 / (1 + j omega mu_0 sigma a Delta / 3)

    Valid limits:

    * thin wall:    ``Delta << a``
    * quasi-static: ``Delta << skin_depth = sqrt(2 / (omega mu_0 sigma))``

    The shell is non-magnetic (mu_r = 1). Inside the shell the field
    is uniform and equal to ``S * H_applied`` (complex phasor; magnitude
    < 1 indicates shielding, phase indicates lag).

    Parameters
    ----------
    a : float
        Shell radius [m].
    delta : float
        Shell wall thickness [m]; ``Delta`` in the PDF.
    sigma : float
        Conductivity [S/m].
    omega : float
        Angular frequency [rad/s].

    Returns
    -------
    S : complex
        AC shielding factor (complex; ``|S|`` is the amplitude
        attenuation, ``arg(S)`` is the phase shift).
    """
    if a <= 0 or delta <= 0:
        raise ValueError(f"a, delta must be > 0; got a={a}, delta={delta}")
    if sigma < 0 or omega < 0:
        raise ValueError(
            f"sigma, omega must be >= 0; got sigma={sigma}, omega={omega}"
        )
    return 1.0 / (1.0 + 1j * omega * MU_0 * sigma * a * delta / 3.0)


def shielding_factor_cylinder_thin_ac(
    a: float,
    delta: float,
    sigma: float,
    omega: float,
) -> complex:
    """AC shielding factor of an infinite thin cylindrical conducting shell
    in a transverse field (Part 8 eq 5).

        S = 1 / (1 + j omega mu_0 sigma a Delta / 2)

    Same thin-wall + quasi-static assumptions as
    :func:`shielding_factor_sphere_thin_ac`. The "/ 2" factor (vs "/ 3"
    for the sphere) reflects the geometry-only image-dipole strength.
    """
    if a <= 0 or delta <= 0:
        raise ValueError(f"a, delta must be > 0; got a={a}, delta={delta}")
    if sigma < 0 or omega < 0:
        raise ValueError(
            f"sigma, omega must be >= 0; got sigma={sigma}, omega={omega}"
        )
    return 1.0 / (1.0 + 1j * omega * MU_0 * sigma * a * delta / 2.0)


# ---------------------------------------------------------------------------
# Part 6 §2 -- internal field inside a magnetic spherical shell
# ---------------------------------------------------------------------------


def spherical_shell_internal_field(
    H0: float,
    a: float,
    b: float,
    mu_r: float,
) -> dict:
    """Field components inside a magnetic spherical shell in a uniform
    external field (Part 6 §2, eq 5).

    The hollow interior holds a uniform field ``H_h``; the shell wall
    has a uniform internal H ``H_s`` and uniform magnetization ``M_s``;
    outside the shell the dipole moment ``M_i`` is induced. All four
    are returned as a dict.

    Parameters
    ----------
    H0 : float
        External applied field magnitude [A/m].
    a, b : float
        Inner / outer radii (b > a > 0) [m].
    mu_r : float
        Relative permeability of the shell material (>= 1).

    Returns
    -------
    dict with keys:
        ``H_hollow``     : magnitude of uniform field in hollow interior [A/m]
        ``H_shell``      : magnitude of uniform field inside shell wall [A/m]
        ``M_shell``      : magnetisation of the shell wall [A/m]
        ``M_image``      : dipole moment of the equivalent external image [A m^2]
        ``S``            : shielding factor ``H_hollow / H0`` (matches
                           :func:`shielding_factor_sphere`)
    """
    if a <= 0 or b <= a:
        raise ValueError(f"require 0 < a < b, got a={a}, b={b}")
    if mu_r < 1.0:
        raise ValueError(f"mu_r must be >= 1, got {mu_r}")
    denom = 9.0 + 2.0 * (mu_r - 1) ** 2 / mu_r * (1.0 - (a / b) ** 3)
    H_hollow = 9.0 * H0 / denom
    H_shell = (3.0 * (2.0 + mu_r) / (3.0 * mu_r)) * H_hollow
    M_shell = (mu_r - 1.0) * H_shell
    # External induced moment from boundary condition at b
    M_image = -((mu_r - 1) * (mu_r + 2) / (3.0 * mu_r)) * H_hollow * (
        4.0 / 3.0 * math.pi * b ** 3
    )
    return {
        "H_hollow": H_hollow,
        "H_shell": H_shell,
        "M_shell": M_shell,
        "M_image": M_image,
        "S": H_hollow / H0,
    }
