"""Static shielding factor of a magnetic shell in a uniform external field.

Part 1, eq 23-24 of the Wakao-Igarashi-Fujiwara-Kameari series.

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
