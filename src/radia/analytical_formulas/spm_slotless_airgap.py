"""Air-gap field of a slotless surface-magnet (SPM) machine, internal rotor.

2-D subdomain solution of Zhu, Howe, Bolte and Ackermann (1993) for radially
magnetised magnets: rotor yoke at R_r and stator bore at R_s are infinitely
permeable, magnets fill R_r < r < R_m over a pole-arc fraction alpha_p, the
recoil permeability is mu_r and the remanence B_rem (M = B_rem / mu0).
For n odd and k = n p (p pole pairs, k != 1):

    M_n  = 2 M alpha_p sin(n pi alpha_p / 2) / (n pi alpha_p / 2)
    K_B  = mu0 M_n / mu_r * k / (k^2 - 1)
           * [(k - 1) + 2 (R_r/R_m)^(k+1) - (k + 1)(R_r/R_m)^(2k)]
           / [(mu_r+1)/mu_r (1 - (R_r/R_s)^(2k)) - (mu_r-1)/mu_r ((R_m/R_s)^(2k) - (R_r/R_m)^(2k))]
    B_r  = sum K_B [ (r/R_s)^(k-1) (R_m/R_s)^(k+1) + (R_m/r)^(k+1)] cos(k theta)
    B_t  = sum K_B [-(r/R_s)^(k-1) (R_m/R_s)^(k+1) + (R_m/r)^(k+1)] sin(k theta)

theta = 0 is the centre of a north pole.  Valid in the air gap R_m <= r <= R_s.

Reference: Z. Q. Zhu, D. Howe, E. Bolte, B. Ackermann, "Instantaneous magnetic
field distribution in brushless permanent magnet dc motors, Part I: open-circuit
field", IEEE Trans. Magn. 29(1):124-135, 1993.

Checked against a 2-D NGSolve A_z FEM of the same idealised machine (natural
boundary on both iron surfaces): the maximum air-gap difference falls from 7 %
to 0.13 % of max |B_r| under p/h refinement.
"""
from __future__ import annotations

import numpy as np

MU0 = 4.0e-7 * np.pi

__all__ = ["slotless_spm_airgap_field"]


def slotless_spm_airgap_field(r_m, theta_rad, pole_pairs: int, remanence_t: float,
                              recoil_permeability: float, rotor_radius_m: float,
                              magnet_radius_m: float, stator_radius_m: float,
                              pole_arc_fraction: float, max_harmonic: int = 199):
    """Return (B_r, B_theta) [T] at air-gap radius r_m for angles theta_rad.

    Raises:
        ValueError: p < 2 (the k = 1 term needs a separate limit and is not
            implemented), radii out of order, r outside the air gap,
            pole_arc_fraction outside (0, 1], non-positive mu_r, non-integral
            pole_pairs or max_harmonic, or even/negative max_harmonic.
    """
    if isinstance(pole_pairs, bool) or not float(pole_pairs).is_integer():
        raise ValueError(f"pole_pairs must be an integer, got {pole_pairs!r}")
    if isinstance(max_harmonic, bool) or not float(max_harmonic).is_integer():
        raise ValueError(f"max_harmonic must be an integer, got {max_harmonic!r}")
    p = int(pole_pairs)
    max_harmonic = int(max_harmonic)
    Rr, Rm, Rs, ap, mur = (float(rotor_radius_m), float(magnet_radius_m),
                           float(stator_radius_m), float(pole_arc_fraction),
                           float(recoil_permeability))
    r = float(r_m)
    th = np.asarray(theta_rad, float)
    if p < 2:
        raise ValueError("pole_pairs must be >= 2 (k = 1 harmonic not implemented)")
    if not 0.0 < Rr < Rm < Rs:
        raise ValueError(f"radii must satisfy 0 < R_r < R_m < R_s, got {Rr}, {Rm}, {Rs}")
    if not Rm <= r <= Rs:
        raise ValueError(f"r = {r} is outside the air gap [{Rm}, {Rs}]")
    if not 0.0 < ap <= 1.0:
        raise ValueError("pole_arc_fraction must be in (0, 1]")
    if not mur > 0.0:
        raise ValueError("recoil_permeability must be positive")
    if max_harmonic < 1 or max_harmonic % 2 == 0:
        raise ValueError("max_harmonic must be a positive odd integer")
    if not np.all(np.isfinite(th)) or not np.isfinite(remanence_t):
        raise ValueError("theta and remanence must be finite")
    M = remanence_t / MU0
    br = np.zeros_like(th)
    bt = np.zeros_like(th)
    for n in range(1, max_harmonic + 1, 2):
        k = n * p
        x = n * np.pi * ap / 2.0
        Mn = 2.0 * M * ap * np.sin(x) / x
        num = (k - 1) + 2.0 * (Rr / Rm) ** (k + 1) - (k + 1) * (Rr / Rm) ** (2 * k)
        den = ((mur + 1) / mur * (1 - (Rr / Rs) ** (2 * k))
               - (mur - 1) / mur * ((Rm / Rs) ** (2 * k) - (Rr / Rm) ** (2 * k)))
        KB = MU0 * Mn / mur * k / (k * k - 1.0) * num / den
        a = (r / Rs) ** (k - 1) * (Rm / Rs) ** (k + 1)
        b = (Rm / r) ** (k + 1)
        br += KB * (a + b) * np.cos(k * th)
        bt += KB * (b - a) * np.sin(k * th)
    return br, bt
