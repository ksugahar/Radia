"""Eddy-current and AC-shielding showcase.

Combines the Part 6 / 8 / 9 extensions in a single end-to-end story:

1. **Skin depth + planar surface impedance** for a copper plate at 1 kHz.
2. **Cylindrical-conductor AC impedance** sweep (Bessel solution): the
   resistance grows from R_dc as the skin depth shrinks below the
   conductor radius.
3. **Thin-shell AC shielding** of a non-magnetic Cu sphere as a
   function of frequency.
4. **Total Joule dissipation in a thin square plate** via Gauss-Legendre
   integration of the analytic eddy-current field.
5. **Magnetic spherical shell** internal-field decomposition.
6. **Average B over a target box** from a magnetised source box, using
   the numerical-integration path.

Every numerical reference is tagged to the source equation in the
Wakao-Igarashi-Fujiwara-Kameari series (Part 6 / 8 / 9).
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    average_B_in_box,
    cylinder_ac_impedance,
    cylinder_dc_resistance,
    cylinder_internal_inductance,
    planar_surface_impedance,
    plate_eddy_dissipation,
    shielding_factor_sphere_thin_ac,
    skin_depth,
    spherical_shell_internal_field,
)
from radia.analytical_formulas.shielding import MU_0


def main() -> None:
    sigma_cu = 5.96e7
    print("=" * 72)
    print("(1) Copper @ 1 kHz: skin depth + planar surface impedance")
    print("=" * 72)
    omega = 2.0 * math.pi * 1.0e3
    delta = skin_depth(sigma_cu, omega)
    Z_s = planar_surface_impedance(sigma_cu, omega)
    print(f"  skin depth  delta   = {delta * 1e3:.3f} mm")
    print(f"  surface impedance Z_s = {Z_s:.4e} Ohm  (|Z_s|={abs(Z_s):.4e},"
          f" arg={math.degrees(math.atan2(Z_s.imag, Z_s.real)):.2f} deg)")

    print()
    print("=" * 72)
    print("(2) Cylinder Cu, a = 1 mm, AC impedance frequency sweep")
    print("=" * 72)
    a = 1e-3
    R_dc = cylinder_dc_resistance(a, sigma_cu)
    L_int = cylinder_internal_inductance()
    print(f"  R_dc   = {R_dc:.4e} Ohm/m")
    print(f"  L_int  = {L_int * 1e9:.4f} nH/m  (low-freq internal inductance)")
    print(f"  {'f [Hz]':>10}  {'R [Ohm/m]':>12}  {'X [Ohm/m]':>12}  "
          f"{'R / R_dc':>10}")
    for f in (1, 100, 1e3, 1e4, 1e5, 1e6, 1e8):
        Z = cylinder_ac_impedance(a, sigma_cu, 2 * math.pi * f)
        print(f"  {f:>10.1e}  {Z.real:>12.4e}  {Z.imag:>12.4e}  "
              f"{Z.real / R_dc:>10.2f}")

    print()
    print("=" * 72)
    print("(3) Cu spherical shell AC shielding (a = 5 cm, t = 1 mm)")
    print("=" * 72)
    print(f"  {'f [Hz]':>10}  {'|S|':>10}  {'arg(S) [deg]':>14}")
    for f in (1, 10, 60, 1e3, 1e4):
        S = shielding_factor_sphere_thin_ac(0.05, 1e-3, sigma_cu,
                                            2 * math.pi * f)
        phase_deg = math.degrees(math.atan2(S.imag, S.real))
        print(f"  {f:>10.1e}  {abs(S):>10.4f}  {phase_deg:>14.2f}")

    print()
    print("=" * 72)
    print("(4) Square Cu plate (200 mm side, 1 mm thick), dB/dt = 1 T/s")
    print("=" * 72)
    P = plate_eddy_dissipation(0.10, 0.10, 1e-3, sigma_cu, 1.0)
    print(f"  total instantaneous Joule loss P = {P:.4f} W")
    P10 = plate_eddy_dissipation(0.10, 0.10, 1e-3, sigma_cu, 10.0)
    print(f"  at dB/dt = 10 T/s              = {P10:.2f} W "
          f"(10 x faster ramp -> 100 x dissipation, OK)")

    print()
    print("=" * 72)
    print("(5) Magnetic spherical shell (a = 5 cm, b = 6 cm, mu_r = 1000)")
    print("=" * 72)
    out = spherical_shell_internal_field(1.0, 0.05, 0.06, 1000.0)
    print(f"  H_applied   = 1.0 A/m (unit drive)")
    print(f"  H_hollow    = {out['H_hollow']:.4e} A/m")
    print(f"  H_shell     = {out['H_shell']:.4e} A/m")
    print(f"  M_shell     = {out['M_shell']:.4e} A/m")
    print(f"  M_image     = {out['M_image']:.4e} A m^2")
    print(f"  S = H_hollow / H_applied = {out['S']:.4e}")

    print()
    print("=" * 72)
    print("(6) <B> over a 2 mm cube target offset 5 cm in x from a")
    print("    20 x 20 x 10 mm magnetised cube (NdFeB-like, M_x = 9.55e5 A/m)")
    print("=" * 72)
    src_min = np.array([-0.01, -0.01, -0.005])
    src_max = np.array([+0.01, +0.01, +0.005])
    target_centre = np.array([0.05, 0.0, 0.0])
    half_t = np.array([0.001, 0.001, 0.001])
    B_avg = average_B_in_box(
        [9.55e5, 0.0, 0.0],
        src_min, src_max,
        target_centre - half_t, target_centre + half_t,
        n_quad=6,
    )
    # Compare to dipole far-field
    m = 9.55e5 * np.prod(src_max - src_min)
    d = 0.05
    B_dipole_x = (MU_0 / (4 * math.pi)) * 2 * m / d ** 3
    print(f"  <B>_target = {B_avg}  T")
    print(f"  dipole far-field B_x at d = 5 cm = {B_dipole_x:.4e} T")
    print(f"  near-field ratio = {B_avg[0] / B_dipole_x:.4f}")


if __name__ == "__main__":
    main()
