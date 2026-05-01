"""Three-phase line magnetic field: triangle, planar, helical.

Reference: Part 4 §5 (SA-04-1 / RM-04-1, 2004) eq 30-44, and
           Part 5 §3 (SA-04-44 / RM-04-68, 2004) eq 11-28.

Compares the triangle and coplanar arrangement far-field formulas
against direct Biot-Savart sums and shows the exponential decay of a
helical three-phase line.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    ac_locus_axes,
    balanced_three_phase_currents,
    field_xy,
    helical_far_field_amplitude,
    planar_far_field_amplitude,
    planar_positions,
    triangle_far_field_amplitude,
    triangle_positions,
)


def main() -> None:
    a = 0.05            # 50 mm radius (triangle / spacing for planar = 50 mm)
    I_rms = 1000.0      # 1 kA RMS
    I_peak = math.sqrt(2.0) * I_rms
    rs = np.geomspace(0.1, 50.0, 30)

    # ----- Triangle vs planar ---------------------------------------
    print(f"Three-phase line: a = d = {a*1000:.0f} mm, "
          f"I_rms = {I_rms:.0f} A, I_peak = {I_peak:.1f} A")
    print()
    print(f"{'r [m]':>8} {'|B|_tri rigorous':>18} {'tri formula':>14} "
          f"{'|B|_plan rigorous':>18} {'plan formula':>14}")
    bmax_tri = np.zeros_like(rs)
    bmax_plan = np.zeros_like(rs)
    f_tri = np.zeros_like(rs)
    f_plan = np.zeros_like(rs)
    pos_tri = triangle_positions(a)
    pos_plan = planar_positions(a)
    currents = balanced_three_phase_currents(I_peak)
    for k, r in enumerate(rs):
        B_tri = field_xy(pos_tri, currents, [r, 0.0])
        bmax_tri[k], _ = ac_locus_axes([B_tri[0], B_tri[1], 0.0])
        B_plan = field_xy(pos_plan, currents, [r, 0.0])
        bmax_plan[k], _ = ac_locus_axes([B_plan[0], B_plan[1], 0.0])
        f_tri[k] = triangle_far_field_amplitude(r, I_peak, a)
        f_plan[k] = planar_far_field_amplitude(r, I_peak, a)
    for k in (0, 5, 10, 15, 20, 25, 29):
        print(f"  {rs[k]:6.2f}    {bmax_tri[k]:.4e}      {f_tri[k]:.4e}   "
              f"     {bmax_plan[k]:.4e}     {f_plan[k]:.4e}")

    # ----- Helical exp decay ----------------------------------------
    print()
    print("Helical three-phase, p = 0.5 m  (a << p << r regime):")
    print(f"{'r [m]':>8} {'|B|_far (helical)':>18}")
    for r in (0.1, 0.5, 1.0, 2.0, 4.0):
        B_helix = helical_far_field_amplitude(r, I_peak, a=a, p=0.5)
        print(f"  {r:6.2f}    {B_helix:.4e}")

    # ----- Plot -----------------------------------------------------
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(rs, bmax_tri, "-", lw=1.7, label="triangle (rigorous)")
    ax.loglog(rs, f_tri, "--", lw=1.0, label="triangle (formula)")
    ax.loglog(rs, bmax_plan, "-", lw=1.7, label="planar (rigorous)")
    ax.loglog(rs, f_plan, "--", lw=1.0, label="planar (formula)")
    rs_fine = np.geomspace(0.05, 5.0, 100)
    Bhel = np.array([helical_far_field_amplitude(r, I_peak, a, 0.5) for r in rs_fine])
    ax.loglog(rs_fine, Bhel, "-", lw=1.7,
              label="helical (p = 0.5 m, exp decay)")
    ax.set_xlabel("Observation distance r [m]")
    ax.set_ylabel("|B|_max  [T]")
    ax.set_title("Three-phase line: far-field decay")
    ax.legend()
    ax.grid(True, ls=":", alpha=0.5, which="both")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "three_phase_line_field.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
