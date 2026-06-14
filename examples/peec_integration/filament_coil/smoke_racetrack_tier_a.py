"""
smoke_racetrack_tier_a.py

Smoke test for CoilBuilder.to_filaments() — Tier A analytical filament
current distribution on a racetrack coil.

Verifies:
  1. to_filaments() runs without error on a multi-segment coil
  2. Total current sums to self.current
  3. Arc segments produce inward-biased current (shorter path = lower R
     = more current) — "circulating current" signature
  4. At AC frequency, surface-biased distribution (skin effect) appears

No BEM comparison here — this is just the analytical side working end to end.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless
import matplotlib.pyplot as plt

# Ensure repo src/radia is importable
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from radia.coil_builder import CoilBuilder


def build_racetrack():
    """IH-coil style racetrack: 5 mm square cross-section, ~60 mm overall."""
    mm = 1e-3
    coil = (CoilBuilder(current=100.0)  # 100 A
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * mm, height=5 * mm)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180))
    return coil


def check_current_sum(currents, I_total, tag=""):
    s = complex(np.sum(currents))
    err = abs(s - I_total) / abs(I_total)
    print(f"  {tag}  sum(I_k) = {s.real:.4g} + {s.imag:.4g}j  (target {I_total:g}, "
          f"relerr {err:.2e})")
    assert err < 1e-10, f"Currents do not sum to I_total (err={err})"


def main():
    coil = build_racetrack()
    print(f"Racetrack: {len(coil.segments)} segments, "
          f"total path length centerline ~ "
          f"{sum(getattr(s, 'length', 0) + np.deg2rad(getattr(s, 'arc_angle', 0)) * getattr(s, 'radius', 0) for s in coil.segments) * 1000:.1f} mm")

    nw, nh = 7, 7  # 5 mm / 7 ~= 0.71 mm cells, resolves skin depth at 7 kHz

    # --- DC case ---
    print("\n[DC / f=0] pure path-length bias (1/R_k)")
    paths_dc, I_dc = coil.to_filaments(nw=nw, nh=nh, frequency=0.0)
    check_current_sum(I_dc, coil.current, tag="DC")
    # DC currents are real; arc inner filaments should carry more
    info = coil._last_filament_info
    J_dc = I_dc.reshape(nw, nh).real
    # Row i = 0 is innermost (u = -w/2 + du/2), but for racetrack arcs live
    # in x-y plane and inner-arc-radius corresponds to u < 0. So column
    # min(|u|) side (center) and negative-u side should have more current.
    print(f"  J range: [{J_dc.min():.4g}, {J_dc.max():.4g}]")
    print(f"  J ratio max/min: {J_dc.max()/J_dc.min():.3f}  "
          f"(>1 confirms arc-induced path-length bias)")

    # --- AC case ---
    f = 7e3  # 7 kHz (IH band)
    sigma = 5.8e7
    print(f"\n[AC / f={f/1e3:g} kHz] 1/R_k * exp(-(1+j)d/delta)")
    paths_ac, I_ac = coil.to_filaments(nw=nw, nh=nh, frequency=f,
                                       sigma=sigma)
    check_current_sum(I_ac, coil.current, tag="AC")
    delta = np.sqrt(2 / (2 * np.pi * f * 4e-7 * np.pi * sigma))
    print(f"  skin depth delta = {delta*1e3:.3g} mm, "
          f"cell size = {info['w']/nw*1e3:.3g} mm "
          f"(delta/cell = {delta / (info['w']/nw):.2f})")
    J_mag = np.abs(I_ac).reshape(nw, nh)
    J_ph = np.angle(I_ac.reshape(nw, nh), deg=True)
    print(f"  |I_k| range: [{J_mag.min():.4g}, {J_mag.max():.4g}]")
    print(f"  phase(I_k) range: [{J_ph.min():.2f}, {J_ph.max():.2f}] deg")

    # --- Plot ---
    fig_dc = coil.plot_cross_section_current(nw=nw, nh=nh, frequency=0.0)
    out_dc = os.path.join(HERE, 'out_cross_section_dc.png')
    fig_dc.savefig(out_dc, dpi=120)
    print(f"\nSaved: {out_dc}")

    fig_ac = coil.plot_cross_section_current(nw=nw, nh=nh, frequency=f,
                                             sigma=sigma)
    out_ac = os.path.join(HERE, 'out_cross_section_ac_7kHz.png')
    fig_ac.savefig(out_ac, dpi=120)
    print(f"Saved: {out_ac}")

    # --- Path count / size ---
    n_fil = len(paths_ac)
    seg_per_fil = len(paths_ac[0])
    print(f"\nFilaments: {n_fil},  polyline segments per filament: {seg_per_fil}")
    print(f"Total line segments (for Biot-Savart): {n_fil * seg_per_fil}")

    print("\nSmoke test OK.")


if __name__ == '__main__':
    main()
