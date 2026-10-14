"""
compare_filament_vs_centerline.py

Tier A sanity-check: compare external B field from
  (a) filament model (nw x nh > 1, with Tier A distribution)
  (b) centerline-only model (nw = nh = 1, all current on centerline)

The difference quantifies how much the "circulating current" + skin effect
distribution changes the external field vs. the naive centerline Biot-Savart.

Physical expectation
--------------------
- Far from coil (distance >> cross-section): (a) and (b) must agree;
  the field depends only on the total current, not distribution.
- Near the coil surface: (a) can differ by a few % due to current
  redistribution. At a future workpiece position this difference matters.

This script is NOT the BEM validation (see compare_tier_a_vs_ngbem.py
for that, which is a template requiring NGSolve setup).

Usage
-----
    python compare_filament_vs_centerline.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from radia.coil_builder import CoilBuilder


MU_0 = 4.0 * np.pi * 1e-7


def biot_savart_segments(segments_with_currents, probe_points):
    """Biot-Savart sum of straight-line segments carrying (complex) currents.

    Args:
        segments_with_currents: list of (p1, p2, I_complex)
            p1, p2: length-3 tuples or ndarrays [m]
            I_complex: complex current [A]
        probe_points: ndarray (M, 3) probe positions [m]

    Returns:
        B: ndarray (M, 3) complex B field [T]. Real part for DC.
    """
    M = probe_points.shape[0]
    B = np.zeros((M, 3), dtype=complex)
    factor = MU_0 / (4 * np.pi)
    for p1, p2, I in segments_with_currents:
        p1 = np.asarray(p1, dtype=float)
        p2 = np.asarray(p2, dtype=float)
        a = probe_points - p1  # (M, 3)
        b = probe_points - p2  # (M, 3)
        la = np.linalg.norm(a, axis=1)  # (M,)
        lb = np.linalg.norm(b, axis=1)
        ab = p2 - p1            # (3,)
        lab = np.linalg.norm(ab)
        if lab < 1e-18:
            continue
        # Field of finite line segment (closed form)
        #   B = (mu_0 I / (4 pi)) * (a x b) / |a x b|^2
        #                        * ((a.ab)/la - (b.ab)/lb)
        cross = np.cross(a, b)                       # (M, 3)
        cross2 = np.sum(cross * cross, axis=1)       # (M,)
        # Guard against probes on the line
        safe = cross2 > 1e-24
        a_dot_ab = a @ ab                            # (M,)
        b_dot_ab = b @ ab                            # (M,)
        scale = np.zeros(M)
        scale[safe] = factor * (
            a_dot_ab[safe] / la[safe] - b_dot_ab[safe] / lb[safe]
        ) / cross2[safe]
        B += I * (cross * scale[:, None])
    return B


def coil_to_segments(coil, nw, nh, frequency, sigma):
    """CoilBuilder -> list of (p1, p2, I_complex) for Biot-Savart."""
    paths, currents = coil.to_filaments(nw=nw, nh=nh, frequency=frequency,
                                        sigma=sigma)
    flat = []
    for path_k, I_k in zip(paths, currents):
        for p1, p2 in path_k:
            flat.append((p1, p2, I_k))
    return flat


def main():
    mm = 1e-3

    # Same racetrack as smoke test
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * mm, height=5 * mm)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180))

    f = 7e3
    sigma = 5.8e7

    # Probe line: run parallel to the coil's long axis, offset in -Z by
    # various distances (simulating workpiece positions).
    # Racetrack lies in x-y plane (y is flow along straight, arcs in x-y).
    # Put probe line at y = 20 mm (middle of one straight), vary z below.
    probe_z = np.linspace(-50 * mm, -2.6 * mm, 25)  # 2.6 mm ~ just below 5/2
    probes = np.column_stack([
        np.zeros_like(probe_z),
        np.full_like(probe_z, 20 * mm),
        probe_z,
    ])

    # (a) filament model
    segs_a = coil_to_segments(coil, nw=7, nh=7, frequency=f, sigma=sigma)
    B_a = biot_savart_segments(segs_a, probes)

    # (b) centerline-only model (same total current, but concentrated)
    segs_b = coil_to_segments(coil, nw=1, nh=1, frequency=f, sigma=sigma)
    B_b = biot_savart_segments(segs_b, probes)

    # Magnitudes
    mag_a = np.linalg.norm(np.abs(B_a), axis=1)
    mag_b = np.linalg.norm(np.abs(B_b), axis=1)
    rel_diff = (mag_a - mag_b) / mag_b

    # --- Report ---
    print(f"{'|z| (mm)':>10}  {'|B|_filament (T)':>18}  "
          f"{'|B|_centerline (T)':>18}  {'rel diff':>10}")
    for z, ma, mb, rd in zip(probe_z, mag_a, mag_b, rel_diff):
        print(f"{-z*1e3:>10.2f}  {ma:>18.6e}  {mb:>18.6e}  {rd*100:>9.3f}%")

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(-probe_z * 1e3, mag_a, 'o-', label='filament (7x7)')
    ax1.plot(-probe_z * 1e3, mag_b, 's--', label='centerline')
    ax1.set_xlabel('distance below coil (mm)')
    ax1.set_ylabel('|B| (T)')
    ax1.set_title(f'|B| vs probe distance @ {f/1e3:g} kHz')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    ax2.plot(-probe_z * 1e3, rel_diff * 100, 'o-', color='C3')
    ax2.axhline(0, color='k', lw=0.5)
    ax2.set_xlabel('distance below coil (mm)')
    ax2.set_ylabel('(|B|_fil - |B|_cl) / |B|_cl  [%]')
    ax2.set_title('filament vs centerline')
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, 'out_filament_vs_centerline.png')
    fig.savefig(out, dpi=120)
    print(f"\nSaved: {out}")


if __name__ == '__main__':
    main()
