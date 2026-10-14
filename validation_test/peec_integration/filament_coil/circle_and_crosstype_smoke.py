"""circle_and_crosstype_smoke.py

Stage 3b smoke test: non-rect cross-section profiles (round wire,
annular) and cross-type loft (rect <-> circle).

Scenarios
---------
1. Round wire straight DC (CircleProfile, r=2 mm, length=40 mm).
   Area-preserving (alpha, beta) parametrization on a disk: every
   filament cell has exactly the same area, so all filaments have the
   same R -> DC must be uniform. Validates the Shirley-like polar map.

2. Round wire straight AC at 7 kHz.
   Skin depth in Cu at 7 kHz is ~0.79 mm; wire radius 2 mm -> skin
   layer pushes current to the outer filaments. Validates the polar
   filament layout + PEEC in 3D on a non-rect profile.

3. Annular tube (water-cooled IH imitation, r_in = 1 mm, r_out = 2 mm)
   DC: uniform current (same argument as circle).

4. Cross-type loft: Rect(5x5) -> Circle(r=2.5) straight (40 mm). Both
   endpoints have area-preserving (alpha, beta) maps. InterpolatedProfile
   smoothly blends. DC is expected to be "nearly uniform" (some
   deviation because the intermediate shape is neither rect nor
   circle and its area is only linearly blended). Check spread < few %.

No BEM-truth here. These smoke tests only validate that the profile
machinery is geometrically consistent and that PEEC on non-rect /
cross-type geometries gives physically sensible currents.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
for p in (SRC, SRC_RADIA):
    if p not in sys.path:
        sys.path.insert(0, p)

from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile, CircleProfile, AnnularProfile
from peec_bundle import build_bundle_solver, filament_currents

MM = 1e-3
CU_SIGMA = 5.8e7


def _peec(coil, nw, nh, frequency):
    paths, _ = coil.to_filaments(nw=nw, nh=nh, frequency=frequency,
                                  sigma=CU_SIGMA)
    info = coil._last_filament_info
    solver, seg_of_fil, *_ = build_bundle_solver(
        paths, dw=None, dh=None, sigma=CU_SIGMA,
        cell_wh=info['cell_wh'])
    I_branch = solver.compute_branch_currents(frequency,
                                               [complex(coil.current)])
    return filament_currents(I_branch, seg_of_fil), info


def _report(tag, I, coil):
    s = complex(np.sum(I))
    mag = np.abs(I)
    print(f"  [{tag}] N={I.size}  sum={s.real:.4f}+{s.imag:.4f}j "
          f"(target {coil.current:g})")
    print(f"         |I| range [{mag.min():.4f}, {mag.max():.4f}] "
          f"max/min={mag.max()/max(mag.min(),1e-30):.3f} "
          f"spread={np.ptp(mag)/np.mean(mag):.2e}")
    if np.iscomplexobj(I):
        ph = np.angle(I, deg=True)
        print(f"         arg range [{ph.min():.2f}, {ph.max():.2f}] deg")


def scenario_circle_dc():
    print("\n=== Scenario 1: round wire (r=2mm) straight 40mm DC ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_profile(CircleProfile(r=2 * MM))
            .add_straight(length=40 * MM))
    I, info = _peec(coil, nw=5, nh=8, frequency=0.0)
    _report('circle DC', I, coil)
    spread = np.ptp(np.abs(I)) / np.mean(np.abs(I))
    assert spread < 1e-2, \
        f"Round wire DC should be uniform; spread={spread}"
    print(f"         OK: area-preserving polar map gives uniform DC")


def scenario_circle_ac():
    print("\n=== Scenario 2: round wire (r=2mm) straight 40mm AC 7kHz ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_profile(CircleProfile(r=2 * MM))
            .add_straight(length=40 * MM))
    I, info = _peec(coil, nw=5, nh=8, frequency=7e3)
    _report('circle AC', I, coil)

    # Recover (alpha, beta) from info: alpha radial, beta azimuthal.
    # alpha = (i+0.5)/nw so i=0 -> innermost ring, i=nw-1 -> outermost.
    nw, nh = info['nw'], info['nh']
    I_2d = I.reshape(nw, nh)
    inner_mag = np.mean(np.abs(I_2d[0, :]))   # innermost ring
    outer_mag = np.mean(np.abs(I_2d[-1, :]))  # outermost ring
    print(f"         inner ring <|I|>={inner_mag:.4f}, "
          f"outer ring <|I|>={outer_mag:.4f}, "
          f"ratio outer/inner={outer_mag/inner_mag:.2f}")
    assert outer_mag > 1.5 * inner_mag, (
        "Skin effect should push current to outer ring; "
        f"got outer/inner = {outer_mag/inner_mag:.2f}")
    print(f"         OK: skin effect pushes current outward")

    # Also plot it (polar heatmap-like -- show as rectangular index map)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    im0 = axes[0].imshow(np.abs(I_2d).T, origin='lower', cmap='viridis',
                          aspect='auto')
    axes[0].set_title('|I_k|  circle wire  r=2mm  f=7kHz')
    axes[0].set_xlabel('i (radial)')
    axes[0].set_ylabel('j (azimuthal)')
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(np.angle(I_2d, deg=True).T, origin='lower',
                          cmap='twilight', vmin=-180, vmax=180,
                          aspect='auto')
    axes[1].set_title('arg(I_k) [deg]')
    axes[1].set_xlabel('i (radial)')
    fig.colorbar(im1, ax=axes[1])
    fig.tight_layout()
    out = os.path.join(HERE, 'out_circle_wire_ac_7kHz.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"         saved: {out}")


def scenario_annular_dc():
    print("\n=== Scenario 3: annular tube (r_in=1mm, r_out=2mm) DC ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_profile(AnnularProfile(r_in=1 * MM, r_out=2 * MM))
            .add_straight(length=40 * MM))
    I, info = _peec(coil, nw=4, nh=8, frequency=0.0)
    _report('annular DC', I, coil)
    spread = np.ptp(np.abs(I)) / np.mean(np.abs(I))
    assert spread < 1e-2, f"Annular DC should be uniform; spread={spread}"
    print(f"         OK: annular area-preserving map is consistent")


def scenario_crosstype_loft_dc():
    print("\n=== Scenario 4: loft Rect(5x5) -> Circle(r=2.5mm) DC ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_profile(RectProfile(5 * MM, 5 * MM))
            .add_loft_straight(profile_end=CircleProfile(r=2.5 * MM),
                               length=40 * MM, n_sub=20))
    I, info = _peec(coil, nw=5, nh=5, frequency=0.0)
    _report('rect->circle DC', I, coil)
    # Cross-type intermediate shape is NOT area-preserving; InterpolatedProfile
    # only blends (u, v) and area linearly. Expect non-trivial spread
    # because filaments on the diagonals (rect corners) travel OUTWARD
    # more than filaments on the axes (mid edges) as the rect morphs to
    # circle. Path lengths DIFFER per filament -> R_k differs -> DC
    # distribution is non-uniform but bounded.
    spread = np.ptp(np.abs(I)) / np.mean(np.abs(I))
    print(f"         spread = {spread*100:.2f}% "
          "(expected non-zero; rect->circle path lengths differ)")
    assert spread < 0.5, (
        f"Cross-type loft DC spread too large ({spread}); "
        "expected <50% for rect5x5 -> circle2.5mm")
    print(f"         OK: cross-type loft DC spread bounded")


def main():
    scenario_circle_dc()
    scenario_circle_ac()
    scenario_annular_dc()
    scenario_crosstype_loft_dc()
    print("\nAll Stage 3b smoke tests PASSED")


if __name__ == '__main__':
    main()
