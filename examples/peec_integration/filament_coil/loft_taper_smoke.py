"""loft_taper_smoke.py

Smoke test for Stage 3a: non-uniform cross-section (rect loft) in the
filament-approximation PEEC pipeline.

Scenarios
---------
1. Baseline uniform 5x5 straight: filaments have equal R, DC current
   must be uniform.
2. Rect taper 5x5 -> 3x5 straight (40 mm): along-path cell area shrinks
   but shrinks uniformly across filaments. DC currents should STILL be
   uniform (all filaments see the same 1/A(s) integral). This validates
   the loft bookkeeping -- any R_k discrepancy reveals a bug in the
   per-piece area tracking.
3. AC at 7 kHz on the same taper: PEEC finds the self-consistent skin
   + proximity pattern. Expected qualitative behavior: current pushes
   toward the NARROW end (local J increases there because A shrank)
   and toward the edges of the cross-section at each s (skin effect).

No BEM-truth comparison here -- that is Phase 2 territory. This smoke
test only validates that (i) loft paths are geometrically correct,
(ii) per-piece cell_wh flow to PEEC produces a well-conditioned matrix,
(iii) DC symmetry is preserved for pure-taper lofts.
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
from radia.coil_profile import RectProfile
from peec_bundle import build_bundle_solver, filament_currents

MM = 1e-3
CU_SIGMA = 5.8e7


def _peec_currents(coil, nw, nh, frequency):
    """Run the PEEC bundle on the coil's filament discretization."""
    paths, _ = coil.to_filaments(nw=nw, nh=nh, frequency=frequency,
                                  sigma=CU_SIGMA)
    info = coil._last_filament_info
    cell_wh = info['cell_wh']

    solver, seg_of_fil, _p, _m = build_bundle_solver(
        paths, dw=None, dh=None, sigma=CU_SIGMA, cell_wh=cell_wh)
    I_branch = solver.compute_branch_currents(frequency,
                                               [complex(coil.current)])
    return filament_currents(I_branch, seg_of_fil), info


def _print_filament_summary(tag, I, coil, info):
    s = complex(np.sum(I))
    mag = np.abs(I)
    print(f"  [{tag}] N={I.size}  sum(I) = {s.real:.4f} + {s.imag:.4f}j "
          f"(target {coil.current:g})")
    print(f"         |I_k| range [{mag.min():.4f}, {mag.max():.4f}] "
          f"max/min={mag.max()/mag.min():.3f}")
    if np.iscomplexobj(I):
        ph = np.angle(I, deg=True)
        print(f"         arg(I_k) range [{ph.min():.2f}, {ph.max():.2f}] deg")


def scenario_1_uniform_straight():
    print("\n=== Scenario 1: baseline uniform 5x5 straight, DC ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * MM, height=5 * MM)
            .add_straight(length=40 * MM))
    I, info = _peec_currents(coil, nw=5, nh=5, frequency=0.0)
    _print_filament_summary('uniform DC', I, coil, info)
    spread = np.ptp(np.abs(I)) / np.mean(np.abs(I))
    print(f"         relative spread (ptp/mean) = {spread:.2e}")
    assert spread < 1e-3, "Uniform straight DC should give uniform I_k"
    print("         OK: DC uniform as expected")


def scenario_2_taper_dc():
    print("\n=== Scenario 2: rect taper 5x5 -> 3x5, 40 mm straight, DC ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * MM, height=5 * MM)
            .add_loft_straight(profile_end=RectProfile(3 * MM, 5 * MM),
                               length=40 * MM, n_sub=20))
    I, info = _peec_currents(coil, nw=5, nh=5, frequency=0.0)
    _print_filament_summary('taper DC', I, coil, info)
    spread = np.ptp(np.abs(I)) / np.mean(np.abs(I))
    print(f"         relative spread (ptp/mean) = {spread:.2e}")
    # For a pure rect->rect taper, all filaments see the SAME A(s)
    # profile along s, so R_k is identical -> DC must be uniform.
    assert spread < 1e-2, ("Pure rect taper DC should give uniform "
                            f"I_k (spread={spread:.3e}); bug in cell "
                            "area tracking")
    print("         OK: pure taper DC is uniform -> cell_wh plumbing OK")


def scenario_3_taper_ac():
    print("\n=== Scenario 3: rect taper 5x5 -> 3x5, 40 mm, AC 7 kHz ===")
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * MM, height=5 * MM)
            .add_loft_straight(profile_end=RectProfile(3 * MM, 5 * MM),
                               length=40 * MM, n_sub=20))
    I, info = _peec_currents(coil, nw=5, nh=5, frequency=7e3)
    _print_filament_summary('taper AC', I, coil, info)
    # Qualitative: max|I| on edges (skin effect), min|I| at center
    # Cross-section frame: alpha along width (nw=5), beta along height
    # (nh=5). Center filament has alpha=beta=0.5.
    J = np.abs(I).reshape(5, 5)
    mid_ij = (2, 2)  # center
    edge_ij = (0, 0)  # corner
    print(f"         |I| center[{mid_ij}]={J[mid_ij]:.4f}  "
          f"corner[{edge_ij}]={J[edge_ij]:.4f}  "
          f"ratio corner/center = {J[edge_ij]/J[mid_ij]:.2f}")
    assert J[edge_ij] > J[mid_ij], \
        "Skin effect should push current toward the corners"
    print("         OK: corner > center (skin redistribution present)")

    return I


def scenario_4_visualize(I_ac):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    J = np.abs(I_ac.reshape(5, 5))
    ph = np.angle(I_ac.reshape(5, 5), deg=True)
    im0 = axes[0].imshow(J.T, origin='lower', cmap='viridis', aspect='equal')
    axes[0].set_title('|I_k|  taper 5x5 -> 3x5, 7 kHz')
    axes[0].set_xlabel('i (width idx)')
    axes[0].set_ylabel('j (height idx)')
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(ph.T, origin='lower', cmap='twilight',
                          aspect='equal', vmin=-180, vmax=180)
    axes[1].set_title('arg(I_k) [deg]')
    axes[1].set_xlabel('i')
    fig.colorbar(im1, ax=axes[1])
    out = os.path.join(HERE, 'out_loft_taper_ac_7kHz.png')
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n  Saved: {out}")


def main():
    scenario_1_uniform_straight()
    scenario_2_taper_dc()
    I_ac = scenario_3_taper_ac()
    scenario_4_visualize(I_ac)
    print("\nAll loft-taper smoke tests PASSED")


if __name__ == '__main__':
    main()
