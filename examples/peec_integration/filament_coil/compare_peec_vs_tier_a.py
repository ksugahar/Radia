"""
compare_peec_vs_tier_a.py

Phase 1b validation: full filament-PEEC (R + jwM with full mutual coupling
solved via MNA) as AC truth versus the Tier A analytical filament model.

Both methods use the SAME filament discretization from
CoilBuilder.to_filaments(nw, nh, ...). The only difference is how per-
filament currents are assigned:

  Tier A (analytical): I_k = I_total * w_k / sum(w_j),
                       w_k = (1/R_k) * exp(-(1+j) d_k / delta)
                       -- captures path-length bias + a priori skin factor.
                       -- no coupling between filaments.

  PEEC (AC truth):     All N filaments wired in parallel between port+ and
                       port-. Solve (R + jwM) I = V * 1, sum(I) = I_total
                       via MNA. Captures the full proximity coupling
                       between filaments (M_kj off-diagonal terms) and
                       therefore the true AC current distribution.

If the two agree to within a few percent at IH frequencies with no
nearby conductor, Tier A is a validated low-cost model for the LLM-
optimizer inner loop.

Runs on the same racetrack coil as smoke_racetrack_tier_a.py.
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
from peec_bundle import build_bundle_solver as _build_bundle_solver_lib


# ---------------------------------------------------------------------------
# Bundle-solver helper  (thin wrapper around the library implementation
# in src/radia/peec_bundle.py, kept here for backward compatibility with
# callers importing build_bundle_solver from this file)
# ---------------------------------------------------------------------------

build_bundle_solver = _build_bundle_solver_lib


# ---------------------------------------------------------------------------
# Racetrack test
# ---------------------------------------------------------------------------

def build_racetrack():
    """Same geometry as smoke_racetrack_tier_a.py."""
    mm = 1e-3
    coil = (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * mm, height=5 * mm)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180)
            .add_straight(length=40 * mm)
            .add_arc(radius=15 * mm, arc_angle=180))
    return coil


def filament_currents_from_branch(I_branch, seg_indices):
    """Aggregate per-filament current from the per-segment branch current array.

    Since each filament is a pure series chain (no taps), all segments of a
    filament carry the same current. We average to suppress numerical noise.
    """
    N = len(seg_indices)
    I_fil = np.zeros(N, dtype=complex)
    for k in range(N):
        I_fil[k] = np.mean(I_branch[seg_indices[k]])
    return I_fil


def compare(nw, nh, frequency, sigma, n_arc, label):
    print(f"\n{'=' * 70}")
    print(f"Compare  nw={nw} nh={nh}  f={frequency:g} Hz  n_arc={n_arc}  ({label})")
    print('=' * 70)

    coil = build_racetrack()
    I_total = coil.current

    # Tier A
    paths, I_tier_a = coil.to_filaments(nw=nw, nh=nh, frequency=frequency,
                                        sigma=sigma, n_arc=n_arc)
    info = coil._last_filament_info
    dw = info['w'] / nw
    dh = info['h'] / nh

    # PEEC bundle
    solver, seg_indices, port_plus, port_minus = build_bundle_solver(
        paths, dw, dh, sigma)

    print(f"  Filaments (N): {len(paths)}")
    print(f"  PEEC n_loop:   {solver.n_loop}")
    print(f"  Segs/filament: {len(seg_indices[0])}")

    I_branch = solver.compute_branch_currents(frequency, [complex(I_total)])
    I_peec = filament_currents_from_branch(I_branch, seg_indices)

    # Sanity: sum of filament currents == I_total
    s_ta = complex(np.sum(I_tier_a))
    s_pe = complex(np.sum(I_peec))
    print(f"  Sum(I_TierA) = {s_ta.real:.4f} + {s_ta.imag:.4f}j  "
          f"(target {I_total})")
    print(f"  Sum(I_PEEC)  = {s_pe.real:.4f} + {s_pe.imag:.4f}j  "
          f"(target {I_total})")

    # Per-filament comparison
    if np.iscomplexobj(I_tier_a):
        I_ta = I_tier_a
    else:
        I_ta = I_tier_a.astype(complex)

    diff = I_peec - I_ta
    rel = np.max(np.abs(diff)) / np.max(np.abs(I_ta))
    print(f"  max |I_peec - I_ta| / max |I_ta| = {rel:.3%}")
    print(f"  mean rel diff                    = "
          f"{np.mean(np.abs(diff)) / np.mean(np.abs(I_ta)):.3%}")

    # |I| range
    print(f"  |I_TierA| range: [{np.min(np.abs(I_ta)):.4g}, "
          f"{np.max(np.abs(I_ta)):.4g}]")
    print(f"  |I_PEEC|  range: [{np.min(np.abs(I_peec)):.4g}, "
          f"{np.max(np.abs(I_peec)):.4g}]")

    if frequency > 0:
        ph_ta = np.angle(I_ta, deg=True)
        ph_pe = np.angle(I_peec, deg=True)
        print(f"  phase(I_TierA) range: [{ph_ta.min():.2f}, {ph_ta.max():.2f}] deg")
        print(f"  phase(I_PEEC)  range: [{ph_pe.min():.2f}, {ph_pe.max():.2f}] deg")

    return I_ta, I_peec, info


def plot_cross_section(I_ta, I_pe, nw, nh, info, out_path, title):
    w, h = info['w'], info['h']
    A_cell = (w / nw) * (h / nh)
    J_ta = np.abs(I_ta.reshape(nw, nh)) / A_cell
    J_pe = np.abs(I_pe.reshape(nw, nh)) / A_cell

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    extent = [-w / 2 * 1e3, w / 2 * 1e3, -h / 2 * 1e3, h / 2 * 1e3]

    vmax = max(J_ta.max(), J_pe.max())
    vmin = min(J_ta.min(), J_pe.min())

    im0 = axes[0].imshow(J_ta.T, origin='lower', extent=extent,
                          aspect='equal', cmap='viridis', vmin=vmin, vmax=vmax)
    axes[0].set_title('Tier A  |J|')
    axes[0].set_xlabel('u [mm]')
    axes[0].set_ylabel('v [mm]')
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(J_pe.T, origin='lower', extent=extent,
                          aspect='equal', cmap='viridis', vmin=vmin, vmax=vmax)
    axes[1].set_title('PEEC  |J|')
    axes[1].set_xlabel('u [mm]')
    fig.colorbar(im1, ax=axes[1])

    diff = np.abs(I_pe.reshape(nw, nh) - I_ta.reshape(nw, nh)) / A_cell
    im2 = axes[2].imshow(diff.T, origin='lower', extent=extent,
                          aspect='equal', cmap='inferno')
    axes[2].set_title('|PEEC - TierA|')
    axes[2].set_xlabel('u [mm]')
    fig.colorbar(im2, ax=axes[2])

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    sigma = 5.8e7

    # DC: pure path-length bias. PEEC should agree exactly since M has no
    # effect at omega=0 and the problem reduces to a resistor network.
    I_ta_dc, I_pe_dc, info = compare(nw=5, nh=5, frequency=0.0, sigma=sigma,
                                      n_arc=10, label='DC')
    plot_cross_section(I_ta_dc, I_pe_dc, 5, 5, info,
                       os.path.join(HERE, 'out_peec_vs_tiera_dc.png'),
                       'DC: Tier A vs PEEC-bundle')

    # 100 Hz: delta >> w_coil (delta ~ 6.6 mm vs 5 mm coil). Skin effect is
    # weak; both methods should collapse onto the DC distribution.
    I_ta, I_pe, info = compare(nw=5, nh=5, frequency=100.0, sigma=sigma,
                                n_arc=10, label='100 Hz (weak skin)')
    plot_cross_section(I_ta, I_pe, 5, 5, info,
                       os.path.join(HERE, 'out_peec_vs_tiera_100Hz.png'),
                       '100 Hz: Tier A vs PEEC-bundle')

    # AC: 7 kHz (IH band). Skin effect kicks in (delta ~ 0.79 mm, cell ~ 1 mm).
    # Tier A's "a priori" exp(-d/delta) factor vs PEEC's true self-consistent
    # solve on the same discretization.
    I_ta_ac, I_pe_ac, info = compare(nw=5, nh=5, frequency=7e3, sigma=sigma,
                                      n_arc=10, label='AC 7 kHz')
    plot_cross_section(I_ta_ac, I_pe_ac, 5, 5, info,
                       os.path.join(HERE, 'out_peec_vs_tiera_ac_7kHz.png'),
                       'AC 7 kHz: Tier A vs PEEC-bundle')

    # 7 kHz, higher resolution: check mesh sensitivity. delta/cell shrinks.
    I_ta9, I_pe9, info9 = compare(nw=9, nh=9, frequency=7e3, sigma=sigma,
                                    n_arc=10, label='AC 7 kHz nw=nh=9')
    plot_cross_section(I_ta9, I_pe9, 9, 9, info9,
                       os.path.join(HERE, 'out_peec_vs_tiera_ac_7kHz_n9.png'),
                       'AC 7 kHz (n=9): Tier A vs PEEC-bundle')

    print("\nDone.")


if __name__ == '__main__':
    main()
