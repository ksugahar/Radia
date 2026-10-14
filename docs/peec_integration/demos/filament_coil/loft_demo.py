"""loft_demo.py

Demonstrate the Stage 3a loft capability end-to-end:

  (1) Build a coil with a variable cross-section using fluent API
      (two loft segments + a straight mid-section) so the shape
      actually morphs in 3D.
  (2) Export the CAD solid as a STEP file via OCC ThruSections loft.
  (3) Plot the filament paths in 3D matplotlib to visualize how
      individual filaments converge / diverge through the loft.
  (4) Run PEEC on the loft coil at 7 kHz and show the current
      magnitude map on the narrow cross-section.

The coil used here is a "bow-tie" straight: 5x5 -> 2x5 -> 5x5 over
120 mm. It is not a real IH geometry; the point is to exercise the
loft plumbing on a geometry where the taper is visually obvious.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- registers 3d proj

HERE = os.path.dirname(os.path.abspath(__file__))

from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile
from radia.peec_bundle import build_bundle_solver, filament_currents

MM = 1e-3
CU_SIGMA = 5.8e7


def build_bowtie():
    """5x5 -> 2x5 (40mm taper) -> 2x5 (40mm straight) -> 5x5 (40mm taper).

    Total length 120 mm along Y axis.
    """
    return (CoilBuilder(current=100.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(width=5 * MM, height=5 * MM)
            .add_loft_straight(profile_end=RectProfile(2 * MM, 5 * MM),
                               length=40 * MM, n_sub=20)
            .add_straight(length=40 * MM)
            .add_loft_straight(profile_end=RectProfile(5 * MM, 5 * MM),
                               length=40 * MM, n_sub=20))


def step_out_step(coil, path):
    """Export CAD solid via OCC loft. Returns path on success."""
    try:
        shape = coil.to_occ()
        shape.WriteStep(path)
        print(f"  STEP saved: {path}")
        return True
    except Exception as e:
        print(f"  STEP export FAILED: {type(e).__name__}: {e}")
        return False


def plot_filaments_3d(coil, nw, nh, out_path):
    paths, _ = coil.to_filaments(nw=nw, nh=nh, frequency=0.0,
                                  sigma=CU_SIGMA)
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    N = len(paths)

    # Color by filament (corner vs center)
    nw_ = nw
    nh_ = nh
    for k, path in enumerate(paths):
        i = k // nh_
        j = k % nh_
        # corner emphasis
        is_corner = (i in (0, nw_ - 1)) and (j in (0, nh_ - 1))
        is_center = (i == nw_ // 2) and (j == nh_ // 2)
        if is_corner:
            color = 'tab:red'
            lw = 1.0
            alpha = 0.9
        elif is_center:
            color = 'tab:blue'
            lw = 1.5
            alpha = 1.0
        else:
            color = 'lightgray'
            lw = 0.4
            alpha = 0.5

        # Flatten polyline
        waypoints = [path[0][0]]
        for (_p1, p2) in path:
            waypoints.append(p2)
        wp = np.array(waypoints) * 1e3  # mm for plotting
        ax.plot(wp[:, 0], wp[:, 1], wp[:, 2],
                color=color, lw=lw, alpha=alpha)

    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    ax.set_title(f'Loft filaments (nw={nw}, nh={nh}, N={N})\n'
                  'red=corners, blue=center, gray=interior')
    ax.view_init(elev=22, azim=-60)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  3D filament plot: {out_path}")


def run_peec_and_plot(coil, nw, nh, freq, out_path):
    paths, _ = coil.to_filaments(nw=nw, nh=nh, frequency=freq,
                                  sigma=CU_SIGMA)
    info = coil._last_filament_info
    solver, seg_of_fil, _p, _m = build_bundle_solver(
        paths, dw=None, dh=None, sigma=CU_SIGMA,
        cell_wh=info['cell_wh'])
    I_branch = solver.compute_branch_currents(freq, [complex(coil.current)])
    I = filament_currents(I_branch, seg_of_fil)
    J = np.abs(I).reshape(nw, nh)
    phase = np.angle(I.reshape(nw, nh), deg=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    im0 = axes[0].imshow(J.T, origin='lower', cmap='viridis', aspect='equal')
    axes[0].set_title(f'|I_k|  bow-tie loft, f={freq/1e3:g} kHz')
    axes[0].set_xlabel('i (width idx)')
    axes[0].set_ylabel('j (height idx)')
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(phase.T, origin='lower', cmap='twilight',
                          vmin=-180, vmax=180, aspect='equal')
    axes[1].set_title('arg(I_k) [deg]')
    axes[1].set_xlabel('i')
    fig.colorbar(im1, ax=axes[1])
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)

    print(f"  |I_k| range [{J.min():.4f}, {J.max():.4f}]  "
          f"max/min = {J.max()/J.min():.2f}")
    print(f"  arg range [{phase.min():.1f}, {phase.max():.1f}] deg")
    print(f"  PEEC heatmap: {out_path}")
    print(f"  PEEC n_loop = {solver.n_loop}")


def main():
    coil = build_bowtie()
    print(f"Segments: {len(coil.segments)}")
    for i, s in enumerate(coil.segments):
        print(f"  [{i}] {type(s).__name__}  w={s.width*1e3:.2f} "
              f"h={s.height*1e3:.2f} mm  "
              f"len={getattr(s, 'length', 0)*1e3:.1f} mm")

    print("\n[1] STEP export")
    step_path = os.path.join(HERE, 'out_loft_bowtie.step')
    step_out_step(coil, step_path)

    print("\n[2] 3D filament plot")
    plot_filaments_3d(coil, nw=7, nh=7,
                      out_path=os.path.join(HERE, 'out_loft_bowtie_fils.png'))

    print("\n[3] PEEC at 7 kHz")
    run_peec_and_plot(coil, nw=5, nh=5, freq=7e3,
                      out_path=os.path.join(HERE, 'out_loft_bowtie_peec.png'))

    print("\nDone.")


if __name__ == '__main__':
    main()
