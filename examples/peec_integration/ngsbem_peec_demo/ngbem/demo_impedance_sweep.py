"""
ngbem PEEC: Impedance sweep on a flat-plate conductor.

Demonstrates:
  1. BEM matrix assembly (LaplaceSL on HDivSurface x SurfaceL2)
  2. MQS impedance sweep (loop-only: Z = R + jw*L)
  3. Full Loop-Star impedance (includes capacitive P block)
  4. Stabilized mode (Weggler block system, O(1) conditioning)

Requires: ngsolve>=6.2.2601, ngsolve-ngsbem, matplotlib
"""
import sys, os
import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from ngbem_peec import (NGBEMPEECSolver, create_plate_mesh,
                        make_dowell_Zs_func, MU_0, EPS_0)


def main():
    # --- Geometry ---
    width = 0.01       # 10 mm
    height = 0.01      # 10 mm
    maxh = 0.003       # ~3 mm element size
    thickness = 35e-6   # 35 um copper
    sigma = 5.8e7       # Cu conductivity [S/m]

    print("=== ngbem PEEC: Flat-plate impedance sweep ===")
    print(f"Conductor: {width*1e3:.0f} x {height*1e3:.0f} mm, "
          f"t = {thickness*1e6:.0f} um Cu")

    # --- Mesh and assembly ---
    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    print(f"Mesh: {mesh.ne} surface elements")

    t0 = time.perf_counter()
    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=sigma, thickness=thickness,
                              order=0, intorder=5)
    solver.assemble()
    t_assemble = time.perf_counter() - t0

    print(f"DOFs: {solver.n_loop} loop + {solver.n_star} star")
    print(f"Assembly time: {t_assemble*1e3:.1f} ms")

    # --- Frequency sweep: 3 modes ---
    freqs = np.logspace(1, 6, 50)  # 10 Hz to 1 MHz

    t0 = time.perf_counter()
    Z_mqs = solver.solve_frequency(freqs, mode='mqs')
    t_mqs = time.perf_counter() - t0

    t0 = time.perf_counter()
    Z_full = solver.solve_frequency(freqs, mode='full')
    t_full = time.perf_counter() - t0

    t0 = time.perf_counter()
    Z_stab = solver.solve_frequency(freqs, mode='stabilized')
    t_stab = time.perf_counter() - t0

    print(f"\nSweep times ({len(freqs)} freqs):")
    print(f"  mqs:        {t_mqs*1e3:.1f} ms")
    print(f"  full:       {t_full*1e3:.1f} ms")
    print(f"  stabilized: {t_stab*1e3:.1f} ms")

    # --- With Dowell skin effect ---
    Zs_dowell = make_dowell_Zs_func(
        n_loop=solver.n_loop,
        R_dc_per_edge=solver.R_loop,
        thickness=thickness,
        sigma=sigma,
    )
    Z_mqs_dowell = solver.solve_frequency(freqs, mode='mqs', Zs_func=Zs_dowell)

    # --- Results table ---
    print(f"\n{'Freq [Hz]':>12s}  {'L_mqs [nH]':>12s}  {'L_full [nH]':>12s}  "
          f"{'L_stab [nH]':>12s}  {'R_mqs [mOhm]':>13s}  {'R_dowell [mOhm]':>15s}")
    print("-" * 85)
    for f, Zm, Zf, Zs, Zd in list(zip(freqs, Z_mqs, Z_full, Z_stab, Z_mqs_dowell))[::10]:
        omega = 2 * np.pi * f
        Lm = np.imag(Zm) / omega * 1e9
        Lf = np.imag(Zf) / omega * 1e9
        Ls = np.imag(Zs) / omega * 1e9
        Rm = np.real(Zm) * 1e3
        Rd = np.real(Zd) * 1e3
        print(f"{f:12.1e}  {Lm:12.2f}  {Lf:12.2f}  {Ls:12.2f}  {Rm:13.4f}  {Rd:15.4f}")

    # --- Validation ---
    rel_diff_L = np.abs(np.imag(Z_full) - np.imag(Z_stab)) / np.abs(np.imag(Z_stab)) * 100
    print(f"\nfull vs stabilized: max L diff = {np.max(rel_diff_L):.2f}%")
    print(f"  => Both modes give same physics (Schur complement = block system)")

    # --- Plot ---
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

        for Z, label, ls in [(Z_mqs, 'mqs (DC)', '--'),
                              (Z_mqs_dowell, 'mqs (Dowell)', '-'),
                              (Z_full, 'full', ':')]:
            axes[0].loglog(freqs, np.abs(Z), ls, linewidth=1.5, label=label)
        axes[0].set_ylabel('|Z| [Ohm]')
        axes[0].set_title('ngbem PEEC: 10x10 mm Cu plate impedance')
        axes[0].legend()
        axes[0].grid(True, which='both', alpha=0.3)

        for Z, label, ls in [(Z_mqs, 'mqs (DC)', '--'),
                              (Z_mqs_dowell, 'mqs (Dowell)', '-'),
                              (Z_full, 'full', ':')]:
            axes[1].semilogx(freqs, np.angle(Z, deg=True), ls, linewidth=1.5, label=label)
        axes[1].set_xlabel('Frequency [Hz]')
        axes[1].set_ylabel('Phase [deg]')
        axes[1].legend()
        axes[1].grid(True, which='both', alpha=0.3)

        plt.tight_layout()
        out_path = os.path.join(os.path.dirname(__file__), 'demo_impedance_sweep.png')
        plt.savefig(out_path, dpi=150)
        print(f"\nFigure saved: {out_path}")
        plt.close()
    except ImportError:
        print("\nmatplotlib not available, skipping plot.")


if __name__ == '__main__':
    main()
