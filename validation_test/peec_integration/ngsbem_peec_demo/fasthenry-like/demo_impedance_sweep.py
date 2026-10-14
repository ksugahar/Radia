"""
FastHenry-like PEEC: Impedance sweep on a rectangular plate.

Models a 10x10mm copper plate (35 um thick) as a 4-segment rectangular
loop using Radia's FastHenry parser and PEEC filament solver.

Requires: radia (peec_matrices.pyd), fasthenry_parser, peec_topology
"""
import sys
import os
import time
import numpy as np

# Path setup
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_script_dir, '..'))
sys.path.insert(0, os.path.join(_script_dir, '../../../../src/radia'))

# MKL DLL directory for peec_matrices.pyd
_mkl_bin = os.path.join(sys.prefix, "Library", "bin")
if os.path.isdir(_mkl_bin):
    os.add_dll_directory(_mkl_bin)

MU_0 = 4.0 * np.pi * 1e-7

from fasthenry_parser import FastHenryParser
from peec_topology import PEECCircuitSolver


def dowell_F_R(xi):
    """Dowell skin effect factor F_R(xi), xi = d / (2 * delta)."""
    if xi < 1e-6:
        return 1.0
    s2 = np.sinh(2 * xi)
    c2 = np.cosh(2 * xi)
    sn = np.sin(2 * xi)
    cs = np.cos(2 * xi)
    return xi * (s2 + sn) / (c2 - cs)


def main():
    # --- Geometry (same as ngbem/demo_impedance_sweep.py) ---
    width = 0.01        # 10 mm
    height = 0.01       # 10 mm
    thickness = 35e-6   # 35 um copper
    sigma = 5.8e7       # Cu conductivity [S/m]

    w_mm = width * 1e3
    h_mm = height * 1e3
    t_mm = thickness * 1e3

    print("=== FastHenry PEEC: Flat-plate impedance sweep ===")
    print(f"Conductor: {w_mm:.0f} x {h_mm:.0f} mm, t = {thickness*1e6:.0f} um Cu")

    # --- FastHenry model: 4-segment rectangular loop ---
    inp_text = f""".Units mm
.default sigma={sigma:.2e}

N1 x=0       y=0       z=0
N2 x={w_mm}  y=0       z=0
N3 x={w_mm}  y={h_mm}  z=0
N4 x=0       y={h_mm}  z=0

E1 N1 N2 w={w_mm} h={t_mm}
E2 N2 N3 w={h_mm} h={t_mm}
E3 N3 N4 w={w_mm} h={t_mm}
E4 N4 N1 w={h_mm} h={t_mm}

.external N1 N1
.freq fmin=10 fmax=1e6 ndec=10
.end"""

    t0 = time.perf_counter()
    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    t_build = time.perf_counter() - t0

    n_seg = topo['n_loop']
    R_dc = np.array(topo['R'])
    L = np.array(topo['L'])

    print(f"DOFs: {n_seg} segments (filaments)")
    print(f"Build time: {t_build*1e3:.1f} ms")
    print(f"R_dc (per seg): [{', '.join(f'{r*1e3:.4f}' for r in R_dc)}] mOhm")
    print(f"L diagonal [nH]: [{', '.join(f'{L[i,i]*1e9:.2f}' for i in range(n_seg))}]")

    solver = PEECCircuitSolver(topo)

    # --- Frequency sweep ---
    freqs = np.logspace(1, 6, 50)  # 10 Hz to 1 MHz

    # DC resistance only
    t0 = time.perf_counter()
    Z_dc = solver.frequency_sweep(freqs)
    t_dc = time.perf_counter() - t0

    # With Dowell skin effect
    def dowell_Zs(freq):
        if freq <= 0:
            return np.zeros(n_seg, dtype=complex)
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        xi = thickness / (2.0 * delta)
        F = dowell_F_R(xi)
        return R_dc * (F - 1.0) + 0j

    t0 = time.perf_counter()
    Z_dowell = solver.frequency_sweep(freqs, Zs_func=dowell_Zs)
    t_dowell = time.perf_counter() - t0

    print(f"\nSweep times ({len(freqs)} freqs):")
    print(f"  DC:     {t_dc*1e3:.1f} ms")
    print(f"  Dowell: {t_dowell*1e3:.1f} ms")

    # --- Results table ---
    print(f"\n{'Freq [Hz]':>12s}  {'L_dc [nH]':>12s}  {'L_dw [nH]':>12s}  "
          f"{'R_dc [mOhm]':>13s}  {'R_dw [mOhm]':>13s}")
    print("-" * 70)
    for f, Zd, Zdw in list(zip(freqs, Z_dc, Z_dowell))[::10]:
        omega = 2 * np.pi * f
        Ld = np.imag(Zd) / omega * 1e9
        Ldw = np.imag(Zdw) / omega * 1e9
        Rd = np.real(Zd) * 1e3
        Rdw = np.real(Zdw) * 1e3
        print(f"{f:12.1e}  {Ld:12.2f}  {Ldw:12.2f}  {Rd:13.4f}  {Rdw:13.4f}")

    # --- Physics check ---
    omega_lo = 2 * np.pi * freqs[0]
    omega_hi = 2 * np.pi * freqs[-1]
    L_lo = np.imag(Z_dc[0]) / omega_lo * 1e9
    L_hi = np.imag(Z_dc[-1]) / omega_hi * 1e9
    R_dc_lo = np.real(Z_dc[0]) * 1e3
    R_dw_hi = np.real(Z_dowell[-1]) * 1e3
    print(f"\nPhysics check:")
    print(f"  L constant (DC):    L_10Hz = {L_lo:.2f} nH,  L_1MHz = {L_hi:.2f} nH")
    print(f"  R_ac > R_dc (1MHz): R_dc = {R_dc_lo:.4f} mOhm,  R_dowell = {R_dw_hi:.4f} mOhm")

    # --- Plot ---
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

        for Z, label, ls in [(Z_dc, 'DC only', '--'),
                              (Z_dowell, 'Dowell', '-')]:
            axes[0].loglog(freqs, np.abs(Z), ls, linewidth=1.5, label=label)
        axes[0].set_ylabel('|Z| [Ohm]')
        axes[0].set_title('FastHenry PEEC: 10x10 mm Cu plate impedance')
        axes[0].legend()
        axes[0].grid(True, which='both', alpha=0.3)

        for Z, label, ls in [(Z_dc, 'DC only', '--'),
                              (Z_dowell, 'Dowell', '-')]:
            axes[1].semilogx(freqs, np.angle(Z, deg=True), ls, linewidth=1.5, label=label)
        axes[1].set_xlabel('Frequency [Hz]')
        axes[1].set_ylabel('Phase [deg]')
        axes[1].legend()
        axes[1].grid(True, which='both', alpha=0.3)

        plt.tight_layout()
        out_path = os.path.join(_script_dir, 'demo_impedance_sweep.png')
        plt.savefig(out_path, dpi=150)
        print(f"\nFigure saved: {out_path}")
        plt.close()
    except ImportError:
        print("\nmatplotlib not available, skipping plot.")


if __name__ == '__main__':
    main()
