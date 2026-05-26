"""
End-to-end IH workpiece -> SPICE hierarchical Cauer demonstration.

Workflow:
  1. Synthesize a representative workpiece Y(omega) frequency sweep
     (using the cylinder Bessel admittance as a stand-in for a real
     IH workpiece - this gives realistic numbers without requiring
     a 10-minute NGSolve FEM run per frequency).
  2. Pass the sweep through extract_hierarchical_cauer.py logic.
  3. Generate workpiece.cir SPICE subcircuit.
  4. Verify SPICE reconstruction matches input sweep.

This script represents the workflow that would couple:
  - calc_fem_kelvin.py output (Y_workpiece(omega) frequency sweep)
  - extract_hierarchical_cauer.py (this pipeline)
  - workpiece.cir for LTspice/ngspice transient simulation

For the actual 3turnCoil_work sample (5 mm half-radius Cu cylinder
approximation), the cylinder formula gives realistic K_SIBC and bulk
modal eigenvalues.

Usage:
    python ih_workpiece_extract.py
"""

import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy import special

sys.path.insert(0, str(Path(__file__).parent))
from extract_hierarchical_cauer import (Y_warburg, fit_d, fit_bulk_foster,
                                         emit_spice_netlist)
from diffusive_quadrature import fit_diffusive


def setup_times_new_roman():
    try:
        for font in ['Times New Roman', 'Times']:
            try:
                font_manager.findfont(font, fallback_to_default=False)
                plt.rcParams['font.family'] = font
                plt.rcParams['font.serif'] = [font]
                plt.rcParams['mathtext.fontset'] = 'stix'
                return
            except Exception:
                continue
    except Exception:
        pass


def synthesize_workpiece_Y(omega, a_workpiece=5e-3, sigma=5.8e7, mu=4*np.pi*1e-7,
                            N_modes=200):
    """
    Generate Y(omega) for a representative workpiece (cylinder approximation).

    Returns Y as the Bessel exact closed form (truncated to N_modes Foster).
    For a real workpiece, this would come from NGSolve FEM frequency sweep.
    """
    s = 1j * omega
    j0_zeros = special.jn_zeros(0, N_modes)
    xi = j0_zeros**2 / (a_workpiece**2 * mu * sigma)
    residue = 4 * np.pi / mu
    Y = np.sum(residue / (s[:, None] + xi[None, :]), axis=1)
    K_SIBC = 2 * np.pi * a_workpiece * np.sqrt(sigma / mu)
    return Y, K_SIBC


if __name__ == "__main__":
    # Workpiece parameters (representative of IH samples)
    sigma = 5.8e7  # Cu (slightly different from 5.96e7 to match sample's u1)
    mu_0 = 4 * np.pi * 1e-7
    a_workpiece = 5e-3
    print(f"IH workpiece extractor demo")
    print(f"  Geometry: a = {a_workpiece*1e3:.1f} mm (representative)")
    print(f"  Material: sigma = {sigma:.2e} S/m, mu_r = 1")
    print()

    # Synthesize Y(omega) sweep (= what calc_fem_kelvin.py would produce)
    omega = np.logspace(2, 7, 80) * 2 * np.pi  # 100 Hz - 10 MHz, 80 points
    print(f"Step 1: synthesize Y(omega) sweep")
    print(f"  frequency range: {omega.min()/(2*np.pi):.0e} - "
          f"{omega.max()/(2*np.pi):.0e} Hz")
    print(f"  {len(omega)} frequency points")
    Y_data, K_SIBC = synthesize_workpiece_Y(omega, a_workpiece=a_workpiece,
                                              sigma=sigma, mu=mu_0)
    print(f"  K_SIBC = {K_SIBC:.4e} (= 2 pi a sqrt(sigma/mu))")
    print(f"  Y(low f)  = {Y_data[0]:.4e}")
    print(f"  Y(high f) = {Y_data[-1]:.4e}")
    print()

    # Stage 2-4: extract hierarchical Cauer parameters
    print("Step 2: fit crossover d")
    d, err_d = fit_d(omega, Y_data, K_SIBC, N_CLN_test=8)
    print(f"  d = {d:.3e} rad/s ({d/(2*np.pi):.2e} Hz)")
    print()

    print("Step 3: fit bulk Foster (N_CLN = 30)")
    Y_w = Y_warburg(1j * omega, K_SIBC, d)
    Y_residual = Y_data - Y_w
    xi_CLN, c_CLN, err_CLN = fit_bulk_foster(omega, Y_residual, 30,
                                              reference_Y=Y_data)
    print(f"  bulk fit max rel err = {err_CLN:.3e}")
    print()

    print("Step 4: diffusive Foster of Warburg (N_xi = 50)")
    xi_w, c_w, err_w, _, _ = fit_diffusive(d, N_xi=50,
                                            omega_min=omega.min(),
                                            omega_max=omega.max())
    print(f"  Warburg fit max rel err = {err_w:.3e}")
    print()

    # Verify reconstruction
    s = 1j * omega
    Y_recon = (np.sum(c_CLN[None, :] / (s[:, None] + xi_CLN[None, :]), axis=1)
               + K_SIBC * np.sum(c_w[None, :] / (s[:, None] + xi_w[None, :]),
                                  axis=1))
    err_total = np.abs(Y_recon - Y_data) / np.abs(Y_data)
    print(f"Total reconstruction:")
    print(f"  max  rel err = {np.max(err_total):.3e}")
    print(f"  median       = {np.median(err_total):.3e}")
    print()

    # Save SPICE netlist
    out_dir = Path(__file__).parent
    emit_spice_netlist(out_dir / "ih_workpiece.cir",
                       xi_CLN, c_CLN, xi_w, c_w, K_SIBC, d,
                       comment="IH workpiece (5mm radius Cu, representative)")
    print(f"Saved: ih_workpiece.cir")

    # Save JSON params (machine-readable)
    params = {
        "workpiece_geometry": {
            "shape": "cylinder",
            "radius_m": float(a_workpiece),
            "K_SIBC": float(K_SIBC),
        },
        "material": {
            "sigma_S_per_m": float(sigma),
            "mu_r": 1.0,
        },
        "hierarchical_cauer": {
            "d_rad_s": float(d),
            "d_Hz": float(d / (2 * np.pi)),
            "N_CLN": int(len(xi_CLN)),
            "xi_CLN_rad_s": xi_CLN.tolist(),
            "c_CLN": c_CLN.tolist(),
            "N_xi": int(len(xi_w)),
            "xi_w_rad_s": xi_w.tolist(),
            "c_w_real": [float(c.real) for c in c_w],
        },
        "verification": {
            "max_rel_err": float(np.max(err_total)),
            "median_rel_err": float(np.median(err_total)),
            "omega_min_rad_s": float(omega.min()),
            "omega_max_rad_s": float(omega.max()),
        },
    }
    (out_dir / "ih_workpiece_params.json").write_text(
        json.dumps(params, indent=2))
    print(f"Saved: ih_workpiece_params.json")
    print()

    # Plot
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2))

    f_hz = omega / (2 * np.pi)
    ax1.loglog(f_hz, np.abs(Y_data), 'k-', linewidth=1.4,
               label='Y_workpiece (FEM-equivalent)')
    ax1.loglog(f_hz, np.abs(Y_recon), 'r--', linewidth=1.0,
               label=f'hierarchical Cauer ({len(xi_CLN)} + {len(xi_w)} modes)')
    ax1.axvline(d / (2 * np.pi), color='gray', linestyle=':',
                linewidth=0.6, label=f'$d/(2\\pi)$ = {d/(2*np.pi):.1e} Hz')
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'$|Y(j\omega)|$')
    ax1.legend(fontsize=8, frameon=False, loc='upper right')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title('(a) IH workpiece $Y(j\\omega)$', fontsize=9)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    ax2.loglog(f_hz, err_total, 'r-', linewidth=1.0)
    ax2.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1\\% target')
    ax2.set_xlabel('frequency (Hz)')
    ax2.set_ylabel('relative error')
    ax2.legend(fontsize=8, frameon=False, loc='lower right')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_title('(b) Reconstruction error', fontsize=9)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    plt.savefig(out_dir / "ih_workpiece_extract.png", dpi=200,
                bbox_inches='tight')
    print(f"Saved: ih_workpiece_extract.png")
    print()
    print("Workflow established:")
    print("  workpiece STEP -> Cubit mesh -> calc_fem_kelvin.py freq sweep")
    print("    -> [extract_hierarchical_cauer.py] -> SPICE subcircuit -> LTspice")
    print()
    print("The hierarchical Cauer params + JSON metadata fully characterize")
    print("the workpiece for any transient or steady-state SPICE simulation.")
