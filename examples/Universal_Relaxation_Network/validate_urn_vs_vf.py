#!/usr/bin/env python
"""
URN vs Vector Fitting Validation Script

This script validates the Universal Relaxation Network (URN) against
Vector Fitting on real-world measurement data:
1. Li-ion battery EIS data
2. MnZn ferrite complex permeability

Metrics compared:
- Fitting accuracy (NRMSE)
- Model complexity (number of parameters)
- Physical interpretability
- SPICE netlist generation

Requirements:
    pip install numpy scipy torch matplotlib

Usage:
    python validate_urn_vs_vf.py

Output:
    - Fitting comparison plots (PNG)
    - SPICE netlists
    - Validation results JSON
"""

import sys
import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)


def load_battery_data():
    """Load NASA Li-ion battery EIS data (real measurement)."""
    # Try real-world NASA data first
    data_path = Path(__file__).parent / 'data' / 'real_world' / 'nasa_battery' / 'nasa_18650_eis.csv'
    if not data_path.exists():
        # Fall back to synthetic data
        data_path = Path(__file__).parent / 'data' / 'synthetic' / 'liion_battery_eis.csv'
        print(f"  Note: Using synthetic data ({data_path.name})")
        skiprows = 7  # Synthetic data has 7 header lines
    else:
        print(f"  Using real NASA 18650 battery EIS data")
        skiprows = 24  # NASA data has 24 header lines
    data = np.loadtxt(data_path, delimiter=',', skiprows=skiprows)
    freq = data[:, 0]
    Z_real = data[:, 1]
    Z_imag = data[:, 2]
    return freq, Z_real + 1j * Z_imag


def load_ferrite_data():
    """Load TDK MnZn ferrite impedance data (real measurement)."""
    # Try real-world TDK data first
    data_path = Path(__file__).parent / 'data' / 'real_world' / 'tdk_ferrite' / 'tdk_pc50_impedance.csv'
    if not data_path.exists():
        # Fall back to synthetic data
        data_path = Path(__file__).parent / 'data' / 'synthetic' / 'mnzn_ferrite_impedance.csv'
        print(f"  Note: Using synthetic data ({data_path.name})")
        skiprows = 7  # Synthetic data has 7 header lines
    else:
        print(f"  Using real TDK PC50 ferrite impedance data")
        skiprows = 19  # TDK data has 19 header lines
    data = np.loadtxt(data_path, delimiter=',', skiprows=skiprows)
    freq = data[:, 0]
    Z_real = data[:, 1]  # Z_real is column 1
    Z_imag = data[:, 2]  # Z_imag is column 2
    return freq, Z_real + 1j * Z_imag


def compute_nrmse(Z_true, Z_pred):
    """Compute Normalized Root Mean Square Error."""
    error = np.abs(Z_true - Z_pred)
    nrmse = np.sqrt(np.mean(error**2)) / np.mean(np.abs(Z_true))
    return nrmse


def vector_fitting(freq, Z_data, n_poles=10, n_iter=5):
    """
    Simple Vector Fitting implementation.

    Uses rational function approximation with real and complex conjugate poles.
    """
    omega = 2 * np.pi * freq
    s = 1j * omega

    # Initialize poles (logarithmically spaced)
    f_min, f_max = freq.min(), freq.max()
    pole_freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_poles // 2)

    # Real poles (for low-frequency behavior)
    poles_real = -2 * np.pi * pole_freqs[:n_poles//4]

    # Complex conjugate poles (for resonances)
    poles_complex = []
    for f in pole_freqs[n_poles//4:]:
        Q = 2.0  # Quality factor
        omega_p = 2 * np.pi * f
        sigma_p = -omega_p / (2 * Q)
        omega_d = omega_p * np.sqrt(1 - 1/(4*Q**2))
        poles_complex.append(sigma_p + 1j * omega_d)
        poles_complex.append(sigma_p - 1j * omega_d)

    poles = np.concatenate([poles_real, np.array(poles_complex)])

    # Iterative pole relocation
    for iteration in range(n_iter):
        # Build basis matrix
        n_freq = len(freq)
        n_p = len(poles)

        # Residue matrix
        A = np.zeros((2 * n_freq, 2 * n_p + 2), dtype=float)
        b = np.zeros(2 * n_freq, dtype=float)

        for i, (si, Zi) in enumerate(zip(s, Z_data)):
            basis = 1.0 / (si - poles)

            # Real part equation
            A[2*i, :n_p] = basis.real
            A[2*i, n_p:2*n_p] = -basis.imag
            A[2*i, -2] = 1.0  # d (constant term)
            A[2*i, -1] = si.imag  # e (linear term)
            b[2*i] = Zi.real

            # Imaginary part equation
            A[2*i+1, :n_p] = basis.imag
            A[2*i+1, n_p:2*n_p] = basis.real
            A[2*i+1, -2] = 0.0
            A[2*i+1, -1] = -si.real
            b[2*i+1] = Zi.imag

        # Solve least squares
        try:
            x, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
            residues = x[:n_p] + 1j * x[n_p:2*n_p]
            d = x[-2]
            e = x[-1]
        except np.linalg.LinAlgError:
            break

    # Compute fitted response
    Z_fit = np.zeros_like(s, dtype=complex)
    for i, si in enumerate(s):
        Z_fit[i] = d + e * si
        for r, p in zip(residues, poles):
            Z_fit[i] += r / (si - p)

    return Z_fit, poles, residues, d, e


def run_urn_fitting(freq, Z_data, config=None):
    """Run URN fitting on impedance data."""
    if config is None:
        config = URNConfig(
            n_debye=4,
            n_cole_cole=2,
            n_cole_davidson=1,
            n_havriliak_negami=1,
            n_cpe=2,
            n_warburg=2,
            n_gerischer=1,
            n_rlc=1,
            n_skin_effect=1,
            n_debye_parallel=2,
            n_cole_cole_parallel=1,
            n_cpe_parallel=1,
            n_warburg_parallel=1,
            sparsity_weight=0.01,
            lr=0.02,
            n_epochs=4000,
            n_restarts=3
        )

    start_time = time.time()
    model = train_urn(freq, Z_data, config)  # train_urn returns only model
    train_time = time.time() - start_time

    # Get fitted response
    import torch
    omega = 2 * np.pi * freq
    omega_tensor = torch.tensor(omega, dtype=torch.float64)
    Z_fit = model(omega_tensor).detach().numpy()

    # Count active parameters via get_active_components
    active_components = model.get_active_components()
    active_params = sum(len(v) for v in active_components.values())

    # Get discovered mechanisms (renamed API)
    mechanisms = active_components

    return Z_fit, model, active_params, mechanisms, train_time


def generate_comparison_plot(freq, Z_data, Z_urn, Z_vf, title, filename):
    """Generate comparison plot."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Nyquist plot
    ax1 = axes[0, 0]
    ax1.plot(Z_data.real, -Z_data.imag, 'ko', markersize=4, label='Measured')
    ax1.plot(Z_urn.real, -Z_urn.imag, 'b-', linewidth=2, label='URN')
    ax1.plot(Z_vf.real, -Z_vf.imag, 'r--', linewidth=2, label='VF')
    ax1.set_xlabel('Re(Z) [Ohm]')
    ax1.set_ylabel('-Im(Z) [Ohm]')
    ax1.set_title('Nyquist Plot')
    ax1.legend()
    ax1.set_aspect('equal')
    ax1.grid(True)

    # Bode magnitude
    ax2 = axes[0, 1]
    ax2.loglog(freq, np.abs(Z_data), 'ko', markersize=4, label='Measured')
    ax2.loglog(freq, np.abs(Z_urn), 'b-', linewidth=2, label='URN')
    ax2.loglog(freq, np.abs(Z_vf), 'r--', linewidth=2, label='VF')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('|Z| [Ohm]')
    ax2.set_title('Bode Magnitude')
    ax2.legend()
    ax2.grid(True)

    # Bode phase
    ax3 = axes[1, 0]
    ax3.semilogx(freq, np.angle(Z_data, deg=True), 'ko', markersize=4, label='Measured')
    ax3.semilogx(freq, np.angle(Z_urn, deg=True), 'b-', linewidth=2, label='URN')
    ax3.semilogx(freq, np.angle(Z_vf, deg=True), 'r--', linewidth=2, label='VF')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Phase [deg]')
    ax3.set_title('Bode Phase')
    ax3.legend()
    ax3.grid(True)

    # Error comparison
    ax4 = axes[1, 1]
    error_urn = np.abs(Z_data - Z_urn) / np.abs(Z_data) * 100
    error_vf = np.abs(Z_data - Z_vf) / np.abs(Z_data) * 100
    ax4.semilogx(freq, error_urn, 'b-', linewidth=2, label=f'URN (mean={np.mean(error_urn):.2f}%)')
    ax4.semilogx(freq, error_vf, 'r--', linewidth=2, label=f'VF (mean={np.mean(error_vf):.2f}%)')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Relative Error [%]')
    ax4.set_title('Fitting Error')
    ax4.legend()
    ax4.grid(True)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved plot: {filename}")


def main():
    """Run validation suite."""
    print("=" * 70)
    print("URN vs Vector Fitting Validation Suite")
    print("=" * 70)

    results = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'datasets': {}
    }

    output_dir = Path(__file__).parent / 'validation_output'
    output_dir.mkdir(exist_ok=True)

    # =========================================================================
    # Dataset 1: Li-ion Battery EIS
    # =========================================================================
    print("\n[1/2] Li-ion Battery EIS Data")
    print("-" * 40)

    freq_bat, Z_bat = load_battery_data()
    print(f"  Loaded {len(freq_bat)} frequency points")
    print(f"  Frequency range: {freq_bat.min():.4f} Hz - {freq_bat.max():.0f} Hz")

    # Vector Fitting
    print("  Running Vector Fitting (10 poles)...")
    t0 = time.time()
    Z_vf_bat, poles, residues, d, e = vector_fitting(freq_bat, Z_bat, n_poles=10)
    vf_time = time.time() - t0
    nrmse_vf_bat = compute_nrmse(Z_bat, Z_vf_bat)
    print(f"    NRMSE: {nrmse_vf_bat:.4f}, Time: {vf_time:.2f}s")

    # URN Fitting
    print("  Running URN fitting...")
    config_bat = URNConfig(
        n_debye=3, n_cole_cole=2, n_cpe=2, n_warburg=2,
        n_debye_parallel=2, n_warburg_parallel=1,
        sparsity_weight=0.01, n_epochs=3000, n_restarts=3
    )
    Z_urn_bat, model_bat, n_params_bat, mechanisms_bat, urn_time = run_urn_fitting(
        freq_bat, Z_bat, config_bat
    )
    nrmse_urn_bat = compute_nrmse(Z_bat, Z_urn_bat)
    print(f"    NRMSE: {nrmse_urn_bat:.4f}, Time: {urn_time:.2f}s")
    print(f"    Active parameters: {n_params_bat}")
    print(f"    Discovered mechanisms: {list(mechanisms_bat.keys())}")

    # Generate comparison plot
    generate_comparison_plot(
        freq_bat, Z_bat, Z_urn_bat, Z_vf_bat,
        "Li-ion Battery EIS: URN vs Vector Fitting",
        output_dir / 'battery_eis_comparison.png'
    )

    # Generate SPICE netlist
    spice_bat = generate_spice_netlist(model_bat, "BATTERY_URN")
    spice_path_bat = output_dir / 'battery_urn.sp'
    with open(spice_path_bat, 'w') as f:
        f.write(spice_bat)
    print(f"  Saved SPICE netlist: {spice_path_bat}")

    results['datasets']['battery_eis'] = {
        'n_points': len(freq_bat),
        'freq_range': [float(freq_bat.min()), float(freq_bat.max())],
        'vf': {
            'nrmse': float(nrmse_vf_bat),
            'n_poles': 10,
            'time_s': vf_time
        },
        'urn': {
            'nrmse': float(nrmse_urn_bat),
            'n_active_params': n_params_bat,
            'mechanisms': list(mechanisms_bat.keys()),
            'time_s': urn_time
        }
    }

    # =========================================================================
    # Dataset 2: MnZn Ferrite
    # =========================================================================
    print("\n[2/2] MnZn Ferrite Impedance Data")
    print("-" * 40)

    freq_fer, Z_fer = load_ferrite_data()
    print(f"  Loaded {len(freq_fer)} frequency points")
    print(f"  Frequency range: {freq_fer.min():.0f} Hz - {freq_fer.max()/1e6:.1f} MHz")

    # Vector Fitting
    print("  Running Vector Fitting (12 poles)...")
    t0 = time.time()
    Z_vf_fer, poles, residues, d, e = vector_fitting(freq_fer, Z_fer, n_poles=12)
    vf_time = time.time() - t0
    nrmse_vf_fer = compute_nrmse(Z_fer, Z_vf_fer)
    print(f"    NRMSE: {nrmse_vf_fer:.4f}, Time: {vf_time:.2f}s")

    # URN Fitting
    print("  Running URN fitting...")
    config_fer = URNConfig(
        n_debye=4, n_cole_cole=3, n_cpe=1, n_warburg=1,
        n_skin_effect=2, n_rlc=2,
        n_debye_parallel=2, n_cole_cole_parallel=1,
        sparsity_weight=0.01, n_epochs=4000, n_restarts=3
    )
    Z_urn_fer, model_fer, n_params_fer, mechanisms_fer, urn_time = run_urn_fitting(
        freq_fer, Z_fer, config_fer
    )
    nrmse_urn_fer = compute_nrmse(Z_fer, Z_urn_fer)
    print(f"    NRMSE: {nrmse_urn_fer:.4f}, Time: {urn_time:.2f}s")
    print(f"    Active parameters: {n_params_fer}")
    print(f"    Discovered mechanisms: {list(mechanisms_fer.keys())}")

    # Generate comparison plot
    generate_comparison_plot(
        freq_fer, Z_fer, Z_urn_fer, Z_vf_fer,
        "MnZn Ferrite: URN vs Vector Fitting",
        output_dir / 'ferrite_impedance_comparison.png'
    )

    # Generate SPICE netlist
    spice_fer = generate_spice_netlist(model_fer, "FERRITE_URN")
    spice_path_fer = output_dir / 'ferrite_urn.sp'
    with open(spice_path_fer, 'w') as f:
        f.write(spice_fer)
    print(f"  Saved SPICE netlist: {spice_path_fer}")

    results['datasets']['ferrite_impedance'] = {
        'n_points': len(freq_fer),
        'freq_range': [float(freq_fer.min()), float(freq_fer.max())],
        'vf': {
            'nrmse': float(nrmse_vf_fer),
            'n_poles': 12,
            'time_s': vf_time
        },
        'urn': {
            'nrmse': float(nrmse_urn_fer),
            'n_active_params': n_params_fer,
            'mechanisms': list(mechanisms_fer.keys()),
            'time_s': urn_time
        }
    }

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"\n{'Dataset':<25} {'Method':<10} {'NRMSE':<10} {'Params':<10} {'Time':<10}")
    print("-" * 65)

    for name, data in results['datasets'].items():
        print(f"{name:<25} {'VF':<10} {data['vf']['nrmse']:.4f}    "
              f"{data['vf']['n_poles']:<10} {data['vf']['time_s']:.2f}s")
        print(f"{'':<25} {'URN':<10} {data['urn']['nrmse']:.4f}    "
              f"{data['urn']['n_active_params']:<10} {data['urn']['time_s']:.2f}s")
        print("-" * 65)

    # Save results
    results_path = output_dir / 'validation_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    print("\nValidation complete!")
    return results


if __name__ == '__main__':
    main()
