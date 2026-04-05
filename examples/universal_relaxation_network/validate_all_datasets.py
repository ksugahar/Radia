#!/usr/bin/env python
"""
Comprehensive URN vs Vector Fitting Validation on All Real-World Datasets

This script validates URN on all available real-world datasets:
1. NASA 18650 Battery EIS
2. TDK PC47 Ferrite
3. TDK PC50 Ferrite
4. TDK PC95 Ferrite
5. TDK PC200 Ferrite

Output:
- Comprehensive comparison table
- Individual plots for each dataset
- JSON results file

Usage:
    python validate_all_datasets.py
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
import torch


def load_data(data_path, skiprows):
    """Load impedance data from CSV file."""
    data = np.loadtxt(data_path, delimiter=',', skiprows=skiprows)
    freq = data[:, 0]
    Z_real = data[:, 1]
    Z_imag = data[:, 2]
    return freq, Z_real + 1j * Z_imag


def compute_nrmse(Z_true, Z_pred):
    """Compute Normalized Root Mean Square Error."""
    error = np.abs(Z_true - Z_pred)
    nrmse = np.sqrt(np.mean(error**2)) / np.mean(np.abs(Z_true))
    return nrmse


def vector_fitting(freq, Z_data, n_poles=10, n_iter=5):
    """Simple Vector Fitting implementation."""
    omega = 2 * np.pi * freq
    s = 1j * omega

    # Initialize poles (logarithmically spaced)
    f_min, f_max = freq.min(), freq.max()
    pole_freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_poles // 2)

    # Real poles
    poles_real = -2 * np.pi * pole_freqs[:n_poles//4]

    # Complex conjugate poles
    poles_complex = []
    for f in pole_freqs[n_poles//4:]:
        Q = 2.0
        omega_p = 2 * np.pi * f
        sigma_p = -omega_p / (2 * Q)
        omega_d = omega_p * np.sqrt(1 - 1/(4*Q**2))
        poles_complex.append(sigma_p + 1j * omega_d)
        poles_complex.append(sigma_p - 1j * omega_d)

    poles = np.concatenate([poles_real, np.array(poles_complex)])

    # Iterative pole relocation
    for iteration in range(n_iter):
        # Build basis matrix
        A = np.zeros((len(s), len(poles) + 1), dtype=complex)
        A[:, 0] = 1  # d term
        for i, p in enumerate(poles):
            A[:, i + 1] = 1 / (s - p)

        # Solve least squares
        coeffs, residuals, rank, sv = np.linalg.lstsq(A, Z_data, rcond=None)

        # Compute fitted response
        Z_fit = A @ coeffs

    return Z_fit


def run_urn_fitting(freq, Z_data, n_restarts=3):
    """Run URN fitting on impedance data."""
    config = URNConfig(
        n_debye=3,
        n_cole_cole=2,
        n_warburg=2,
        sparsity_weight=0.01,
        lr=0.02,
        n_epochs=4000,
        n_restarts=n_restarts
    )

    start_time = time.time()
    model = train_urn(freq, Z_data, config)
    train_time = time.time() - start_time

    # Get fitted response
    omega = 2 * np.pi * freq
    omega_tensor = torch.tensor(omega, dtype=torch.float64)
    Z_fit = model(omega_tensor).detach().numpy()

    # Count active parameters
    active_components = model.get_active_components()
    active_params = sum(len(v) for v in active_components.values())

    return Z_fit, model, active_params, train_time


def main():
    print("=" * 70)
    print("Comprehensive URN vs Vector Fitting Validation")
    print("=" * 70)

    # Define all datasets
    datasets = {
        'NASA 18650 Battery': {
            'path': 'data/real_world/nasa_battery/nasa_18650_eis.csv',
            'skiprows': 24,
            'n_poles': 10
        },
        'TDK PC47 Ferrite': {
            'path': 'data/real_world/tdk_ferrite/tdk_pc47_impedance.csv',
            'skiprows': 19,
            'n_poles': 12
        },
        'TDK PC50 Ferrite': {
            'path': 'data/real_world/tdk_ferrite/tdk_pc50_impedance.csv',
            'skiprows': 19,
            'n_poles': 12
        },
        'TDK PC95 Ferrite': {
            'path': 'data/real_world/tdk_ferrite/tdk_pc95_impedance.csv',
            'skiprows': 19,
            'n_poles': 12
        },
        'TDK PC200 Ferrite': {
            'path': 'data/real_world/tdk_ferrite/tdk_pc200_impedance.csv',
            'skiprows': 19,
            'n_poles': 12
        },
    }

    # Create output directory
    output_dir = Path(__file__).parent / 'validation_output'
    output_dir.mkdir(exist_ok=True)

    results = []

    for name, info in datasets.items():
        data_path = Path(__file__).parent / info['path']

        if not data_path.exists():
            print(f"\n[SKIP] {name}: Data file not found")
            continue

        print(f"\n[{len(results)+1}/{len(datasets)}] {name}")
        print("-" * 40)

        # Load data
        freq, Z_data = load_data(data_path, info['skiprows'])
        print(f"  Loaded {len(freq)} frequency points")
        print(f"  Frequency range: {freq.min():.4f} Hz - {freq.max()/1e6:.2f} MHz")

        # Run Vector Fitting
        print(f"  Running Vector Fitting ({info['n_poles']} poles)...")
        vf_start = time.time()
        Z_vf = vector_fitting(freq, Z_data, n_poles=info['n_poles'])
        vf_time = time.time() - vf_start
        vf_nrmse = compute_nrmse(Z_data, Z_vf)
        print(f"    NRMSE: {vf_nrmse:.4f}, Time: {vf_time:.2f}s")

        # Run URN fitting
        print(f"  Running URN fitting...")
        Z_urn, model, n_params, urn_time = run_urn_fitting(freq, Z_data)
        urn_nrmse = compute_nrmse(Z_data, Z_urn)
        improvement = (1 - urn_nrmse / vf_nrmse) * 100
        print(f"    NRMSE: {urn_nrmse:.4f}, Time: {urn_time:.2f}s")
        print(f"    Active parameters: {n_params}")
        print(f"    Improvement over VF: {improvement:.1f}%")

        # Store results
        results.append({
            'dataset': name,
            'vf_nrmse': float(vf_nrmse),
            'urn_nrmse': float(urn_nrmse),
            'improvement': float(improvement),
            'vf_time': float(vf_time),
            'urn_time': float(urn_time),
            'n_params': int(n_params),
            'n_points': int(len(freq))
        })

        # Generate plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Nyquist plot
        ax1 = axes[0]
        ax1.plot(Z_data.real, -Z_data.imag, 'ko', markersize=4, label='Measured')
        ax1.plot(Z_urn.real, -Z_urn.imag, 'b-', linewidth=2, label=f'URN (NRMSE={urn_nrmse:.4f})')
        ax1.plot(Z_vf.real, -Z_vf.imag, 'r--', linewidth=2, label=f'VF (NRMSE={vf_nrmse:.4f})')
        ax1.set_xlabel('Z\' (Ohm)')
        ax1.set_ylabel('-Z\'\' (Ohm)')
        ax1.set_title(f'{name} - Nyquist Plot')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal', adjustable='box')

        # Bode magnitude plot
        ax2 = axes[1]
        ax2.loglog(freq, np.abs(Z_data), 'ko', markersize=4, label='Measured')
        ax2.loglog(freq, np.abs(Z_urn), 'b-', linewidth=2, label='URN')
        ax2.loglog(freq, np.abs(Z_vf), 'r--', linewidth=2, label='VF')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('|Z| (Ohm)')
        ax2.set_title(f'{name} - Bode Magnitude')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_name = name.lower().replace(' ', '_')
        plt.savefig(output_dir / f'{plot_name}_comparison.png', dpi=150)
        plt.close()
        print(f"  Saved plot: {output_dir / f'{plot_name}_comparison.png'}")

    # Print summary table
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"\n{'Dataset':<25} {'VF NRMSE':>10} {'URN NRMSE':>10} {'Improvement':>12} {'URN Time':>10}")
    print("-" * 70)

    for r in results:
        print(f"{r['dataset']:<25} {r['vf_nrmse']:>10.4f} {r['urn_nrmse']:>10.4f} {r['improvement']:>11.1f}% {r['urn_time']:>9.1f}s")

    print("-" * 70)

    # Calculate averages
    avg_improvement = np.mean([r['improvement'] for r in results])
    avg_time = np.mean([r['urn_time'] for r in results])
    print(f"{'Average':<25} {'':<10} {'':<10} {avg_improvement:>11.1f}% {avg_time:>9.1f}s")

    # Save results to JSON
    results_file = output_dir / 'comprehensive_validation_results.json'
    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'results': results,
            'summary': {
                'n_datasets': len(results),
                'avg_improvement': float(avg_improvement),
                'avg_urn_time': float(avg_time)
            }
        }, f, indent=2)
    print(f"\nResults saved to: {results_file}")

    print("\nValidation complete!")


if __name__ == '__main__':
    main()
