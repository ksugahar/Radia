#!/usr/bin/env python3
"""
Validate URN with TDK MnZn Ferrite Real Measurement Data

This script validates the Universal Relaxation Network (URN) using REAL manufacturer
data from TDK Corporation's MnZn ferrite datasheet.

Data Source:
    TDK Corporation "Mn-Zn Ferrite Material characteristics" (May 2022)
    Document: ferrite_mn-zn_material_characteristics_en.pdf
    Materials: PC50, PC200 (high-frequency power supply ferrites)

This addresses the "Lack of real-world validation" critique in REVIEW.md by using
actual manufacturer measurement data rather than synthetic benchmarks.

Author: Radia Project
Date: 2026-01-19
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

import torch

script_dir = os.path.dirname(os.path.abspath(__file__))

from radia.urn import (
    URNConfig,
    generate_spice_netlist,
    train_urn,
)


def load_tdk_impedance_data(material: str = 'PC50') -> Tuple[np.ndarray, np.ndarray]:
    """
    Load TDK ferrite impedance data from CSV.

    Args:
        material: 'PC50' or 'PC200'

    Returns:
        freq: Frequency array [Hz]
        Z: Complex impedance array [Ohm]
    """
    data_dir = os.path.join(script_dir, 'data', 'real_world', 'tdk_ferrite')
    csv_file = os.path.join(data_dir, f'tdk_{material.lower()}_impedance.csv')

    if not os.path.exists(csv_file):
        raise FileNotFoundError(
            f"TDK impedance data not found: {csv_file}\n"
            f"Run convert_permeability_to_impedance.py first."
        )

    df = pd.read_csv(csv_file, comment='#')
    freq = df['frequency_Hz'].values
    Z_real = df['Z_real_Ohm'].values
    Z_imag = df['Z_imag_Ohm'].values
    Z = Z_real + 1j * Z_imag

    return freq, Z


def plot_nyquist_comparison(freq: np.ndarray, Z_data: np.ndarray, Z_urn: np.ndarray,
                            material: str, output_dir: str):
    """Plot Nyquist diagram comparing data and URN fit."""
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(Z_data.real, -Z_data.imag, 'bo-', label=f'TDK {material} (measured)', markersize=8)
    ax.plot(Z_urn.real, -Z_urn.imag, 'r-', label='URN fit', linewidth=2)

    ax.set_xlabel('Re(Z) [Ohm]', fontsize=12)
    ax.set_ylabel('-Im(Z) [Ohm]', fontsize=12)
    ax.set_title(f'Nyquist Plot: TDK {material} Ferrite\n(Real Manufacturer Data)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'tdk_{material.lower()}_nyquist.png'), dpi=150)
    plt.close()


def plot_bode_comparison(freq: np.ndarray, Z_data: np.ndarray, Z_urn: np.ndarray,
                         material: str, output_dir: str):
    """Plot Bode diagram comparing data and URN fit."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Magnitude
    ax1.loglog(freq, np.abs(Z_data), 'bo-', label=f'TDK {material} (measured)', markersize=6)
    ax1.loglog(freq, np.abs(Z_urn), 'r-', label='URN fit', linewidth=2)
    ax1.set_ylabel('|Z| [Ohm]', fontsize=12)
    ax1.set_title(f'Bode Plot: TDK {material} Ferrite\n(Real Manufacturer Data)', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, which='both')

    # Phase
    phase_data = np.degrees(np.angle(Z_data))
    phase_urn = np.degrees(np.angle(Z_urn))
    ax2.semilogx(freq, phase_data, 'bo-', label=f'TDK {material} (measured)', markersize=6)
    ax2.semilogx(freq, phase_urn, 'r-', label='URN fit', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]', fontsize=12)
    ax2.set_ylabel('Phase [deg]', fontsize=12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'tdk_{material.lower()}_bode.png'), dpi=150)
    plt.close()


def plot_mu_comparison(freq: np.ndarray, Z_data: np.ndarray, Z_urn: np.ndarray,
                       material: str, output_dir: str):
    """
    Plot complex permeability comparison (convert back from impedance).

    This shows how well URN captures the mu'(f) and mu"(f) behavior.
    """
    # Core parameters (same as conversion script)
    MU0 = 4 * np.pi * 1e-7
    OD, ID, TH = 31e-3, 19e-3, 8e-3
    N_TURNS = 10
    Ae = TH * (OD - ID) / 2
    le = np.pi * (OD + ID) / 2
    L0 = MU0 * N_TURNS**2 * Ae / le

    omega = 2 * np.pi * freq

    # Convert impedance back to permeability
    # Z = omega*L0*mu" + j*omega*L0*mu'
    # mu' = Im(Z) / (omega*L0)
    # mu" = Re(Z) / (omega*L0)
    mu_real_data = Z_data.imag / (omega * L0)
    mu_imag_data = Z_data.real / (omega * L0)
    mu_real_urn = Z_urn.imag / (omega * L0)
    mu_imag_urn = Z_urn.real / (omega * L0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # mu' (real part)
    ax1.semilogx(freq, mu_real_data, 'bo-', label=f"TDK {material} mu' (measured)", markersize=6)
    ax1.semilogx(freq, mu_real_urn, 'r-', label="URN fit mu'", linewidth=2)
    ax1.set_ylabel("mu' (real part)", fontsize=12)
    ax1.set_title(f"Complex Permeability: TDK {material} Ferrite\n(Real Manufacturer Data)", fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, which='both')

    # mu" (imaginary part)
    ax2.semilogx(freq, mu_imag_data, 'bs-', label=f'TDK {material} mu" (measured)', markersize=6)
    ax2.semilogx(freq, mu_imag_urn, 'r--', label='URN fit mu"', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]', fontsize=12)
    ax2.set_ylabel('mu" (imaginary part)', fontsize=12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'tdk_{material.lower()}_permeability.png'), dpi=150)
    plt.close()


def compute_error_metrics(Z_data: np.ndarray, Z_urn: np.ndarray) -> Dict:
    """Compute error metrics between data and URN fit."""
    # Relative error
    rel_error = np.abs(Z_urn - Z_data) / np.abs(Z_data)

    # Phase error
    phase_data = np.angle(Z_data)
    phase_urn = np.angle(Z_urn)
    phase_error = np.abs(phase_urn - phase_data)

    return {
        'rel_error_mean': float(np.mean(rel_error)),
        'rel_error_max': float(np.max(rel_error)),
        'rel_error_std': float(np.std(rel_error)),
        'phase_error_mean_deg': float(np.degrees(np.mean(phase_error))),
        'phase_error_max_deg': float(np.degrees(np.max(phase_error))),
    }


def validate_tdk_ferrite(material: str = 'PC50', output_dir: str = None):
    """
    Run URN validation on TDK ferrite data.

    Args:
        material: 'PC50' or 'PC200'
        output_dir: Output directory for results
    """
    print("="*70)
    print(f"URN Validation with TDK {material} Ferrite (Real Manufacturer Data)")
    print("="*70)
    print()
    print("Data Source:")
    print("  TDK Corporation 'Mn-Zn Ferrite Material characteristics' (May 2022)")
    print("  Document: ferrite_mn-zn_material_characteristics_en.pdf")
    print()

    # Setup output directory
    if output_dir is None:
        repo_root = Path(script_dir).parents[1]
        output_dir = str(
            repo_root / 'validation_test' / 'universal_relaxation_network'
            / 'tdk_ferrite'
        )
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    print(f"Loading TDK {material} impedance data...")
    freq, Z_data = load_tdk_impedance_data(material)
    print(f"  Frequency range: {freq[0]/1e3:.1f} kHz - {freq[-1]/1e6:.1f} MHz")
    print(f"  Number of points: {len(freq)}")
    print(f"  |Z| range: {np.abs(Z_data).min():.2f} - {np.abs(Z_data).max():.1f} Ohm")
    print()

    # Configure URN for ferrite (magnetic material)
    # Ferrite exhibits magnetic relaxation -> use magnetic basis functions
    config = URNConfig(
        # Magnetic relaxation basis functions (primary)
        n_debye=3,           # Debye relaxation
        n_cole_cole=2,       # Cole-Cole (distributed relaxation)
        n_havriliak_negami=1,  # Havriliak-Negami (asymmetric)

        # Skin effect (secondary, for high frequency)
        n_skin_effect=1,

        # Disable non-magnetic basis
        n_cpe=1,             # CPE for non-ideal behavior
        n_warburg=0,         # No diffusion in ferrite
        n_gerischer=0,
        n_rlc=1,             # Resonance possible at very high freq

        # Parallel branch (minimal)
        n_debye_parallel=1,
        n_cole_cole_parallel=1,
        n_cpe_parallel=0,
        n_warburg_parallel=0,

        # Training parameters
        sparsity_weight=0.005,  # Lower sparsity for better fit
        lr=0.02,
        n_epochs=8000,
        n_restarts=5,
    )

    print("URN Configuration:")
    print(f"  Series branch: {config.n_debye} Debye + {config.n_cole_cole} Cole-Cole + "
          f"{config.n_havriliak_negami} HN + {config.n_skin_effect} Skin effect")
    print(f"  Parallel branch: {config.n_debye_parallel} Debye + {config.n_cole_cole_parallel} Cole-Cole")
    print(f"  Training: {config.n_epochs} epochs x {config.n_restarts} restarts")
    print()

    # Fit URN
    print("Fitting URN model...")
    urn_model = train_urn(freq, Z_data, config)

    # Evaluate
    omega = 2 * np.pi * freq
    with torch.no_grad():
        Z_urn = urn_model(
            torch.tensor(omega, dtype=torch.float64)
        ).detach().cpu().numpy()

    # Compute metrics
    print("\nComputing error metrics...")
    metrics = compute_error_metrics(Z_data, Z_urn)
    print(f"  Mean relative error: {metrics['rel_error_mean']*100:.2f}%")
    print(f"  Max relative error: {metrics['rel_error_max']*100:.2f}%")
    print(f"  Mean phase error: {metrics['phase_error_mean_deg']:.2f} deg")
    print()

    # Generate plots
    print("Generating comparison plots...")
    plot_nyquist_comparison(freq, Z_data, Z_urn, material, output_dir)
    plot_bode_comparison(freq, Z_data, Z_urn, material, output_dir)
    plot_mu_comparison(freq, Z_data, Z_urn, material, output_dir)
    print(f"  Plots saved to: {output_dir}")

    # Generate SPICE netlist
    print("\nGenerating SPICE netlist...")
    netlist = generate_spice_netlist(urn_model, f"TDK_{material}_URN")
    netlist_file = os.path.join(output_dir, f'tdk_{material.lower()}_urn.sp')
    with open(netlist_file, 'w') as f:
        f.write(netlist)
    print(f"  SPICE netlist: {netlist_file}")

    # Get active mechanisms
    print("\nDiscovered relaxation mechanisms:")
    components = urn_model.get_active_components(threshold=0.01)
    active = sorted(
        (
            (f"{family}[{index}]", float(component['weight_magnitude']))
            for family, family_components in components.items()
            for index, component in enumerate(family_components)
        ),
        key=lambda item: item[1],
        reverse=True,
    )
    for name, weight in active[:5]:
        print(f"  - {name}: weight = {weight:.4f}")

    # Save results JSON
    results = {
        'material': material,
        'data_source': 'TDK Corporation datasheet (May 2022)',
        'document': 'ferrite_mn-zn_material_characteristics_en.pdf',
        'data_type': 'REAL MANUFACTURER DATA',
        'timestamp': datetime.now().isoformat(),
        'frequency_range_Hz': [float(freq[0]), float(freq[-1])],
        'n_points': len(freq),
        'metrics': metrics,
        'active_mechanisms': [(name, float(w)) for name, w in active],
        'config': {
            'n_debye': config.n_debye,
            'n_cole_cole': config.n_cole_cole,
            'n_epochs': config.n_epochs,
            'n_restarts': config.n_restarts,
            'sparsity_weight': config.sparsity_weight,
        }
    }
    results_file = os.path.join(output_dir, f'tdk_{material.lower()}_results.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results JSON: {results_file}")

    print("\n" + "="*70)
    print("Validation complete!")
    print("="*70)

    return results


def main():
    """Run validation on both TDK ferrite materials."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Validate URN with TDK MnZn Ferrite real measurement data')
    parser.add_argument('--material', type=str, default='both',
                        choices=['PC50', 'PC200', 'both'],
                        help='TDK ferrite material to validate')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for results')

    args = parser.parse_args()

    materials = ['PC50', 'PC200'] if args.material == 'both' else [args.material]

    all_results = {}
    for material in materials:
        try:
            results = validate_tdk_ferrite(material, args.output_dir)
            all_results[material] = results
            print()
        except Exception as e:
            print(f"Error validating {material}: {e}")
            import traceback
            traceback.print_exc()

    # Summary
    print("\n" + "="*70)
    print("SUMMARY: URN Validation with TDK Real Manufacturer Data")
    print("="*70)
    for material, results in all_results.items():
        metrics = results['metrics']
        print(f"\n{material}:")
        print(f"  Mean relative error: {metrics['rel_error_mean']*100:.2f}%")
        print(f"  Max relative error: {metrics['rel_error_max']*100:.2f}%")
        print(f"  Top mechanism: {results['active_mechanisms'][0][0]}")

    print("\n" + "-"*70)
    print("This validation uses REAL MANUFACTURER DATA from TDK Corporation,")
    print("addressing the 'Lack of real-world validation' critique in REVIEW.md.")
    print("-"*70)


if __name__ == '__main__':
    main()
