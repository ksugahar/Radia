"""
URN Validation with Real-World Measured Data

This script validates the Universal Relaxation Network (URN) against real
measured impedance data from publicly available datasets.

Supported datasets:
1. NASA Li-ion Battery Aging Dataset
   - URL: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets
   - Format: MATLAB .mat files
   - Contains: EIS measurements at different aging stages

2. Mendeley SoC EIS Dataset (2024)
   - DOI: 10.17632/cb887gkmxw.2
   - Format: CSV files
   - Contains: EIS at 19 SoC levels for 11 LiFePO4 batteries

Usage:
    python validate_real_data.py --nasa-path data/real_world/nasa_battery/
    python validate_real_data.py --mendeley-path data/real_world/mendeley_eis/

Author: K. Sugahara, Y. Sato
Date: 2026-01-19
"""

import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from radia.urn import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)


def load_nasa_battery_eis(mat_file_path):
    """
    Load EIS data from NASA Battery Dataset MATLAB file.

    NASA Battery Dataset contains .mat files with structure:
    - data.Impedance.Re: Real part of impedance
    - data.Impedance.Im: Imaginary part of impedance
    - data.Impedance.freq: Frequency points

    Returns:
        freq: Frequency array (Hz)
        Z: Complex impedance array (Ohm)
        metadata: Dict with battery info
    """
    try:
        from scipy.io import loadmat
    except ImportError:
        print("ERROR: scipy is required to load .mat files")
        print("Install with: pip install scipy")
        sys.exit(1)

    print(f"Loading NASA battery data from: {mat_file_path}")

    # Load MATLAB file
    mat_data = loadmat(mat_file_path, struct_as_record=False, squeeze_me=True)

    # NASA dataset structure varies by file, try common patterns
    if 'data' in mat_data:
        data = mat_data['data']
        # Try to extract impedance data
        if hasattr(data, 'Impedance'):
            imp = data.Impedance
            freq = np.array(imp.freq)
            Z_re = np.array(imp.Re)
            Z_im = np.array(imp.Im)
            Z = Z_re + 1j * Z_im
            metadata = {
                'source': 'NASA Battery Dataset',
                'file': mat_file_path,
                'n_points': len(freq)
            }
            return freq, Z, metadata

    # Alternative structure (some NASA files)
    if 'freq' in mat_data and 'Zreal' in mat_data and 'Zimag' in mat_data:
        freq = np.array(mat_data['freq']).flatten()
        Z_re = np.array(mat_data['Zreal']).flatten()
        Z_im = np.array(mat_data['Zimag']).flatten()
        Z = Z_re + 1j * Z_im
        metadata = {
            'source': 'NASA Battery Dataset',
            'file': mat_file_path,
            'n_points': len(freq)
        }
        return freq, Z, metadata

    # Print available keys for debugging
    print(f"Available keys in .mat file: {list(mat_data.keys())}")
    raise ValueError(f"Could not parse NASA battery .mat file: {mat_file_path}")


def load_mendeley_eis(csv_file_path):
    """
    Load EIS data from Mendeley SoC EIS Dataset CSV file.

    Mendeley dataset CSV format:
    - Column 0: Frequency (Hz)
    - Column 1: Z_real (Ohm)
    - Column 2: Z_imag (Ohm) [typically negative for capacitive]

    Returns:
        freq: Frequency array (Hz)
        Z: Complex impedance array (Ohm)
        metadata: Dict with measurement info
    """
    print(f"Loading Mendeley EIS data from: {csv_file_path}")

    # Read CSV, skip header rows
    data = np.genfromtxt(csv_file_path, delimiter=',', skip_header=1)

    if data.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns, got {data.shape[1]}")

    freq = data[:, 0]
    Z_re = data[:, 1]
    Z_im = data[:, 2]
    Z = Z_re + 1j * Z_im

    metadata = {
        'source': 'Mendeley SoC EIS Dataset',
        'file': csv_file_path,
        'n_points': len(freq),
        'doi': '10.17632/cb887gkmxw.2'
    }

    return freq, Z, metadata


def run_urn_validation(freq, Z, metadata, output_dir=None):
    """
    Run URN validation on real measured data.

    Args:
        freq: Frequency array (Hz)
        Z: Complex impedance array (Ohm)
        metadata: Dict with data source info
        output_dir: Directory to save results (optional)

    Returns:
        Dict with validation results
    """
    print("\n" + "="*70)
    print("URN Validation on Real Measured Data")
    print("="*70)
    print(f"Data source: {metadata.get('source', 'Unknown')}")
    print(f"Number of points: {len(freq)}")
    print(f"Frequency range: {freq.min():.2e} - {freq.max():.2e} Hz")
    print(f"Impedance range: |Z| = {np.abs(Z).min():.4f} - {np.abs(Z).max():.4f} Ohm")
    print()

    # Configure URN for battery EIS (typical configuration)
    config = URNConfig(
        n_debye=3,
        n_cole_cole=2,
        n_warburg=2,
        n_cpe=2,
        n_gerischer=1,
        sparsity_weight=0.01,
        n_epochs=8000,
        n_restarts=8  # Increased for robustness as per reviewer request
    )

    # Train URN
    print("Training URN model...")
    model, history = train_urn(freq, Z, config, verbose=True)

    # Evaluate fit
    omega = 2 * np.pi * freq
    Z_pred = model.forward(omega)

    if isinstance(Z_pred, np.ndarray):
        Z_fit = Z_pred
    else:
        Z_fit = Z_pred.detach().cpu().numpy()

    # Calculate errors
    error_abs = np.abs(Z - Z_fit)
    error_rel = error_abs / np.abs(Z) * 100

    results = {
        'metadata': metadata,
        'freq': freq,
        'Z_measured': Z,
        'Z_fit': Z_fit,
        'error_mean': np.mean(error_rel),
        'error_max': np.max(error_rel),
        'error_std': np.std(error_rel),
        'n_mechanisms': len(model.get_discovered_mechanisms()),
        'mechanisms': model.get_discovered_mechanisms(),
        'training_history': history
    }

    print("\n" + "-"*50)
    print("Validation Results:")
    print("-"*50)
    print(f"Mean relative error: {results['error_mean']:.2f}%")
    print(f"Max relative error:  {results['error_max']:.2f}%")
    print(f"Std relative error:  {results['error_std']:.2f}%")
    print(f"Number of discovered mechanisms: {results['n_mechanisms']}")
    print("\nDiscovered mechanisms:")
    for name, params in results['mechanisms'].items():
        print(f"  {name}: {params}")

    # Generate SPICE netlist
    netlist = generate_spice_netlist(model, "REAL_DATA")
    results['spice_netlist'] = netlist
    print("\nSPICE Netlist Generated (first 20 lines):")
    for line in netlist.split('\n')[:20]:
        print(f"  {line}")

    # Plot results
    if output_dir:
        plot_validation_results(results, output_dir)

    return results


def plot_validation_results(results, output_dir):
    """Plot Nyquist and Bode diagrams for validation results."""
    freq = results['freq']
    Z = results['Z_measured']
    Z_fit = results['Z_fit']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Nyquist plot
    ax = axes[0, 0]
    ax.plot(Z.real, -Z.imag, 'b.-', label='Measured', markersize=3)
    ax.plot(Z_fit.real, -Z_fit.imag, 'r--', label='URN Fit', linewidth=2)
    ax.set_xlabel('Z_real (Ohm)')
    ax.set_ylabel('-Z_imag (Ohm)')
    ax.set_title('Nyquist Plot')
    ax.legend()
    ax.axis('equal')
    ax.grid(True, alpha=0.3)

    # Bode magnitude
    ax = axes[0, 1]
    ax.loglog(freq, np.abs(Z), 'b.-', label='Measured', markersize=3)
    ax.loglog(freq, np.abs(Z_fit), 'r--', label='URN Fit', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('|Z| (Ohm)')
    ax.set_title('Bode Magnitude')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Bode phase
    ax = axes[1, 0]
    phase_meas = np.angle(Z, deg=True)
    phase_fit = np.angle(Z_fit, deg=True)
    ax.semilogx(freq, phase_meas, 'b.-', label='Measured', markersize=3)
    ax.semilogx(freq, phase_fit, 'r--', label='URN Fit', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Phase (deg)')
    ax.set_title('Bode Phase')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Error plot
    ax = axes[1, 1]
    error_rel = np.abs(Z - Z_fit) / np.abs(Z) * 100
    ax.semilogx(freq, error_rel, 'g.-', markersize=3)
    ax.axhline(y=results['error_mean'], color='r', linestyle='--',
               label=f'Mean: {results["error_mean"]:.2f}%')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Relative Error (%)')
    ax.set_title('Fitting Error')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    source = results['metadata'].get('source', 'unknown').replace(' ', '_')
    fig_path = os.path.join(output_dir, f'urn_validation_{source}.png')
    plt.savefig(fig_path, dpi=150)
    print(f"\nValidation plot saved to: {fig_path}")
    plt.close()


def run_convergence_analysis(freq, Z, output_dir=None):
    """
    Analyze convergence vs number of restarts (reviewer request 2.3).

    Tests n_restarts = 1, 3, 5, 8, 12, 16 and reports:
    - Convergence probability (% runs achieving < 5% error)
    - Error variance across runs
    - Best/worst case errors
    """
    print("\n" + "="*70)
    print("Convergence Analysis: n_restarts Sensitivity")
    print("="*70)

    restart_values = [1, 3, 5, 8, 12, 16]
    results_by_restart = {}

    for n_restart in restart_values:
        print(f"\nTesting n_restarts = {n_restart}...")

        config = URNConfig(
            n_debye=2,
            n_cole_cole=2,
            n_warburg=1,
            sparsity_weight=0.01,
            n_epochs=5000,
            n_restarts=n_restart
        )

        # Run multiple trials
        n_trials = 5
        errors = []

        for trial in range(n_trials):
            model, _ = train_urn(freq, Z, config, verbose=False)
            omega = 2 * np.pi * freq
            Z_pred = model.forward(omega)
            if not isinstance(Z_pred, np.ndarray):
                Z_pred = Z_pred.detach().cpu().numpy()

            error = np.mean(np.abs(Z - Z_pred) / np.abs(Z)) * 100
            errors.append(error)
            print(f"  Trial {trial+1}: {error:.2f}% error")

        errors = np.array(errors)
        results_by_restart[n_restart] = {
            'mean_error': np.mean(errors),
            'std_error': np.std(errors),
            'min_error': np.min(errors),
            'max_error': np.max(errors),
            'convergence_prob': np.mean(errors < 5.0) * 100
        }

    # Print summary table
    print("\n" + "-"*70)
    print("Convergence Analysis Summary")
    print("-"*70)
    print(f"{'n_restarts':>12} {'Mean Err(%)':>12} {'Std(%)':>10} {'Min(%)':>10} {'Max(%)':>10} {'Conv.Prob':>12}")
    print("-"*70)

    for n_restart in restart_values:
        r = results_by_restart[n_restart]
        print(f"{n_restart:>12} {r['mean_error']:>12.2f} {r['std_error']:>10.2f} "
              f"{r['min_error']:>10.2f} {r['max_error']:>10.2f} {r['convergence_prob']:>11.0f}%")

    return results_by_restart


def run_sparsity_sensitivity(freq, Z, output_dir=None):
    """
    Analyze sparsity weight sensitivity (reviewer request 2.3).

    Tests lambda = 0.001, 0.005, 0.01, 0.02, 0.05, 0.1 and reports:
    - Fitting error
    - Number of discovered mechanisms
    - Discovered mechanism types
    """
    print("\n" + "="*70)
    print("Hyperparameter Sensitivity: sparsity_weight")
    print("="*70)

    lambda_values = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    results_by_lambda = {}

    for lam in lambda_values:
        print(f"\nTesting sparsity_weight = {lam}...")

        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_warburg=2,
            n_cpe=2,
            sparsity_weight=lam,
            n_epochs=6000,
            n_restarts=5
        )

        model, _ = train_urn(freq, Z, config, verbose=False)
        omega = 2 * np.pi * freq
        Z_pred = model.forward(omega)
        if not isinstance(Z_pred, np.ndarray):
            Z_pred = Z_pred.detach().cpu().numpy()

        error = np.mean(np.abs(Z - Z_pred) / np.abs(Z)) * 100
        mechanisms = model.get_discovered_mechanisms()

        results_by_lambda[lam] = {
            'error': error,
            'n_mechanisms': len(mechanisms),
            'mechanisms': list(mechanisms.keys())
        }

        print(f"  Error: {error:.2f}%, Mechanisms: {len(mechanisms)}")
        for name in mechanisms:
            print(f"    - {name}")

    # Print summary table
    print("\n" + "-"*70)
    print("Sparsity Weight Sensitivity Summary")
    print("-"*70)
    print(f"{'lambda':>12} {'Error(%)':>12} {'# Mechs':>10} {'Types'}")
    print("-"*70)

    for lam in lambda_values:
        r = results_by_lambda[lam]
        types_str = ', '.join(r['mechanisms'][:3])
        if len(r['mechanisms']) > 3:
            types_str += '...'
        print(f"{lam:>12.3f} {r['error']:>12.2f} {r['n_mechanisms']:>10} {types_str}")

    return results_by_lambda


def run_noise_robustness_test(freq, Z, output_dir=None):
    """
    Test noise robustness (reviewer request 2.4.A.3).

    Adds Gaussian noise at SNR = inf, 40, 20, 10, 0 dB and evaluates:
    - Parameter recovery accuracy
    - Topology consistency (same mechanisms discovered)
    """
    print("\n" + "="*70)
    print("Noise Robustness Test")
    print("="*70)

    snr_values = [float('inf'), 40, 20, 10, 0]
    results_by_snr = {}

    # Reference: clean data fit
    config = URNConfig(
        n_debye=2,
        n_cole_cole=2,
        n_warburg=1,
        sparsity_weight=0.01,
        n_epochs=6000,
        n_restarts=5
    )

    for snr in snr_values:
        print(f"\nTesting SNR = {snr} dB...")

        # Add noise
        if snr == float('inf'):
            Z_noisy = Z.copy()
        else:
            signal_power = np.mean(np.abs(Z)**2)
            noise_power = signal_power / (10**(snr/10))
            noise = np.sqrt(noise_power/2) * (np.random.randn(len(Z)) + 1j*np.random.randn(len(Z)))
            Z_noisy = Z + noise

        model, _ = train_urn(freq, Z_noisy, config, verbose=False)
        omega = 2 * np.pi * freq
        Z_pred = model.forward(omega)
        if not isinstance(Z_pred, np.ndarray):
            Z_pred = Z_pred.detach().cpu().numpy()

        # Error vs clean data (not noisy)
        error_vs_clean = np.mean(np.abs(Z - Z_pred) / np.abs(Z)) * 100
        mechanisms = model.get_discovered_mechanisms()

        results_by_snr[snr] = {
            'error_vs_clean': error_vs_clean,
            'n_mechanisms': len(mechanisms),
            'mechanisms': list(mechanisms.keys())
        }

        print(f"  Error vs clean: {error_vs_clean:.2f}%, Mechanisms: {len(mechanisms)}")

    # Check topology consistency
    ref_mechanisms = set(results_by_snr[float('inf')]['mechanisms'])

    print("\n" + "-"*70)
    print("Noise Robustness Summary")
    print("-"*70)
    print(f"{'SNR (dB)':>12} {'Error(%)':>12} {'# Mechs':>10} {'Topology Match'}")
    print("-"*70)

    for snr in snr_values:
        r = results_by_snr[snr]
        snr_str = 'inf' if snr == float('inf') else f'{snr}'
        current_mechs = set(r['mechanisms'])
        match = 'Yes' if current_mechs == ref_mechanisms else 'Partial'
        print(f"{snr_str:>12} {r['error_vs_clean']:>12.2f} {r['n_mechanisms']:>10} {match}")

    return results_by_snr


def main():
    parser = argparse.ArgumentParser(description='URN Validation with Real-World Data')
    parser.add_argument('--nasa-path', type=str, help='Path to NASA battery .mat file')
    parser.add_argument('--mendeley-path', type=str, help='Path to Mendeley EIS .csv file')
    parser.add_argument('--synthetic', action='store_true',
                        help='Use synthetic data for testing (no real data required)')
    parser.add_argument('--output-dir', type=str, default='results',
                        help='Directory for output files')
    parser.add_argument('--convergence-test', action='store_true',
                        help='Run convergence analysis (n_restarts sensitivity)')
    parser.add_argument('--sparsity-test', action='store_true',
                        help='Run sparsity weight sensitivity analysis')
    parser.add_argument('--noise-test', action='store_true',
                        help='Run noise robustness test')
    parser.add_argument('--all-tests', action='store_true',
                        help='Run all theoretical analyses')

    args = parser.parse_args()

    # Load data
    if args.nasa_path:
        freq, Z, metadata = load_nasa_battery_eis(args.nasa_path)
    elif args.mendeley_path:
        freq, Z, metadata = load_mendeley_eis(args.mendeley_path)
    elif args.synthetic:
        # Use synthetic data for testing
        print("Using synthetic battery EIS data for demonstration...")
        data_path = Path(__file__).parent / 'data' / 'synthetic' / 'liion_battery_eis.csv'
        if not data_path.exists():
            print(f"ERROR: Synthetic data not found at {data_path}")
            sys.exit(1)
        data = np.genfromtxt(data_path, delimiter=',', skip_header=18)
        freq = data[:, 0]
        Z = data[:, 1] + 1j * data[:, 2]
        metadata = {
            'source': 'Synthetic Battery EIS',
            'file': str(data_path),
            'n_points': len(freq),
            'note': 'Physics-based synthetic data for algorithm validation'
        }
    else:
        print("ERROR: Please specify data source:")
        print("  --nasa-path <path>     : NASA battery .mat file")
        print("  --mendeley-path <path> : Mendeley EIS .csv file")
        print("  --synthetic            : Use bundled synthetic data")
        sys.exit(1)

    # Run main validation
    results = run_urn_validation(freq, Z, metadata, args.output_dir)

    # Run additional theoretical analyses if requested
    if args.all_tests or args.convergence_test:
        run_convergence_analysis(freq, Z, args.output_dir)

    if args.all_tests or args.sparsity_test:
        run_sparsity_sensitivity(freq, Z, args.output_dir)

    if args.all_tests or args.noise_test:
        run_noise_robustness_test(freq, Z, args.output_dir)

    print("\n" + "="*70)
    print("Validation Complete")
    print("="*70)

    return results


if __name__ == '__main__':
    main()
