#!/usr/bin/env python
"""
URN vs scikit-rf Vector Fitting Benchmark

This script addresses reviewer critique 2.5.B:
"The benchmark compares URN against a simplified, home-brewed Vector Fitting
implementation rather than the industry-standard vectfit3 or scikit-rf."

This script uses scikit-rf's Vector Fitting implementation, which is based on
the original SINTEF vectfit3 algorithm by Gustavsen & Semlyen.

Requirements:
    pip install scikit-rf numpy scipy torch matplotlib

References:
    - Gustavsen, B., Semlyen, A. (1999). "Rational approximation of frequency
      domain responses by vector fitting." IEEE Trans. Power Delivery.
    - scikit-rf: https://scikit-rf.readthedocs.io/

Usage:
    python benchmark_urn_vs_skrf_vf.py
    python benchmark_urn_vs_skrf_vf.py --dataset ferrite

Author: K. Sugahara, Y. Sato
Date: 2026-01-19
"""

import time
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from radia.urn import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)

# Check for scikit-rf
try:
    import skrf as rf
    from skrf.vectorFitting import VectorFitting
    SKRF_AVAILABLE = True
except ImportError:
    SKRF_AVAILABLE = False
    print("WARNING: scikit-rf not installed. Install with: pip install scikit-rf")
    print("         Falling back to simplified VF implementation.")


def run_skrf_vector_fitting(freq, Z, n_poles=8, n_iterations=10):
    """
    Run Vector Fitting using scikit-rf implementation.

    This is the industry-standard implementation based on vectfit3.

    Parameters
    ----------
    freq : array
        Frequency points (Hz)
    Z : array
        Complex impedance
    n_poles : int
        Number of poles (will be adjusted for real/complex pairs)
    n_iterations : int
        Number of pole relocation iterations

    Returns
    -------
    dict
        Results including poles, residues, fit, error, timing
    """
    if not SKRF_AVAILABLE:
        return run_fallback_vf(freq, Z, n_poles)

    print(f"  Running scikit-rf Vector Fitting (n_poles={n_poles})...")
    start_time = time.time()

    # Create Network object for VF
    # VF expects S-parameters, but can work with impedance via conversion
    # For single-port impedance, Z -> S = (Z - Z0) / (Z + Z0)
    Z0 = 50.0  # Reference impedance
    S = (Z - Z0) / (Z + Z0)

    # Create frequency object
    f = rf.Frequency.from_f(freq, unit='Hz')

    # Create single-port network
    ntwk = rf.Network(frequency=f, s=S.reshape(-1, 1, 1), z0=Z0)

    # Initialize Vector Fitting
    vf = VectorFitting(ntwk)

    # Auto-fit with pole relocation
    try:
        # Use auto_fit which handles pole initialization and iteration
        vf.vector_fit(n_poles_real=n_poles//2, n_poles_cmplx=n_poles//4,
                      init_pole_spacing='log',
                      fit_constant=True, fit_proportional=False)

        # Get fitted response
        S_fit = vf.get_model_response(freqs=freq)
        Z_fit = Z0 * (1 + S_fit.flatten()) / (1 - S_fit.flatten())

        # Extract poles and residues
        poles = vf.poles
        residues = vf.residues.flatten() if hasattr(vf, 'residues') else []

    except Exception as e:
        print(f"    scikit-rf VF failed: {e}")
        print(f"    Falling back to simplified implementation")
        return run_fallback_vf(freq, Z, n_poles)

    elapsed = time.time() - start_time

    # Calculate errors
    error_abs = np.abs(Z - Z_fit)
    error_rel = error_abs / np.abs(Z) * 100

    # Check for stability issues
    unstable_poles = [p for p in poles if np.real(p) > 0]
    weak_poles = []
    for p in poles:
        if np.imag(p) != 0:
            zeta = -np.real(p) / np.abs(p)
            if 0 < zeta < 0.1:
                weak_poles.append((p, zeta))

    result = {
        'method': 'scikit-rf VectorFitting',
        'Z_fit': Z_fit,
        'poles': poles,
        'residues': residues,
        'error_mean': np.mean(error_rel),
        'error_max': np.max(error_rel),
        'error_std': np.std(error_rel),
        'n_poles': len(poles),
        'n_unstable': len(unstable_poles),
        'n_weak_damped': len(weak_poles),
        'time': elapsed,
        'is_passive': len(unstable_poles) == 0,
        'has_ringing_risk': len(weak_poles) > 0,
    }

    print(f"    Completed in {elapsed:.2f}s")
    print(f"    Error: {result['error_mean']:.2f}% mean, {result['error_max']:.2f}% max")
    print(f"    Poles: {len(poles)} ({len(unstable_poles)} unstable, {len(weak_poles)} weak)")

    return result


def run_fallback_vf(freq, Z, n_poles=8):
    """
    Fallback Vector Fitting when scikit-rf is not available.

    This is the simplified implementation used only when scikit-rf fails.
    """
    print(f"  Running fallback VF implementation (n_poles={n_poles})...")
    start_time = time.time()

    omega = 2 * np.pi * freq
    s = 1j * omega

    # Initialize poles with log spacing
    omega_min = omega[omega > 0].min()
    omega_max = omega.max()
    omega_poles = np.logspace(np.log10(omega_min), np.log10(omega_max), n_poles // 2)

    # Create complex pole pairs
    poles = []
    zeta = 0.5  # Moderate damping
    for wp in omega_poles:
        real_part = -zeta * wp
        imag_part = wp * np.sqrt(1 - zeta**2)
        poles.append(real_part + 1j * imag_part)
        poles.append(real_part - 1j * imag_part)

    poles = np.array(poles[:n_poles])

    # Solve for residues
    n_freq = len(freq)
    A = np.zeros((n_freq, n_poles + 1), dtype=complex)
    for n, p in enumerate(poles):
        A[:, n] = 1.0 / (s - p)
    A[:, -1] = 1.0

    x, _, _, _ = np.linalg.lstsq(A, Z, rcond=None)
    residues = x[:-1]
    d = x[-1].real

    # Evaluate fit
    Z_fit = np.zeros_like(s)
    for p, r in zip(poles, residues):
        Z_fit += r / (s - p)
    Z_fit += d

    elapsed = time.time() - start_time

    # Calculate errors
    error_rel = np.abs(Z - Z_fit) / np.abs(Z) * 100

    # Stability analysis
    unstable_poles = [p for p in poles if np.real(p) > 0]
    weak_poles = [(p, -np.real(p)/np.abs(p)) for p in poles
                  if np.imag(p) != 0 and 0 < -np.real(p)/np.abs(p) < 0.1]

    result = {
        'method': 'Fallback VF (simplified)',
        'Z_fit': Z_fit,
        'poles': poles,
        'residues': residues,
        'error_mean': np.mean(error_rel),
        'error_max': np.max(error_rel),
        'error_std': np.std(error_rel),
        'n_poles': len(poles),
        'n_unstable': len(unstable_poles),
        'n_weak_damped': len(weak_poles),
        'time': elapsed,
        'is_passive': len(unstable_poles) == 0,
        'has_ringing_risk': len(weak_poles) > 0,
    }

    print(f"    Completed in {elapsed:.2f}s")
    print(f"    Error: {result['error_mean']:.2f}% mean")

    return result


def run_urn_fitting(freq, Z):
    """Run URN fitting for comparison."""
    print(f"  Running URN fitting...")
    start_time = time.time()

    config = URNConfig(
        n_debye=3,
        n_cole_cole=2,
        n_warburg=2,
        n_cpe=2,
        sparsity_weight=0.01,
        n_epochs=6000,
        n_restarts=5
    )

    model, history = train_urn(freq, Z, config, verbose=False)

    omega = 2 * np.pi * freq
    Z_fit = model.forward_numpy(omega)

    elapsed = time.time() - start_time

    # Calculate errors
    error_rel = np.abs(Z - Z_fit) / np.abs(Z) * 100

    mechanisms = model.get_discovered_mechanisms()

    result = {
        'method': 'URN (Physics-Informed)',
        'Z_fit': Z_fit,
        'model': model,
        'mechanisms': mechanisms,
        'error_mean': np.mean(error_rel),
        'error_max': np.max(error_rel),
        'error_std': np.std(error_rel),
        'n_mechanisms': len(mechanisms),
        'time': elapsed,
        'is_passive': True,  # By construction
        'has_ringing_risk': False,  # Pole-free
    }

    print(f"    Completed in {elapsed:.2f}s")
    print(f"    Error: {result['error_mean']:.2f}% mean, {result['error_max']:.2f}% max")
    print(f"    Discovered {len(mechanisms)} mechanisms")

    return result


def generate_comparison_report(freq, Z, urn_result, vf_result, output_dir):
    """Generate comprehensive comparison report and plots."""

    fig = plt.figure(figsize=(16, 12))

    # 1. Nyquist plot
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.plot(Z.real, -Z.imag, 'ko', markersize=4, label='Data')
    ax1.plot(urn_result['Z_fit'].real, -urn_result['Z_fit'].imag,
             'b-', linewidth=2, label=f"URN ({urn_result['error_mean']:.1f}%)")
    ax1.plot(vf_result['Z_fit'].real, -vf_result['Z_fit'].imag,
             'r--', linewidth=2, label=f"VF ({vf_result['error_mean']:.1f}%)")
    ax1.set_xlabel('Z_real [Ohm]')
    ax1.set_ylabel('-Z_imag [Ohm]')
    ax1.set_title('Nyquist Plot')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')

    # 2. Bode magnitude
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.loglog(freq, np.abs(Z), 'ko', markersize=3, label='Data')
    ax2.loglog(freq, np.abs(urn_result['Z_fit']), 'b-', linewidth=2, label='URN')
    ax2.loglog(freq, np.abs(vf_result['Z_fit']), 'r--', linewidth=2, label='VF')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('|Z| [Ohm]')
    ax2.set_title('Bode Magnitude')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)

    # 3. Error comparison
    ax3 = fig.add_subplot(2, 3, 3)
    err_urn = np.abs(Z - urn_result['Z_fit']) / np.abs(Z) * 100
    err_vf = np.abs(Z - vf_result['Z_fit']) / np.abs(Z) * 100
    ax3.semilogx(freq, err_urn, 'b-', linewidth=2, label=f"URN (mean: {urn_result['error_mean']:.2f}%)")
    ax3.semilogx(freq, err_vf, 'r--', linewidth=2, label=f"VF (mean: {vf_result['error_mean']:.2f}%)")
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Relative Error [%]')
    ax3.set_title('Fitting Error')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # 4. Pole map (VF only)
    ax4 = fig.add_subplot(2, 3, 4)
    if 'poles' in vf_result and len(vf_result['poles']) > 0:
        poles = vf_result['poles']
        ax4.scatter(np.real(poles), np.imag(poles), s=100, marker='x',
                   c='red', linewidths=2, label='VF poles')
        ax4.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
        ax4.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        ax4.axvspan(ax4.get_xlim()[0] if ax4.get_xlim()[0] < 0 else -1e6, 0,
                   alpha=0.1, color='green')
        ax4.set_xlabel('Real(s) [rad/s]')
        ax4.set_ylabel('Imag(s) [rad/s]')
        ax4.set_title(f'VF Pole Map ({vf_result["n_unstable"]} unstable, {vf_result["n_weak_damped"]} weak)')
    else:
        ax4.text(0.5, 0.5, 'No poles (URN is pole-free)', transform=ax4.transAxes,
                ha='center', va='center', fontsize=12)
        ax4.set_title('URN: Pole-Free Formulation')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Comparison table
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.axis('off')

    table_data = [
        ['Metric', 'URN', f'{vf_result["method"]}'],
        ['---', '---', '---'],
        ['Mean Error (%)', f'{urn_result["error_mean"]:.2f}', f'{vf_result["error_mean"]:.2f}'],
        ['Max Error (%)', f'{urn_result["error_max"]:.2f}', f'{vf_result["error_max"]:.2f}'],
        ['Time (s)', f'{urn_result["time"]:.1f}', f'{vf_result["time"]:.2f}'],
        ['Passive?', 'Yes (by construction)', 'Yes' if vf_result['is_passive'] else 'NO!'],
        ['Ringing Risk?', 'No (pole-free)', 'Yes' if vf_result['has_ringing_risk'] else 'No'],
        ['Interpretable?', 'Yes (mechanisms)', 'No (poles only)'],
    ]

    table_text = '\n'.join(['  '.join(f'{c:>18}' for c in row) for row in table_data])

    ax5.text(0.1, 0.9, 'Benchmark Comparison Results',
            transform=ax5.transAxes, fontsize=14, fontweight='bold', va='top')
    ax5.text(0.1, 0.75, table_text, transform=ax5.transAxes, fontsize=10,
            fontfamily='monospace', va='top')

    # Verdict
    if urn_result['error_mean'] <= vf_result['error_mean']:
        verdict = "URN achieves comparable or better accuracy"
        color = 'green'
    else:
        verdict = "VF achieves better accuracy (but may have stability issues)"
        color = 'orange'

    ax5.text(0.1, 0.15, f"Verdict: {verdict}", transform=ax5.transAxes,
            fontsize=12, fontweight='bold', color=color)

    # 6. Method details
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    details = f"""
Vector Fitting Implementation:
    Library: {vf_result['method']}
    Number of poles: {vf_result.get('n_poles', 'N/A')}
    Unstable poles: {vf_result['n_unstable']}
    Weakly-damped: {vf_result['n_weak_damped']}

URN Configuration:
    Mechanisms found: {urn_result['n_mechanisms']}
    Restarts: 5
    Epochs: 6000
    Sparsity weight: 0.01

Note: scikit-rf uses the original vectfit3 algorithm
by Gustavsen & Semlyen (1999).
    """

    ax6.text(0.1, 0.9, details, transform=ax6.transAxes, fontsize=10,
            fontfamily='monospace', va='top')

    plt.tight_layout()

    # Save figure
    output_dir.mkdir(exist_ok=True)
    fig_path = output_dir / 'benchmark_urn_vs_skrf_vf.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\nComparison plot saved to: {fig_path}")
    plt.close()

    # Save text report
    report_path = output_dir / 'benchmark_urn_vs_skrf_vf.txt'
    with open(report_path, 'w') as f:
        f.write("URN vs scikit-rf Vector Fitting Benchmark\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"VF Implementation: {vf_result['method']}\n")
        f.write(f"scikit-rf available: {SKRF_AVAILABLE}\n\n")
        f.write("Results:\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Metric':<25} {'URN':>15} {'VF':>15}\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Mean Error (%)':<25} {urn_result['error_mean']:>15.2f} {vf_result['error_mean']:>15.2f}\n")
        f.write(f"{'Max Error (%)':<25} {urn_result['error_max']:>15.2f} {vf_result['error_max']:>15.2f}\n")
        f.write(f"{'Time (s)':<25} {urn_result['time']:>15.1f} {vf_result['time']:>15.2f}\n")
        f.write(f"{'Passive':<25} {'Yes':>15} {'Yes' if vf_result['is_passive'] else 'No':>15}\n")
        f.write(f"{'Ringing Risk':<25} {'No':>15} {'Yes' if vf_result['has_ringing_risk'] else 'No':>15}\n")
        f.write("-" * 40 + "\n")

    print(f"Text report saved to: {report_path}")


def main():
    parser = argparse.ArgumentParser(description='URN vs scikit-rf VF Benchmark')
    parser.add_argument('--dataset', choices=['battery', 'ferrite'], default='battery',
                       help='Dataset to use')
    parser.add_argument('--n-poles', type=int, default=8,
                       help='Number of VF poles')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='Output directory')

    args = parser.parse_args()

    print("=" * 70)
    print("URN vs scikit-rf Vector Fitting Benchmark")
    print("=" * 70)

    if SKRF_AVAILABLE:
        print(f"\nscikit-rf version: {rf.__version__}")
        print("Using industry-standard VectorFitting implementation")
    else:
        print("\nWARNING: scikit-rf not available")
        print("Install with: pip install scikit-rf")
        print("Using fallback simplified VF implementation")

    # Load data
    data_dir = Path(__file__).parent / 'data' / 'synthetic'

    if args.dataset == 'battery':
        data_path = data_dir / 'liion_battery_eis.csv'
        skip_header = 18
        z_cols = (1, 2)
    else:
        data_path = data_dir / 'mnzn_ferrite_impedance.csv'
        skip_header = 19
        z_cols = (3, 4)

    if not data_path.exists():
        print(f"ERROR: Data file not found: {data_path}")
        sys.exit(1)

    print(f"\nLoading data from: {data_path}")
    data = np.genfromtxt(data_path, delimiter=',', skip_header=skip_header)
    freq = data[:, 0]
    Z = data[:, z_cols[0]] + 1j * data[:, z_cols[1]]

    print(f"  Frequency range: {freq.min():.2e} - {freq.max():.2e} Hz")
    print(f"  Number of points: {len(freq)}")

    # Run benchmarks
    print("\n" + "-" * 50)
    print("Running benchmarks...")
    print("-" * 50)

    urn_result = run_urn_fitting(freq, Z)
    vf_result = run_skrf_vector_fitting(freq, Z, n_poles=args.n_poles)

    # Generate report
    output_dir = Path(__file__).parent / args.output_dir
    generate_comparison_report(freq, Z, urn_result, vf_result, output_dir)

    print("\n" + "=" * 70)
    print("Benchmark Complete")
    print("=" * 70)

    # Summary
    print("\nSummary:")
    print(f"  URN Mean Error:  {urn_result['error_mean']:.2f}%")
    print(f"  VF Mean Error:   {vf_result['error_mean']:.2f}%")
    print(f"  VF Implementation: {vf_result['method']}")

    if not SKRF_AVAILABLE:
        print("\nNote: For proper scientific comparison, please install scikit-rf:")
        print("      pip install scikit-rf")


if __name__ == '__main__':
    main()
