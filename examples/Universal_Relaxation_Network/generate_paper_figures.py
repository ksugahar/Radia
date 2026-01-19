#!/usr/bin/env python
"""
Generate publication-quality PDF figures for URN paper.

Requirements:
- Times New Roman font
- 8cm width with 10pt font
- Minimal margins
- Inward ticks on all sides

Output:
- docs/paper/figures/fig_nasa_battery_eis.pdf
- docs/paper/figures/fig_tdk_ferrite_impedance.pdf
- docs/paper/figures/fig_urn_vs_vf_comparison.pdf
- docs/paper/figures/fig_ablation_attention.pdf
"""

import sys
import os
from pathlib import Path
import numpy as np

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Publication settings
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 10
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 8
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.02

# Tick settings: inward on all sides
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['xtick.bottom'] = True
rcParams['ytick.left'] = True
rcParams['ytick.right'] = True
rcParams['axes.linewidth'] = 0.5
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['xtick.minor.width'] = 0.3
rcParams['ytick.minor.width'] = 0.3

# Figure size: 8cm width
CM_TO_INCH = 1 / 2.54
FIG_WIDTH = 8 * CM_TO_INCH  # 8cm
FIG_HEIGHT = 6 * CM_TO_INCH  # 6cm (adjustable)


def load_nasa_battery_data():
    """Load NASA 18650 battery EIS data."""
    data_path = Path(__file__).parent / 'data' / 'real_world' / 'nasa_battery' / 'nasa_18650_eis.csv'
    if not data_path.exists():
        print(f"[WARN] NASA data not found: {data_path}")
        return None, None
    data = np.loadtxt(data_path, delimiter=',', skiprows=24)
    freq = data[:, 0]
    Z = data[:, 1] + 1j * data[:, 2]
    return freq, Z


def load_tdk_ferrite_data(material='pc50'):
    """Load TDK ferrite impedance data."""
    data_path = Path(__file__).parent / 'data' / 'real_world' / 'tdk_ferrite' / f'tdk_{material}_impedance.csv'
    if not data_path.exists():
        print(f"[WARN] TDK data not found: {data_path}")
        return None, None
    data = np.loadtxt(data_path, delimiter=',', skiprows=19)
    freq = data[:, 0]
    Z = data[:, 1] + 1j * data[:, 2]
    return freq, Z


def run_urn_fitting(freq, Z_data):
    """Run URN fitting and return fitted impedance."""
    try:
        from universal_relaxation_network import URNConfig, train_urn
        import torch

        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_warburg=2,
            sparsity_weight=0.01,
            lr=0.02,
            n_epochs=3000,
            n_restarts=2
        )

        model = train_urn(freq, Z_data, config)
        omega = 2 * np.pi * freq
        omega_tensor = torch.tensor(omega, dtype=torch.float64)
        Z_fit = model(omega_tensor).detach().numpy()
        return Z_fit
    except Exception as e:
        print(f"[WARN] URN fitting failed: {e}")
        return None


def run_vector_fitting(freq, Z_data, n_poles=10):
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

    # Build basis matrix and solve
    A = np.zeros((len(s), len(poles) + 1), dtype=complex)
    A[:, 0] = 1
    for i, p in enumerate(poles):
        A[:, i + 1] = 1 / (s - p)

    coeffs, _, _, _ = np.linalg.lstsq(A, Z_data, rcond=None)
    Z_fit = A @ coeffs
    return Z_fit


def create_output_dir():
    """Create output directory for figures."""
    output_dir = Path(__file__).parent / 'paper' / 'Figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def fig_nasa_battery_eis(output_dir):
    """Generate NASA battery EIS comparison figure."""
    print("[1/4] Generating NASA battery EIS figure...")

    freq, Z_data = load_nasa_battery_data()
    if freq is None:
        print("  Skipped (data not found)")
        return

    Z_urn = run_urn_fitting(freq, Z_data)
    Z_vf = run_vector_fitting(freq, Z_data, n_poles=10)

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH * 2, FIG_HEIGHT))

    # Nyquist plot
    ax1 = axes[0]
    ax1.plot(Z_data.real * 1000, -Z_data.imag * 1000, 'ko', markersize=3,
             label='Measured', markerfacecolor='none', markeredgewidth=0.5)
    if Z_urn is not None:
        ax1.plot(Z_urn.real * 1000, -Z_urn.imag * 1000, 'b-', linewidth=1, label='URN')
    ax1.plot(Z_vf.real * 1000, -Z_vf.imag * 1000, 'r--', linewidth=1, label='VF')
    ax1.set_xlabel(r"$Z'$ (m$\Omega$)")
    ax1.set_ylabel(r"$-Z''$ (m$\Omega$)")
    ax1.legend(loc='upper right', frameon=False)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title('(a) Nyquist plot')

    # Bode magnitude plot
    ax2 = axes[1]
    ax2.loglog(freq, np.abs(Z_data) * 1000, 'ko', markersize=3,
               label='Measured', markerfacecolor='none', markeredgewidth=0.5)
    if Z_urn is not None:
        ax2.loglog(freq, np.abs(Z_urn) * 1000, 'b-', linewidth=1, label='URN')
    ax2.loglog(freq, np.abs(Z_vf) * 1000, 'r--', linewidth=1, label='VF')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel(r'$|Z|$ (m$\Omega$)')
    ax2.legend(loc='upper right', frameon=False)
    ax2.set_title('(b) Bode magnitude')

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_nasa_battery_eis.pdf'
    plt.savefig(output_path, format='pdf')
    plt.close()
    print(f"  Saved: {output_path}")


def fig_tdk_ferrite_impedance(output_dir):
    """Generate TDK ferrite impedance comparison figure."""
    print("[2/4] Generating TDK ferrite impedance figure...")

    materials = ['pc47', 'pc50', 'pc95', 'pc200']
    fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH * 2, FIG_HEIGHT * 2))
    axes = axes.flatten()

    for idx, material in enumerate(materials):
        freq, Z_data = load_tdk_ferrite_data(material)
        if freq is None:
            axes[idx].text(0.5, 0.5, f'{material.upper()}\n(data not found)',
                          ha='center', va='center', transform=axes[idx].transAxes)
            continue

        Z_urn = run_urn_fitting(freq, Z_data)
        Z_vf = run_vector_fitting(freq, Z_data, n_poles=12)

        ax = axes[idx]
        ax.loglog(freq, np.abs(Z_data), 'ko', markersize=2,
                  label='Measured', markerfacecolor='none', markeredgewidth=0.4)
        if Z_urn is not None:
            ax.loglog(freq, np.abs(Z_urn), 'b-', linewidth=0.8, label='URN')
        ax.loglog(freq, np.abs(Z_vf), 'r--', linewidth=0.8, label='VF')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(r'$|Z|$ ($\Omega$)')
        ax.set_title(f'({chr(97+idx)}) TDK {material.upper()}')
        if idx == 0:
            ax.legend(loc='upper left', frameon=False)

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_tdk_ferrite_impedance.pdf'
    plt.savefig(output_path, format='pdf')
    plt.close()
    print(f"  Saved: {output_path}")


def fig_urn_vs_vf_comparison(output_dir):
    """Generate URN vs VF bar chart comparison."""
    print("[3/4] Generating URN vs VF comparison figure...")

    # Data from validation results
    datasets = ['NASA\nBattery', 'PC47', 'PC50', 'PC95', 'PC200']
    vf_nrmse = [0.2700, 0.0146, 0.0288, 0.0080, 0.0108]
    urn_nrmse = [0.2454, 0.0088, 0.0098, 0.0120, 0.0056]

    x = np.arange(len(datasets))
    width = 0.35

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    bars1 = ax.bar(x - width/2, vf_nrmse, width, label='Vector Fitting',
                   color='#d62728', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, urn_nrmse, width, label='URN',
                   color='#1f77b4', edgecolor='black', linewidth=0.5)

    ax.set_ylabel('NRMSE')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.legend(loc='upper right', frameon=False)
    ax.set_ylim(0, 0.35)

    # Add improvement annotations
    improvements = [9.1, 39.4, 65.9, -48.9, 48.4]
    for i, (v, u, imp) in enumerate(zip(vf_nrmse, urn_nrmse, improvements)):
        y_pos = max(v, u) + 0.01
        if imp > 0:
            ax.annotate(f'+{imp:.0f}%', (x[i], y_pos), ha='center', fontsize=7, color='green')
        else:
            ax.annotate(f'{imp:.0f}%', (x[i], y_pos), ha='center', fontsize=7, color='red')

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_urn_vs_vf_comparison.pdf'
    plt.savefig(output_path, format='pdf')
    plt.close()
    print(f"  Saved: {output_path}")


def fig_ablation_attention(output_dir):
    """Generate ablation study figure for attention mechanism."""
    print("[4/4] Generating ablation study figure...")

    # Data from ablation study
    datasets = ['NASA\nBattery', 'PC47', 'PC50', 'PC200']
    without_attn = [0.5982, 0.0443, 0.0535, 0.0252]
    with_attn = [0.2454, 0.0088, 0.0098, 0.0056]

    x = np.arange(len(datasets))
    width = 0.35

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    bars1 = ax.bar(x - width/2, without_attn, width, label='Without Attention',
                   color='#ff7f0e', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, with_attn, width, label='With Attention',
                   color='#2ca02c', edgecolor='black', linewidth=0.5)

    ax.set_ylabel('NRMSE')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.legend(loc='upper right', frameon=False)
    ax.set_ylim(0, 0.7)

    # Add improvement annotations
    improvements = [59.0, 80.1, 81.7, 77.8]
    for i, (wo, w, imp) in enumerate(zip(without_attn, with_attn, improvements)):
        y_pos = wo + 0.02
        ax.annotate(f'-{imp:.0f}%', (x[i], y_pos), ha='center', fontsize=7, color='green')

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_ablation_attention.pdf'
    plt.savefig(output_path, format='pdf')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("Generating Publication-Quality Figures for URN Paper")
    print("=" * 60)
    print(f"Figure width: 8cm, Font: Times New Roman 10pt")
    print(f"Ticks: inward on all sides")
    print()

    output_dir = create_output_dir()
    print(f"Output directory: {output_dir}")
    print()

    fig_nasa_battery_eis(output_dir)
    fig_tdk_ferrite_impedance(output_dir)
    fig_urn_vs_vf_comparison(output_dir)
    fig_ablation_attention(output_dir)

    print()
    print("=" * 60)
    print("All figures generated successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
