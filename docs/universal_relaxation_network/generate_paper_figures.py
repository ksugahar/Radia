#!/usr/bin/env python
"""
Generate publication-quality PNG figures for URN paper.

Requirements:
- Times New Roman font
- 8cm width with 10pt font
- Minimal margins
- Inward ticks on all sides

Output:
- paper/Figures/fig_*.pdf
"""

import sys
import os
from pathlib import Path
import numpy as np

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib
matplotlib.use('Agg')
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic', 'bf':'serif:bold', 'fontset':'cm'})
import matplotlib.pyplot as plt

# Figure size: 8cm width at 600 dpi
CM_TO_INCH = 1 / 2.54
FIG_WIDTH = 8 * CM_TO_INCH  # 8cm = 3.15 inches
FIG_HEIGHT = 6 * CM_TO_INCH  # 6cm (adjustable)
DPI = 600

# Font settings
FONT_NAME = 'times new roman'
FONT_SIZE = 10
FONT_SIZE_TICK = 9
FONT_SIZE_LEGEND = 8


def setup_axis(ax, xlabel, ylabel, title=None):
    """Setup axis with publication-quality formatting."""
    # Axis labels
    ax.set_xlabel(xlabel, fontname=FONT_NAME, fontsize=FONT_SIZE)
    ax.set_ylabel(ylabel, fontname=FONT_NAME, fontsize=FONT_SIZE)
    if title:
        ax.set_title(title, fontname=FONT_NAME, fontsize=FONT_SIZE)

    # Tick labels
    plt.setp(ax.get_xticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK)
    plt.setp(ax.get_yticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK)

    # Tick direction: inward on all sides
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', top=True, right=True, width=0.5, length=4)
    ax.tick_params(which='minor', direction='in', top=True, right=True, width=0.3, length=2)


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

    # 1x2 panel figure: final width is columnwidth (8cm), so each panel is 4cm
    # Font size must be 2x to appear as 10pt when scaled to 50%
    FONT_SIZE_2X = FONT_SIZE * 2  # 20pt -> 10pt after scaling
    FONT_SIZE_TICK_2X = FONT_SIZE_TICK * 2  # 18pt -> 9pt after scaling
    FONT_SIZE_LEGEND_1_2X = FONT_SIZE_LEGEND * 1.2  # smaller legend to avoid overlap

    freq, Z_data = load_nasa_battery_data()
    if freq is None:
        print("  Skipped (data not found)")
        return

    Z_urn = run_urn_fitting(freq, Z_data)
    Z_vf = run_vector_fitting(freq, Z_data, n_poles=10)

    # Create figure with two subplots
    fig = plt.figure(figsize=(FIG_WIDTH * 2, FIG_HEIGHT * 1.5), dpi=DPI)

    # Nyquist plot - 2x sizes for 1x2 panel figure
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(Z_data.real * 1000, -Z_data.imag * 1000, 'ko', markersize=6,
             label='Measured', markerfacecolor='none', markeredgewidth=1.0)
    if Z_urn is not None:
        ax1.plot(Z_urn.real * 1000, -Z_urn.imag * 1000, 'b-', linewidth=2, label='URN')
    ax1.plot(Z_vf.real * 1000, -Z_vf.imag * 1000, 'r--', linewidth=2, label='VF')
    ax1.set_xlabel(r"$Z'$ (m$\Omega$)", fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    ax1.set_ylabel(r"$-Z''$ (m$\Omega$)", fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    ax1.set_title('(a) Nyquist plot', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    plt.setp(ax1.get_xticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
    plt.setp(ax1.get_yticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
    ax1.minorticks_on()
    ax1.tick_params(which='major', direction='in', top=True, right=True, width=1.0, length=8)
    ax1.tick_params(which='minor', direction='in', top=True, right=True, width=0.6, length=4)
    ax1.legend(loc='upper left', frameon=False, prop={'family': FONT_NAME, 'size': FONT_SIZE_LEGEND_1_2X})

    # Bode magnitude plot - 2x sizes for 1x2 panel figure
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.loglog(freq, np.abs(Z_data) * 1000, 'ko', markersize=6,
               label='Measured', markerfacecolor='none', markeredgewidth=1.0)
    if Z_urn is not None:
        ax2.loglog(freq, np.abs(Z_urn) * 1000, 'b-', linewidth=2, label='URN')
    ax2.loglog(freq, np.abs(Z_vf) * 1000, 'r--', linewidth=2, label='VF')
    ax2.set_xlabel('Frequency (Hz)', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    ax2.set_ylabel(r'$|Z|$ (m$\Omega$)', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    ax2.set_title('(b) Bode magnitude', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
    plt.setp(ax2.get_xticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
    plt.setp(ax2.get_yticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
    ax2.minorticks_on()
    ax2.tick_params(which='major', direction='in', top=True, right=True, width=1.0, length=8)
    ax2.tick_params(which='minor', direction='in', top=True, right=True, width=0.6, length=4)
    ax2.legend(loc='upper right', frameon=False, prop={'family': FONT_NAME, 'size': FONT_SIZE_LEGEND_1_2X})

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_nasa_battery_eis.pdf'
    plt.savefig(output_path, dpi=DPI)
    plt.close()
    print(f"  Saved: {output_path}")


def fig_tdk_ferrite_impedance(output_dir):
    """Generate TDK ferrite impedance comparison figure."""
    print("[2/4] Generating TDK ferrite impedance figure...")

    materials = ['pc47', 'pc50', 'pc95', 'pc200']
    # 2x2 panel figure: final width is columnwidth (8cm), so each panel is 4cm
    # Font size must be 2x to appear as 10pt when scaled to 50%
    FONT_SIZE_2X = FONT_SIZE * 2  # 20pt -> 10pt after scaling
    FONT_SIZE_TICK_2X = FONT_SIZE_TICK * 2  # 18pt -> 9pt after scaling
    FONT_SIZE_LEGEND_1_5X = FONT_SIZE_LEGEND * 1.5  # 12pt -> 6pt after scaling (smaller legend)

    fig = plt.figure(figsize=(FIG_WIDTH * 2, FIG_HEIGHT * 2), dpi=DPI)

    for idx, material in enumerate(materials):
        freq, Z_data = load_tdk_ferrite_data(material)
        ax = fig.add_subplot(2, 2, idx + 1)

        if freq is None:
            ax.text(0.5, 0.5, f'{material.upper()}\n(data not found)',
                    ha='center', va='center', transform=ax.transAxes,
                    fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
            continue

        Z_urn = run_urn_fitting(freq, Z_data)
        Z_vf = run_vector_fitting(freq, Z_data, n_poles=12)

        ax.loglog(freq, np.abs(Z_data), 'ko', markersize=4,
                  label='Measured', markerfacecolor='none', markeredgewidth=0.8)
        if Z_urn is not None:
            ax.loglog(freq, np.abs(Z_urn), 'b-', linewidth=1.6, label='URN')
        ax.loglog(freq, np.abs(Z_vf), 'r--', linewidth=1.6, label='VF')

        # Inline axis setup with 2x font sizes for 2x2 panel
        ax.set_xlabel('Frequency (Hz)', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
        ax.set_ylabel(r'$|Z|$ ($\Omega$)', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
        ax.set_title(f'({chr(97+idx)}) TDK {material.upper()}', fontname=FONT_NAME, fontsize=FONT_SIZE_2X)
        plt.setp(ax.get_xticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
        plt.setp(ax.get_yticklabels(), fontname=FONT_NAME, fontsize=FONT_SIZE_TICK_2X)
        ax.minorticks_on()
        ax.tick_params(which='major', direction='in', top=True, right=True, width=1.0, length=8)
        ax.tick_params(which='minor', direction='in', top=True, right=True, width=0.6, length=4)
        if idx == 0:
            ax.legend(loc='upper left', frameon=False, prop={'family': FONT_NAME, 'size': FONT_SIZE_LEGEND_1_5X})

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_tdk_ferrite_impedance.pdf'
    plt.savefig(output_path, dpi=DPI)
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

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=DPI)
    ax = fig.add_subplot(1, 1, 1)

    bars1 = ax.bar(x - width/2, vf_nrmse, width, label='Vector Fitting',
                   color='#d62728', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, urn_nrmse, width, label='URN',
                   color='#1f77b4', edgecolor='black', linewidth=0.5)

    setup_axis(ax, '', 'NRMSE')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontname=FONT_NAME, fontsize=FONT_SIZE_TICK)
    ax.legend(loc='upper right', frameon=False, prop={'family': FONT_NAME, 'size': FONT_SIZE_LEGEND})
    ax.set_ylim(0, 0.35)

    # Add improvement annotations
    improvements = [9.1, 39.4, 65.9, -48.9, 48.4]
    for i, (v, u, imp) in enumerate(zip(vf_nrmse, urn_nrmse, improvements)):
        y_pos = max(v, u) + 0.01
        color = 'green' if imp > 0 else 'red'
        text = f'+{imp:.0f}%' if imp > 0 else f'{imp:.0f}%'
        ax.annotate(text, (x[i], y_pos), ha='center', fontname=FONT_NAME, fontsize=7, color=color)

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_urn_vs_vf_comparison.pdf'
    plt.savefig(output_path, dpi=DPI)
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

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=DPI)
    ax = fig.add_subplot(1, 1, 1)

    bars1 = ax.bar(x - width/2, without_attn, width, label='Without Attention',
                   color='#ff7f0e', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, with_attn, width, label='With Attention',
                   color='#2ca02c', edgecolor='black', linewidth=0.5)

    setup_axis(ax, '', 'NRMSE')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontname=FONT_NAME, fontsize=FONT_SIZE_TICK)
    ax.legend(loc='upper right', frameon=False, prop={'family': FONT_NAME, 'size': FONT_SIZE_LEGEND})
    ax.set_ylim(0, 0.7)

    # Add improvement annotations
    improvements = [59.0, 80.1, 81.7, 77.8]
    for i, (wo, w, imp) in enumerate(zip(without_attn, with_attn, improvements)):
        y_pos = wo + 0.02
        ax.annotate(f'-{imp:.0f}%', (x[i], y_pos), ha='center', fontname=FONT_NAME, fontsize=7, color='green')

    plt.tight_layout(pad=0.5)

    output_path = output_dir / 'fig_ablation_attention.pdf'
    plt.savefig(output_path, dpi=DPI)
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("Generating Publication-Quality Figures for URN Paper")
    print("=" * 60)
    print(f"Figure width: 8cm, Font: {FONT_NAME} {FONT_SIZE}pt")
    print(f"DPI: {DPI}, Ticks: inward on all sides")
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
