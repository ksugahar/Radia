"""
Paper Figures Generation for PEEC/PRIMA Implementation

Generates publication-quality figures for:
1. PEEC Loop-Star matrix verification
2. PRIMA model order reduction accuracy
3. Adaptive MSC/dipole error analysis
4. WPT coil impedance characteristics at 85 kHz

Output: PDF figures suitable for journal submission
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

try:
    from peec_matrices import PEECBuilder
    print("Using C++ PEEC implementation")
except ImportError:
    print("C++ PEEC module not available, using Python fallback")
    from radia.peec_matrices import PEECBuilder

from lanczos_reduction import LanczosReducer

# Physical constants
MU_0 = 4 * np.pi * 1e-7
SIGMA_COPPER = 5.8e7

# Matplotlib settings for publication quality
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'


def fig1_peec_matrix_verification(output_dir):
    """
    Figure 1: PEEC Loop-Star Matrix Verification
    - (a) L matrix structure (heatmap)
    - (b) Eigenvalue spectrum of L and P
    - (c) Symmetry error vs matrix size
    """
    print("\n" + "=" * 60)
    print("Figure 1: PEEC Matrix Verification")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # Create PEEC matrices for circular loop
    radius = 0.1
    width = 2e-3
    height = 1e-3
    n_segments = 24

    builder = PEECBuilder()
    builder.create_loop([0, 0, 0], radius, [0, 0, 1], width, height, n_segments, SIGMA_COPPER)
    L, R, P, M_LS = builder.build()

    # (a) L matrix heatmap
    ax = axes[0]
    im = ax.imshow(L * 1e9, cmap='RdBu_r', aspect='equal')
    ax.set_xlabel('Column index')
    ax.set_ylabel('Row index')
    ax.set_title('(a) Inductance matrix L [nH]')
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('[nH]')

    # (b) Eigenvalue spectrum
    ax = axes[1]
    eigvals_L = np.linalg.eigvalsh(L)
    ax.semilogy(range(1, len(eigvals_L)+1), np.sort(eigvals_L)[::-1] * 1e9, 'bo-', label='L matrix', markersize=4)

    if P.size > 0:
        eigvals_P = np.linalg.eigvalsh(P)
        ax2 = ax.twinx()
        ax2.semilogy(range(1, len(eigvals_P)+1), np.sort(eigvals_P)[::-1], 'rs--', label='P matrix', markersize=4)
        ax2.set_ylabel('P eigenvalues [1/F]', color='r')
        ax2.tick_params(axis='y', labelcolor='r')

    ax.set_xlabel('Eigenvalue index')
    ax.set_ylabel('L eigenvalues [nH]', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title('(b) Eigenvalue spectrum')
    ax.grid(True, alpha=0.3)

    # (c) Symmetry and positive definiteness vs N
    ax = axes[2]
    N_values = [8, 12, 16, 20, 24, 32]
    sym_errors = []
    min_eigs = []

    for n in N_values:
        builder = PEECBuilder()
        builder.create_loop([0, 0, 0], radius, [0, 0, 1], width, height, n, SIGMA_COPPER)
        L_n, _, _, _ = builder.build()

        sym_err = np.linalg.norm(L_n - L_n.T) / np.linalg.norm(L_n)
        sym_errors.append(sym_err)

        min_eig = np.min(np.linalg.eigvalsh(L_n))
        min_eigs.append(min_eig)

    ax.semilogy(N_values, sym_errors, 'bo-', label='Symmetry error', markersize=5)
    ax.set_xlabel('Number of segments')
    ax.set_ylabel('Relative error', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title('(c) Matrix properties vs N')
    ax.grid(True, alpha=0.3)

    ax2 = ax.twinx()
    ax2.plot(N_values, np.array(min_eigs) * 1e9, 'rs--', label='Min eigenvalue', markersize=5)
    ax2.set_ylabel('Min eigenvalue [nH]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # Add legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='center right')

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'fig1_peec_matrix_verification.pdf')
    plt.savefig(fig_path)
    plt.savefig(fig_path.replace('.pdf', '.png'))
    print(f"Saved: {fig_path}")
    plt.close()

    return fig_path


def fig2_prima_accuracy(output_dir):
    """
    Figure 2: PRIMA Model Order Reduction Accuracy
    - (a) Frequency response: Full vs PRIMA reduced
    - (b) Error vs Lanczos order
    - (c) Computation time vs model size
    """
    print("\n" + "=" * 60)
    print("Figure 2: PRIMA Accuracy Comparison")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # Create PEEC matrices
    length = 0.1
    width = 2e-3
    height = 0.5e-3
    n_segments = 20

    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, SIGMA_COPPER)
    L, R, P, M_LS = builder.build()

    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R

    n = L.shape[0]
    frequencies = np.logspace(3, 8, 51)  # 1 kHz to 100 MHz

    # (a) Frequency response comparison
    ax = axes[0]

    # Full model
    Z_full_list = []
    for f in frequencies:
        omega = 2 * np.pi * f
        Z = R_full + 1j * omega * L
        try:
            Z_port = Z[0, 0] - Z[0, 1:] @ np.linalg.solve(Z[1:, 1:], Z[1:, 0])
        except Exception:
            Z_port = Z[0, 0]
        Z_full_list.append(Z_port)

    # PRIMA reduced (k=5)
    lanczos = LanczosReducer(tol=1e-10)
    result_5 = lanczos.lanczos_symmetric(L, k=5)
    R_red_5 = result_5.Q.T @ R_full @ result_5.Q

    Z_red_5_list = []
    for f in frequencies:
        omega = 2 * np.pi * f
        Z = R_red_5 + 1j * omega * result_5.T
        try:
            Z_port = Z[0, 0] - Z[0, 1:] @ np.linalg.solve(Z[1:, 1:], Z[1:, 0])
        except Exception:
            Z_port = Z[0, 0]
        Z_red_5_list.append(Z_port)

    # PRIMA reduced (k=10)
    result_10 = lanczos.lanczos_symmetric(L, k=10)
    R_red_10 = result_10.Q.T @ R_full @ result_10.Q

    Z_red_10_list = []
    for f in frequencies:
        omega = 2 * np.pi * f
        Z = R_red_10 + 1j * omega * result_10.T
        try:
            Z_port = Z[0, 0] - Z[0, 1:] @ np.linalg.solve(Z[1:, 1:], Z[1:, 0])
        except Exception:
            Z_port = Z[0, 0]
        Z_red_10_list.append(Z_port)

    ax.loglog(frequencies/1e3, np.abs(Z_full_list)*1e3, 'b-', label=f'Full (N={n})', linewidth=1.5)
    ax.loglog(frequencies/1e3, np.abs(Z_red_5_list)*1e3, 'r--', label='PRIMA (k=5)', linewidth=1.5)
    ax.loglog(frequencies/1e3, np.abs(Z_red_10_list)*1e3, 'g:', label='PRIMA (k=10)', linewidth=1.5)
    ax.axvline(x=85, color='gray', linestyle='-.', alpha=0.7, label='85 kHz')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('|Z| [mOhm]')
    ax.set_title('(a) Impedance magnitude')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3, which='both')

    # (b) Error vs Lanczos order
    ax = axes[1]
    k_values = range(2, min(n, 15)+1)
    errors_85k = []
    errors_1M = []

    omega_85k = 2 * np.pi * 85000
    omega_1M = 2 * np.pi * 1e6

    # Reference (full model)
    Z_full_85k = R_full + 1j * omega_85k * L
    Z_port_full_85k = Z_full_85k[0, 0] - Z_full_85k[0, 1:] @ np.linalg.solve(Z_full_85k[1:, 1:], Z_full_85k[1:, 0])

    Z_full_1M = R_full + 1j * omega_1M * L
    Z_port_full_1M = Z_full_1M[0, 0] - Z_full_1M[0, 1:] @ np.linalg.solve(Z_full_1M[1:, 1:], Z_full_1M[1:, 0])

    for k in k_values:
        result_k = lanczos.lanczos_symmetric(L, k=k)
        R_red_k = result_k.Q.T @ R_full @ result_k.Q

        # 85 kHz
        Z_red_85k = R_red_k + 1j * omega_85k * result_k.T
        try:
            Z_port_red_85k = Z_red_85k[0, 0] - Z_red_85k[0, 1:] @ np.linalg.solve(Z_red_85k[1:, 1:], Z_red_85k[1:, 0])
        except Exception:
            Z_port_red_85k = Z_red_85k[0, 0]
        err_85k = np.abs(Z_port_full_85k - Z_port_red_85k) / np.abs(Z_port_full_85k)
        errors_85k.append(err_85k * 100)

        # 1 MHz
        Z_red_1M = R_red_k + 1j * omega_1M * result_k.T
        try:
            Z_port_red_1M = Z_red_1M[0, 0] - Z_red_1M[0, 1:] @ np.linalg.solve(Z_red_1M[1:, 1:], Z_red_1M[1:, 0])
        except Exception:
            Z_port_red_1M = Z_red_1M[0, 0]
        err_1M = np.abs(Z_port_full_1M - Z_port_red_1M) / np.abs(Z_port_full_1M)
        errors_1M.append(err_1M * 100)

    ax.semilogy(list(k_values), errors_85k, 'bo-', label='85 kHz', markersize=5)
    ax.semilogy(list(k_values), errors_1M, 'rs--', label='1 MHz', markersize=5)
    ax.axhline(y=5, color='gray', linestyle=':', alpha=0.7)
    ax.text(max(k_values)-1, 6, '5%', fontsize=8, color='gray')
    ax.set_xlabel('Lanczos order k')
    ax.set_ylabel('Relative error [%]')
    ax.set_title('(b) Error vs reduction order')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (c) Model size comparison
    ax = axes[2]
    n_full = n * n  # Full matrix elements
    k_values_bar = [5, 10, 15]
    n_prima = [2*k - 1 for k in k_values_bar]  # Tridiagonal has 2k-1 non-zeros

    compression = [100 * (1 - n_p / n_full) for n_p in n_prima]

    bars = ax.bar(range(len(k_values_bar)), compression, color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    ax.set_xticks(range(len(k_values_bar)))
    ax.set_xticklabels([f'k={k}' for k in k_values_bar])
    ax.set_ylabel('Compression ratio [%]')
    ax.set_title(f'(c) Model size reduction (N={n})')
    ax.set_ylim([95, 100])

    for i, (bar, comp) in enumerate(zip(bars, compression)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                f'{comp:.1f}%', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'fig2_prima_accuracy.pdf')
    plt.savefig(fig_path)
    plt.savefig(fig_path.replace('.pdf', '.png'))
    print(f"Saved: {fig_path}")
    plt.close()

    return fig_path


def fig3_adaptive_msc_error(output_dir):
    """
    Figure 3: Adaptive MSC/Dipole Error Analysis
    - (a) Error vs distance (normalized)
    - (b) Threshold selection for target error
    - (c) Speedup potential
    """
    print("\n" + "=" * 60)
    print("Figure 3: Adaptive MSC/Dipole Error")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # Theoretical error model: error ~ (a/r)^3
    r_over_a = np.logspace(-0.5, 1.5, 100)  # 0.3 to 30
    error_theory = (1.0 / r_over_a) ** 3 * 100  # Percent

    # (a) Error vs normalized distance
    ax = axes[0]
    ax.loglog(r_over_a, error_theory, 'b-', linewidth=2, label='Dipole error ~ (a/r)^3')
    ax.axhline(y=5, color='r', linestyle='--', label='5% threshold')
    ax.axhline(y=1, color='g', linestyle=':', label='1% threshold')

    # Mark threshold distances
    r_5pct = 1.0 / (0.05 ** (1/3))  # ~2.7
    r_1pct = 1.0 / (0.01 ** (1/3))  # ~4.6
    ax.axvline(x=r_5pct, color='r', linestyle='--', alpha=0.5)
    ax.axvline(x=r_1pct, color='g', linestyle=':', alpha=0.5)

    ax.fill_between(r_over_a, 0.01, error_theory, where=(r_over_a < r_5pct),
                    alpha=0.3, color='yellow', label='MSC region')
    ax.fill_between(r_over_a, 0.01, error_theory, where=(r_over_a >= r_5pct),
                    alpha=0.3, color='blue', label='Dipole region')

    ax.set_xlabel('Normalized distance r/a')
    ax.set_ylabel('Relative error [%]')
    ax.set_title('(a) Dipole approximation error')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim([0.3, 30])
    ax.set_ylim([0.01, 1000])

    # (b) Threshold selection
    ax = axes[1]
    target_errors = np.logspace(-3, -0.5, 50)  # 0.1% to 30%
    threshold_factors = 1.0 / (target_errors ** (1/3))

    ax.loglog(target_errors * 100, threshold_factors, 'b-', linewidth=2)
    ax.set_xlabel('Target error [%]')
    ax.set_ylabel('Distance threshold r/a')
    ax.set_title('(b) Threshold vs target error')

    # Mark common targets
    for err, col, lbl in [(1, 'g', '1%'), (5, 'r', '5%'), (10, 'orange', '10%')]:
        thresh = 1.0 / ((err/100) ** (1/3))
        ax.plot(err, thresh, 'o', color=col, markersize=8)
        ax.annotate(f'{lbl}: r>{thresh:.1f}a', (err, thresh),
                    textcoords='offset points', xytext=(10, 5), fontsize=8)

    ax.grid(True, alpha=0.3, which='both')

    # (c) Speedup potential (conceptual)
    ax = axes[2]
    # Assume uniform element distribution in cube
    N_elements = np.array([100, 500, 1000, 5000, 10000])
    # Fraction of elements in near-field (r < 3a) for typical geometry
    near_field_fraction = 0.1 + 0.3 * np.exp(-N_elements / 2000)

    # Time ratio: MSC is ~10x slower than dipole
    msc_to_dipole_ratio = 10.0
    time_full_msc = N_elements * msc_to_dipole_ratio
    time_adaptive = N_elements * (near_field_fraction * msc_to_dipole_ratio +
                                   (1 - near_field_fraction) * 1.0)
    speedup = time_full_msc / time_adaptive

    ax.semilogx(N_elements, speedup, 'bo-', markersize=6, linewidth=2)
    ax.set_xlabel('Number of elements')
    ax.set_ylabel('Speedup factor')
    ax.set_title('(c) Adaptive method speedup')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([1, 10])

    for n, s in zip(N_elements, speedup):
        ax.annotate(f'{s:.1f}x', (n, s), textcoords='offset points',
                    xytext=(0, 8), ha='center', fontsize=8)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'fig3_adaptive_msc_error.pdf')
    plt.savefig(fig_path)
    plt.savefig(fig_path.replace('.pdf', '.png'))
    print(f"Saved: {fig_path}")
    plt.close()

    return fig_path


def fig4_wpt_coil_characteristics(output_dir):
    """
    Figure 4: WPT Coil Characteristics at 85 kHz
    - (a) Impedance vs frequency
    - (b) Phase angle vs frequency
    - (c) Quality factor vs coil parameters
    """
    print("\n" + "=" * 60)
    print("Figure 4: WPT Coil at 85 kHz")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # WPT coil parameters
    radius = 0.15  # 15 cm
    width = 3e-3
    height = 1e-3
    n_segments = 36

    builder = PEECBuilder()
    builder.create_loop([0, 0, 0], radius, [0, 0, 1], width, height, n_segments, SIGMA_COPPER)
    L, R, P, M_LS = builder.build()

    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R

    # Frequency sweep
    frequencies = np.logspace(3, 6, 101)  # 1 kHz to 1 MHz

    Z_list = []
    for f in frequencies:
        omega = 2 * np.pi * f
        Z = R_full + 1j * omega * L
        Z_total = np.sum(Z)  # Total loop impedance
        Z_list.append(Z_total)

    Z_array = np.array(Z_list)

    # (a) Impedance magnitude
    ax = axes[0]
    ax.loglog(frequencies/1e3, np.abs(Z_array), 'b-', linewidth=1.5)
    ax.axvline(x=85, color='r', linestyle='--', alpha=0.8, label='85 kHz')

    # Mark 85 kHz point
    idx_85k = np.argmin(np.abs(frequencies - 85000))
    ax.plot(85, np.abs(Z_array[idx_85k]), 'ro', markersize=8)
    ax.annotate(f'|Z|={np.abs(Z_array[idx_85k]):.3f} Ohm',
                (85, np.abs(Z_array[idx_85k])),
                textcoords='offset points', xytext=(10, 10), fontsize=9)

    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('(a) Impedance magnitude')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')

    # (b) Phase angle
    ax = axes[1]
    phase = np.angle(Z_array, deg=True)
    ax.semilogx(frequencies/1e3, phase, 'b-', linewidth=1.5)
    ax.axvline(x=85, color='r', linestyle='--', alpha=0.8, label='85 kHz')
    ax.axhline(y=90, color='gray', linestyle=':', alpha=0.5)

    ax.plot(85, phase[idx_85k], 'ro', markersize=8)
    ax.annotate(f'{phase[idx_85k]:.1f} deg',
                (85, phase[idx_85k]),
                textcoords='offset points', xytext=(10, -15), fontsize=9)

    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('(b) Phase angle')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 95])

    # (c) Quality factor vs radius
    ax = axes[2]
    radii = np.linspace(0.05, 0.25, 9)  # 5 to 25 cm
    Q_factors = []

    for r in radii:
        builder = PEECBuilder()
        builder.create_loop([0, 0, 0], r, [0, 0, 1], width, height, n_segments, SIGMA_COPPER)
        L_r, R_r, _, _ = builder.build()

        if R_r.ndim == 1:
            R_r_full = np.diag(R_r)
        else:
            R_r_full = R_r

        omega_85k = 2 * np.pi * 85000
        Z_r = R_r_full + 1j * omega_85k * L_r
        Z_total_r = np.sum(Z_r)
        Q = np.abs(Z_total_r.imag) / np.abs(Z_total_r.real)
        Q_factors.append(Q)

    ax.plot(radii * 100, Q_factors, 'bo-', markersize=6, linewidth=1.5)
    ax.set_xlabel('Coil radius [cm]')
    ax.set_ylabel('Quality factor Q')
    ax.set_title('(c) Q factor at 85 kHz')
    ax.grid(True, alpha=0.3)

    # Add annotation for typical WPT coil
    ax.axvline(x=15, color='r', linestyle='--', alpha=0.5)
    idx_15cm = np.argmin(np.abs(radii - 0.15))
    ax.annotate(f'15cm: Q={Q_factors[idx_15cm]:.0f}',
                (15, Q_factors[idx_15cm]),
                textcoords='offset points', xytext=(10, 5), fontsize=9)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'fig4_wpt_coil_85khz.pdf')
    plt.savefig(fig_path)
    plt.savefig(fig_path.replace('.pdf', '.png'))
    print(f"Saved: {fig_path}")
    plt.close()

    return fig_path


def main():
    """Generate all paper figures."""
    print("=" * 70)
    print("Paper Figures Generation")
    print("=" * 70)

    output_dir = os.path.dirname(__file__)
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    print(f"\nOutput directory: {figures_dir}")

    # Generate all figures
    figures = []

    try:
        figures.append(fig1_peec_matrix_verification(figures_dir))
    except Exception as e:
        print(f"ERROR in Figure 1: {e}")
        import traceback
        traceback.print_exc()

    try:
        figures.append(fig2_prima_accuracy(figures_dir))
    except Exception as e:
        print(f"ERROR in Figure 2: {e}")
        import traceback
        traceback.print_exc()

    try:
        figures.append(fig3_adaptive_msc_error(figures_dir))
    except Exception as e:
        print(f"ERROR in Figure 3: {e}")
        import traceback
        traceback.print_exc()

    try:
        figures.append(fig4_wpt_coil_characteristics(figures_dir))
    except Exception as e:
        print(f"ERROR in Figure 4: {e}")
        import traceback
        traceback.print_exc()

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"\nGenerated {len(figures)} figures:")
    for fig in figures:
        if fig:
            print(f"  - {os.path.basename(fig)}")

    print(f"\nAll figures saved to: {figures_dir}")

    return len(figures) == 4


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
