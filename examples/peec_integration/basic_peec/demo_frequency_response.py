"""
Demonstration of Frequency Response Analysis

This script demonstrates frequency response analysis for both:
1. PEEC conductors (impedance Z(f) with skin effect)
2. MMM magnetic materials (complex permeability mu(f) with H-dependence)

The examples show:
- DC to 10 MHz frequency sweep
- Dowell skin effect correction for conductors
- H-dependent permeability for nonlinear magnetic materials
- Comparison of different operating H levels

Author: Radia Development Team
Date: 2026-01-15
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add build path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

# Import analysis framework
from analysis import (
    UnifiedAnalysis,
    UnifiedMMMAnalysis,
    MU_0,
)


def create_coil_peec_model():
    """
    Create PEEC matrices for a multi-turn coil.

    Returns:
        L: Inductance matrix [H]
        R: Resistance matrix [Ohm]
    """
    # Coil parameters
    n_turns = 20
    r_coil = 0.03       # Coil radius [m]
    wire_length = 2.0 * np.pi * r_coil * n_turns
    wire_radius = 0.3e-3  # Wire radius [m]
    sigma = 5.8e7       # Copper conductivity [S/m]

    # Discretize into segments
    n_segments = 30

    # Resistance per segment
    segment_length = wire_length / n_segments
    segment_area = np.pi * wire_radius**2
    R_segment = segment_length / (sigma * segment_area)
    R = np.ones(n_segments) * R_segment

    # Self-inductance per segment (Neumann formula)
    mu_0 = 4.0 * np.pi * 1e-7
    L_self = (mu_0 * segment_length / (2.0 * np.pi)) * (
        np.log(2.0 * segment_length / wire_radius) - 1.0
    )

    # Build inductance matrix with mutual inductance
    L = np.zeros((n_segments, n_segments))

    for i in range(n_segments):
        L[i, i] = L_self

        for j in range(i + 1, n_segments):
            # Approximate distance based on coil geometry
            angle_diff = 2.0 * np.pi * abs(i - j) / n_segments
            distance = 2.0 * r_coil * np.sin(angle_diff / 2.0)
            distance = max(distance, wire_radius)

            # Neumann formula (simplified)
            M_ij = (mu_0 * segment_length**2) / (4.0 * np.pi * distance)
            L[i, j] = M_ij
            L[j, i] = M_ij

    return L, R, wire_radius * 2


def demo_peec_frequency_response():
    """Demonstrate PEEC frequency response with skin effect."""
    print("\n" + "="*70)
    print("PEEC Conductor Frequency Response")
    print("="*70)

    # Create PEEC model
    L, R, wire_diameter = create_coil_peec_model()
    print(f"Coil segments: {L.shape[0]}")
    print(f"Wire diameter: {wire_diameter*1e3:.2f} mm")
    print(f"DC Resistance: {np.sum(R):.4f} Ohm")
    print(f"DC Inductance: {np.sum(L)*1e6:.3f} uH")

    # Create analysis object
    analysis = UnifiedAnalysis()
    analysis.set_peec_model(L, R, prima_order=5)
    analysis.set_skin_effect(
        sigma=5.8e7,
        use_dowell=True,
        conductor_height=wire_diameter
    )

    # Frequency sweep: 1 Hz to 10 MHz
    frequencies = np.logspace(0, 7, 200)

    result = analysis.frequency_sweep(frequencies)

    print(f"\nFrequency sweep: {frequencies[0]:.0f} Hz to {frequencies[-1]/1e6:.1f} MHz")
    print(f"Computation time: {result.computation_time*1000:.2f} ms")

    # Print key values
    for f_target in [100, 1000, 10000, 100000, 1000000]:
        idx = np.argmin(np.abs(frequencies - f_target))
        f = frequencies[idx]
        R_ac = result.resistance[idx]
        L_ac = result.inductance[idx] * 1e6
        Q = result.quality_factor[idx]
        print(f"  @ {f/1e3:6.1f} kHz: R = {R_ac:.4f} Ohm, L = {L_ac:.3f} uH, Q = {Q:.1f}")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Impedance magnitude
    ax = axes[0, 0]
    Z_mag = np.abs(result.impedance)
    ax.loglog(frequencies, Z_mag)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Impedance Magnitude')
    ax.grid(True, which='both', alpha=0.3)

    # AC Resistance (normalized to DC)
    ax = axes[0, 1]
    R_dc = np.sum(R)
    ax.semilogx(frequencies, result.resistance / R_dc)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('R_ac / R_dc')
    ax.set_title('AC Resistance Ratio (Skin Effect)')
    ax.grid(True, which='both', alpha=0.3)
    ax.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='DC')
    ax.legend()

    # Inductance
    ax = axes[1, 0]
    L_dc = np.sum(L) * 1e6
    ax.semilogx(frequencies[1:], result.inductance[1:] * 1e6)
    ax.axhline(y=L_dc, color='r', linestyle='--', alpha=0.5, label=f'DC = {L_dc:.2f} uH')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Inductance')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()

    # Quality factor
    ax = axes[1, 1]
    ax.semilogx(frequencies, result.quality_factor)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Q')
    ax.set_title('Quality Factor')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_peec_frequency_response.png', dpi=150)
    print("\nSaved: demo_peec_frequency_response.png")

    return result


def demo_mmm_frequency_response():
    """Demonstrate MMM frequency response with H-dependent permeability."""
    print("\n" + "="*70)
    print("MMM Magnetic Material Frequency Response")
    print("="*70)

    # Create MMM analysis
    analysis = UnifiedMMMAnalysis()

    # For this demo, we don't use a Radia object but directly set parameters
    # In real usage, you would create a Radia object first
    analysis._solver.build = lambda: True
    analysis._solver._is_built = True
    analysis._is_configured = True
    analysis._mode = 'direct'

    # Material parameters (soft iron)
    analysis._solver.set_material_parameters(
        mu_r_static=5000,
        mu_r_inf=1.0,
        tau=1e-6,  # 1 us relaxation time
        sigma=1e7  # Iron conductivity
    )

    # H-dependent permeability (saturation curve)
    mu_H_table = np.array([
        [0, 5000],       # Initial permeability
        [50, 4800],
        [100, 4500],
        [200, 4000],
        [500, 3000],
        [1000, 2000],
        [2000, 1000],
        [5000, 300],
        [10000, 100],    # Deep saturation
        [50000, 20],
    ])

    analysis.set_h_dependent_permeability(mu_H_table, tau=1e-6, sigma=1e7)

    # Frequency sweep
    frequencies = np.logspace(0, 7, 200)

    # Compare different operating H levels
    H_levels = [0, 100, 1000, 5000]
    results = {}

    print("\nFrequency response at different H levels:")
    for H in H_levels:
        result = analysis.frequency_sweep(frequencies, H_operating=H)
        results[H] = result

        # Print 1 kHz values
        idx_1k = np.argmin(np.abs(frequencies - 1e3))
        mu_r_real = result.mu_real[idx_1k] / MU_0
        mu_r_imag = result.mu_imag[idx_1k] / MU_0
        print(f"  H = {H:5d} A/m: mu_r' = {mu_r_real:.0f}, mu_r\" = {mu_r_imag:.0f}")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # mu' vs frequency
    ax = axes[0, 0]
    for H, result in results.items():
        ax.semilogx(frequencies, result.mu_real / MU_0, label=f'H = {H} A/m')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("mu_r' (real part)")
    ax.set_title('Real Permeability')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()

    # mu" vs frequency
    ax = axes[0, 1]
    for H, result in results.items():
        ax.semilogx(frequencies, result.mu_imag / MU_0, label=f'H = {H} A/m')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("mu_r\" (imaginary part)")
    ax.set_title('Imaginary Permeability (Loss)')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()

    # Loss tangent vs frequency
    ax = axes[1, 0]
    for H, result in results.items():
        ax.semilogx(frequencies, result.loss_tangent, label=f'H = {H} A/m')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('tan(delta)')
    ax.set_title('Loss Tangent')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()

    # |mu| vs frequency
    ax = axes[1, 1]
    for H, result in results.items():
        mu_mag = np.abs(result.complex_permeability) / MU_0
        ax.loglog(frequencies, mu_mag, label=f'H = {H} A/m')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('|mu_r|')
    ax.set_title('Permeability Magnitude')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig('demo_mmm_frequency_response.png', dpi=150)
    print("\nSaved: demo_mmm_frequency_response.png")

    # Additional plot: mu_r(H) curve at different frequencies
    fig, ax = plt.subplots(figsize=(8, 5))

    H_sweep = np.linspace(0, 10000, 100)
    freq_targets = [100, 1000, 10000, 100000]

    for f_target in freq_targets:
        mu_r_values = []
        for H in H_sweep:
            analysis.set_operating_field(H)
            result = analysis._solver.solve_frequency(np.array([f_target]))
            mu_r_values.append(np.real(result.complex_permeability[0]) / MU_0)

        ax.plot(H_sweep, mu_r_values, label=f'f = {f_target/1e3:.0f} kHz')

    ax.set_xlabel('H [A/m]')
    ax.set_ylabel("mu_r'")
    ax.set_title('Permeability vs H at Different Frequencies')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig('demo_mmm_permeability_vs_H.png', dpi=150)
    print("Saved: demo_mmm_permeability_vs_H.png")

    return results


def main():
    """Main demonstration function."""
    print("="*70)
    print("Frequency Response Analysis Demo")
    print("="*70)

    # PEEC conductor frequency response
    demo_peec_frequency_response()

    # MMM magnetic material frequency response
    demo_mmm_frequency_response()

    print("\n" + "="*70)
    print("Demo completed successfully!")
    print("="*70)

    plt.show()


if __name__ == "__main__":
    main()
