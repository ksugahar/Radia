"""
Demonstration of Loop-Star Decomposition for PEEC Analysis

This script demonstrates the Loop-Star decomposition technique for
frequency response analysis in PEEC (Partial Element Equivalent Circuit).

Loop-Star decomposition separates the current distribution into:
- Loop (solenoidal) currents: Circulating currents, inductive behavior
- Star (irrotational) currents: Divergent currents, capacitive behavior

Benefits:
1. Low-frequency numerical stability (avoids EFIE breakdown)
2. Physical insight into inductive vs capacitive behavior
3. MQS approximation for power electronics (ignore Star currents)

The system equation becomes:
    [Z_LL  Z_LS] [I_L]   [V_L]
    [Z_SL  Z_SS] [I_S] = [V_S]

where:
    Z_LL = R_LL + j*omega*L_LL (Loop-Loop: inductive)
    Z_SS = R_SS + j*omega*L_SS + P_SS/(j*omega) (Star-Star: capacitive)

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
from analysis import UnifiedAnalysis, MU_0


def create_simple_loop_star_model():
    """
    Create a simple Loop-Star PEEC model for a parallel LC circuit.

    This represents an inductor in parallel with parasitic capacitance,
    which is common in real inductors at high frequencies.

    Physical model:
    - 1 Loop (inductance + resistance)
    - 1 Star (parasitic capacitance)
    - Zero coupling (parallel connection: both see same voltage)

    Returns:
        L_LL: Loop-Loop inductance matrix [H]
        R_LL: Loop-Loop resistance matrix [Ohm]
        L_LS: Loop-Star coupling [H] (zero for parallel)
        L_SS: Star-Star inductance [H] (zero)
        P_SS: Star-Star potential coefficient [1/F]
    """
    # Physical parameters
    L = 10e-6          # Inductance [H] (10 uH)
    R = 0.1            # Resistance [Ohm]
    C = 10e-12         # Parasitic capacitance [F] (10 pF)

    # Expected resonant frequency: f = 1/(2*pi*sqrt(L*C)) = 15.92 MHz

    # Loop-Loop inductance (1x1)
    L_LL = np.array([[L]])

    # Loop-Loop resistance (1x1)
    R_LL = np.array([[R]])

    # Loop-Star coupling: zero for parallel connection
    # (both branches see the same terminal voltage)
    L_LS = np.array([[0.0]])

    # Star-Star inductance: zero (pure capacitance)
    L_SS = np.array([[0.0]])

    # Star-Star potential coefficient P = 1/C
    P_SS = np.array([[1.0 / C]])

    return L_LL, R_LL, L_LS, L_SS, P_SS


def create_multi_turn_coil_loop_star(n_turns=5):
    """
    Create Loop-Star model for a multi-turn coil.

    The coil has:
    - n_turns independent loops (one per turn)
    - (n_turns-1) star DOFs (capacitance between adjacent turns)

    Parameters:
        n_turns: Number of coil turns

    Returns:
        L_LL, R_LL, L_LS, L_SS, P_SS matrices
    """
    # Physical parameters per turn
    L_turn = 2e-6       # Inductance per turn [H]
    R_turn = 0.05       # Resistance per turn [Ohm]
    C_inter = 5e-12     # Inter-turn capacitance [F]

    n_loops = n_turns
    n_stars = n_turns - 1

    # Loop-Loop inductance (includes mutual coupling)
    L_LL = np.zeros((n_loops, n_loops))
    for i in range(n_loops):
        L_LL[i, i] = L_turn
        for j in range(i + 1, n_loops):
            # Mutual inductance decreases with distance
            k = 0.8 * np.exp(-0.5 * abs(i - j))  # Coupling coefficient
            M_ij = k * L_turn
            L_LL[i, j] = M_ij
            L_LL[j, i] = M_ij

    # Loop-Loop resistance (diagonal)
    R_LL = np.diag([R_turn] * n_loops)

    # Loop-Star coupling (n_loops x n_stars)
    # Each star connects adjacent loops
    L_LS = np.zeros((n_loops, n_stars))
    for s in range(n_stars):
        # Star s connects loop s and loop s+1
        L_LS[s, s] = 0.01e-6      # Small coupling
        L_LS[s + 1, s] = -0.01e-6  # Opposite sign

    # Star-Star inductance (small diagonal)
    L_SS = np.eye(n_stars) * 0.001e-6

    # Star-Star potential coefficient (capacitance)
    P_SS = np.eye(n_stars) / C_inter

    return L_LL, R_LL, L_LS, L_SS, P_SS


def demo_simple_loop_star():
    """Demonstrate Loop-Star analysis on a parallel LC circuit."""
    print("\n" + "=" * 70)
    print("1. Simple Loop-Star Model (Parallel LC Circuit)")
    print("=" * 70)

    # Create model
    L_LL, R_LL, L_LS, L_SS, P_SS = create_simple_loop_star_model()

    L = L_LL[0, 0]
    C = 1.0 / P_SS[0, 0]
    f_res_expected = 1.0 / (2.0 * np.pi * np.sqrt(L * C))

    print(f"Loop DOFs: {L_LL.shape[0]}")
    print(f"Star DOFs: {P_SS.shape[0]}")
    print(f"L = {L*1e6:.2f} uH")
    print(f"R = {R_LL[0,0]:.3f} Ohm")
    print(f"C = {C*1e12:.2f} pF")
    print(f"Expected resonant frequency: {f_res_expected/1e6:.2f} MHz")

    # Create analysis object
    analysis = UnifiedAnalysis()
    analysis.set_loop_star_model(L_LL, R_LL, L_LS, L_SS, P_SS=P_SS)

    # Frequency sweep: 1 Hz to 100 MHz
    frequencies = np.logspace(0, 8, 300)

    # Compare MQS mode (Loop-only) vs Full Loop-Star
    print("\n--- MQS Mode (Loop-only) ---")
    analysis.set_mqs_mode(True)
    result_mqs = analysis.frequency_sweep_loop_star(frequencies)
    print(f"Computation time: {result_mqs.computation_time*1000:.2f} ms")

    print("\n--- Full Loop-Star Mode ---")
    analysis.set_mqs_mode(False)
    result_full = analysis.frequency_sweep_loop_star(frequencies)
    print(f"Computation time: {result_full.computation_time*1000:.2f} ms")

    # Print some key values
    print("\nImpedance comparison at selected frequencies:")
    for f_target in [1e3, 100e3, 1e6, 10e6, 100e6]:
        idx = np.argmin(np.abs(frequencies - f_target))
        f = frequencies[idx]
        Z_mqs = result_mqs.impedance[idx]
        Z_full = result_full.impedance[idx]
        diff_pct = abs(Z_mqs - Z_full) / abs(Z_full) * 100 if abs(Z_full) > 0 else 0
        print(f"  @ {f/1e6:.1f} MHz: |Z_MQS| = {abs(Z_mqs):.2f} Ohm, "
              f"|Z_full| = {abs(Z_full):.2f} Ohm, diff = {diff_pct:.1f}%")

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Impedance magnitude
    ax = axes[0, 0]
    ax.loglog(frequencies, np.abs(result_mqs.impedance), 'b-', label='MQS (Loop-only)', alpha=0.8)
    ax.loglog(frequencies, np.abs(result_full.impedance), 'r--', label='Full Loop-Star', alpha=0.8)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Impedance Magnitude')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Impedance phase
    ax = axes[0, 1]
    phase_mqs = np.angle(result_mqs.impedance) * 180 / np.pi
    phase_full = np.angle(result_full.impedance) * 180 / np.pi
    ax.semilogx(frequencies, phase_mqs, 'b-', label='MQS', alpha=0.8)
    ax.semilogx(frequencies, phase_full, 'r--', label='Full', alpha=0.8)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Impedance Phase')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.axhline(y=90, color='k', linestyle=':', alpha=0.3)
    ax.axhline(y=-90, color='k', linestyle=':', alpha=0.3)

    # Inductance
    ax = axes[1, 0]
    ax.semilogx(frequencies[1:], result_mqs.inductance[1:] * 1e6, 'b-', label='MQS', alpha=0.8)
    ax.semilogx(frequencies[1:], result_full.inductance[1:] * 1e6, 'r--', label='Full', alpha=0.8)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Effective Inductance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Resistance
    ax = axes[1, 1]
    ax.semilogx(frequencies, result_mqs.resistance, 'b-', label='MQS', alpha=0.8)
    ax.semilogx(frequencies, result_full.resistance, 'r--', label='Full', alpha=0.8)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('R [Ohm]')
    ax.set_title('Effective Resistance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_loop_star_simple.png', dpi=150)
    print("\nSaved: demo_loop_star_simple.png")

    return result_mqs, result_full


def demo_multi_turn_coil():
    """Demonstrate Loop-Star analysis on a multi-turn coil."""
    print("\n" + "=" * 70)
    print("2. Multi-Turn Coil Loop-Star Model")
    print("=" * 70)

    # Create 10-turn coil model
    n_turns = 10
    L_LL, R_LL, L_LS, L_SS, P_SS = create_multi_turn_coil_loop_star(n_turns)

    print(f"Coil turns: {n_turns}")
    print(f"Loop DOFs: {L_LL.shape[0]}")
    print(f"Star DOFs: {P_SS.shape[0]}")
    print(f"Total inductance (L_LL): {np.sum(L_LL)*1e6:.2f} uH")
    print(f"Total resistance (R_LL): {np.sum(R_LL):.4f} Ohm")

    # Create analysis object
    analysis = UnifiedAnalysis()
    analysis.set_loop_star_model(L_LL, R_LL, L_LS, L_SS, P_SS=P_SS)

    # Frequency sweep
    frequencies = np.logspace(0, 9, 400)

    # MQS mode (Loop-only)
    analysis.set_mqs_mode(True)
    result_mqs = analysis.frequency_sweep_loop_star(frequencies)

    # Full Loop-Star
    analysis.set_mqs_mode(False)
    result_full = analysis.frequency_sweep_loop_star(frequencies)

    # Find self-resonant frequency (SRF) where phase crosses zero
    phase_full = np.angle(result_full.impedance) * 180 / np.pi
    srf_idx = None
    for i in range(1, len(phase_full)):
        if phase_full[i-1] > 0 and phase_full[i] < 0:
            srf_idx = i
            break

    if srf_idx is not None:
        srf = frequencies[srf_idx]
        print(f"\nSelf-Resonant Frequency (SRF): {srf/1e6:.2f} MHz")
        print(f"  |Z| at SRF: {np.abs(result_full.impedance[srf_idx]):.1f} Ohm")
    else:
        print("\nSRF not found in frequency range")

    # Print frequency response at key points
    print("\nFrequency response:")
    for f_target in [1e3, 10e3, 100e3, 1e6, 10e6]:
        idx = np.argmin(np.abs(frequencies - f_target))
        f = frequencies[idx]
        L_mqs = result_mqs.inductance[idx] * 1e6
        L_full = result_full.inductance[idx] * 1e6
        Q_mqs = result_mqs.quality_factor[idx]
        Q_full = result_full.quality_factor[idx]
        print(f"  @ {f/1e3:6.0f} kHz: L_MQS = {L_mqs:.2f} uH, L_full = {L_full:.2f} uH, "
              f"Q_MQS = {Q_mqs:.1f}, Q_full = {Q_full:.1f}")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Impedance magnitude (Bode plot)
    ax = axes[0, 0]
    ax.loglog(frequencies, np.abs(result_mqs.impedance), 'b-', label='MQS (inductive only)')
    ax.loglog(frequencies, np.abs(result_full.impedance), 'r-', label='Full (with capacitance)')
    if srf_idx is not None:
        ax.axvline(x=srf, color='g', linestyle='--', alpha=0.5, label=f'SRF = {srf/1e6:.1f} MHz')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title(f'{n_turns}-Turn Coil Impedance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Phase
    ax = axes[0, 1]
    ax.semilogx(frequencies, np.angle(result_mqs.impedance) * 180 / np.pi, 'b-', label='MQS')
    ax.semilogx(frequencies, np.angle(result_full.impedance) * 180 / np.pi, 'r-', label='Full')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Impedance Phase')
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.set_ylim([-100, 100])
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Quality factor
    ax = axes[1, 0]
    ax.semilogx(frequencies, result_mqs.quality_factor, 'b-', label='MQS')
    ax.semilogx(frequencies, result_full.quality_factor, 'r-', label='Full')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Q')
    ax.set_title('Quality Factor')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Effective inductance
    ax = axes[1, 1]
    ax.semilogx(frequencies[1:], result_mqs.inductance[1:] * 1e6, 'b-', label='MQS')
    ax.semilogx(frequencies[1:], result_full.inductance[1:] * 1e6, 'r-', label='Full')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('L_eff [uH]')
    ax.set_title('Effective Inductance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_loop_star_coil.png', dpi=150)
    print("\nSaved: demo_loop_star_coil.png")

    return result_mqs, result_full


def main():
    """Main demonstration function."""
    print("=" * 70)
    print("Loop-Star Decomposition Demo for PEEC Analysis")
    print("=" * 70)

    print("\nLoop-Star decomposition separates current into:")
    print("  - Loop currents (I_L): Solenoidal, inductive (dominant at low freq)")
    print("  - Star currents (I_S): Irrotational, capacitive (significant at high freq)")
    print("\nMQS (Magneto-Quasi-Static) mode sets I_S = 0, valid for power electronics")

    # Demo 1: Simple 2-loop system
    demo_simple_loop_star()

    # Demo 2: Multi-turn coil
    demo_multi_turn_coil()

    print("\n" + "=" * 70)
    print("Demo completed successfully!")
    print("=" * 70)

    plt.show()


if __name__ == "__main__":
    main()
