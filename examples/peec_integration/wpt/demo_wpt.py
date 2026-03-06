"""
Demonstration of Multi-Port WPT (Wireless Power Transfer) Analysis

This script demonstrates:
1. Two-coil WPT system with mutual inductance
2. Coupling coefficient k vs frequency
3. Power transfer efficiency analysis
4. S-parameter computation

Physical model:
    Tx coil (L1, R1) ===k=== Rx coil (L2, R2)

    Mutual inductance: M = k * sqrt(L1 * L2)
    Z-parameters:
        Z11 = R1 + j*omega*L1  (Tx self-impedance)
        Z22 = R2 + j*omega*L2  (Rx self-impedance)
        Z12 = Z21 = j*omega*M  (Transfer impedance)

Author: Radia Development Team
Date: 2026-01-15
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from analysis import UnifiedAnalysis, MultiPortResult


def demo_basic_wpt():
    """
    Basic 2-coil WPT system analysis.

    Parameters:
        L1 = L2 = 10 uH (identical coils)
        R1 = R2 = 0.1 Ohm (copper losses)
        k = 0.3 (typical loosely coupled WPT)
    """
    print("\n" + "=" * 70)
    print("1. Basic 2-Coil WPT System")
    print("=" * 70)

    # Coil parameters
    L1 = 10e-6   # 10 uH
    L2 = 10e-6   # 10 uH
    R1 = 0.1     # 0.1 Ohm
    R2 = 0.1     # 0.1 Ohm
    k = 0.3      # Coupling coefficient

    # Mutual inductance
    M = k * np.sqrt(L1 * L2)

    print(f"Tx coil: L1 = {L1*1e6:.1f} uH, R1 = {R1:.2f} Ohm")
    print(f"Rx coil: L2 = {L2*1e6:.1f} uH, R2 = {R2:.2f} Ohm")
    print(f"Coupling coefficient k = {k:.2f}")
    print(f"Mutual inductance M = {M*1e6:.2f} uH")

    # Build inductance and resistance matrices
    L_matrix = np.array([[L1, M],
                         [M, L2]])
    R_matrix = np.array([[R1, 0],
                         [0, R2]])

    # Create analysis object
    analysis = UnifiedAnalysis()
    analysis.set_multiport_model(L_matrix, R_matrix)

    # Frequency sweep: 10 kHz to 1 MHz
    frequencies = np.logspace(4, 6, 200)

    # Run multi-port analysis
    result = analysis.frequency_sweep_multiport(frequencies)

    print(f"\nAnalysis completed in {result.computation_time:.3f} s")
    print(f"Number of frequency points: {len(frequencies)}")

    # Verify coupling coefficient from Z-parameters
    idx_mid = len(frequencies) // 2
    k_computed = result.coupling_coefficient[idx_mid]
    print(f"\nComputed k at f = {frequencies[idx_mid]/1e3:.1f} kHz: {k_computed:.4f}")
    print(f"Input k: {k:.4f}")
    print(f"Difference: {abs(k_computed - k):.6f}")

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Z11 (Tx self-impedance)
    ax = axes[0, 0]
    Z11 = result.get_Z(0, 0)
    ax.loglog(frequencies / 1e3, np.abs(Z11), label='|Z11|')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('|Z11| [Ohm]')
    ax.set_title('Tx Self-Impedance')
    ax.grid(True, which='both', alpha=0.3)

    # Z21 (Transfer impedance)
    ax = axes[0, 1]
    Z21 = result.get_Z(1, 0)
    ax.loglog(frequencies / 1e3, np.abs(Z21), label='|Z21|', color='orange')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('|Z21| [Ohm]')
    ax.set_title('Transfer Impedance (Mutual Coupling)')
    ax.grid(True, which='both', alpha=0.3)

    # Coupling coefficient vs frequency (should be constant)
    ax = axes[1, 0]
    ax.semilogx(frequencies / 1e3, result.coupling_coefficient)
    ax.axhline(y=k, color='r', linestyle='--', label=f'Input k = {k:.2f}')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Coupling coefficient k')
    ax.set_title('Coupling Coefficient vs Frequency')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Power transfer efficiency
    ax = axes[1, 1]
    ax.semilogx(frequencies / 1e3, result.power_transfer_efficiency * 100)
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Efficiency [%]')
    ax.set_title('Power Transfer Efficiency (Uncompensated)')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_wpt_basic.png', dpi=150)
    print("\nSaved: demo_wpt_basic.png")

    return result


def demo_variable_coupling():
    """
    WPT system with variable coupling coefficient.

    Simulates different coil separations/alignments.
    """
    print("\n" + "=" * 70)
    print("2. Variable Coupling Analysis")
    print("=" * 70)

    # Fixed coil parameters
    L1 = 10e-6
    L2 = 10e-6
    R1 = 0.1
    R2 = 0.1

    # Variable coupling coefficients
    k_values = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]

    # Single frequency (e.g., 85 kHz for EV charging)
    f_target = 85e3  # 85 kHz (SAE J2954 standard)
    frequencies = np.array([f_target])

    print(f"Operating frequency: {f_target/1e3:.0f} kHz")
    print(f"Coils: L1 = L2 = {L1*1e6:.0f} uH, R1 = R2 = {R1:.1f} Ohm")

    results = []

    for k in k_values:
        M = k * np.sqrt(L1 * L2)

        L_matrix = np.array([[L1, M], [M, L2]])
        R_matrix = np.array([[R1, 0], [0, R2]])

        analysis = UnifiedAnalysis()
        analysis.set_multiport_model(L_matrix, R_matrix)
        result = analysis.frequency_sweep_multiport(frequencies)

        eta = result.power_transfer_efficiency[0] * 100
        Z21_mag = np.abs(result.transfer_impedance[0])

        results.append({
            'k': k,
            'M': M,
            'eta': eta,
            'Z21': Z21_mag
        })

        print(f"k = {k:.2f}: M = {M*1e6:.2f} uH, |Z21| = {Z21_mag:.2f} Ohm, eta = {eta:.1f}%")

    # Plot k vs efficiency
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ks = [r['k'] for r in results]
    etas = [r['eta'] for r in results]
    Z21s = [r['Z21'] for r in results]

    ax = axes[0]
    ax.plot(ks, etas, 'o-', markersize=8)
    ax.set_xlabel('Coupling coefficient k')
    ax.set_ylabel('Efficiency [%]')
    ax.set_title(f'WPT Efficiency vs Coupling @ {f_target/1e3:.0f} kHz')
    ax.grid(True)

    ax = axes[1]
    ax.plot(ks, Z21s, 's-', markersize=8, color='orange')
    ax.set_xlabel('Coupling coefficient k')
    ax.set_ylabel('|Z21| [Ohm]')
    ax.set_title('Transfer Impedance vs Coupling')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('demo_wpt_variable_k.png', dpi=150)
    print("\nSaved: demo_wpt_variable_k.png")

    return results


def demo_s_parameters():
    """
    S-parameter analysis for RF WPT (e.g., 6.78 MHz ISM band).
    """
    print("\n" + "=" * 70)
    print("3. S-Parameter Analysis (RF WPT)")
    print("=" * 70)

    # RF WPT at 6.78 MHz (Airfuel standard)
    f0 = 6.78e6  # 6.78 MHz

    # Smaller inductance for higher frequency
    L1 = 1e-6   # 1 uH
    L2 = 1e-6   # 1 uH
    R1 = 0.05   # Lower resistance for high Q
    R2 = 0.05
    k = 0.2

    M = k * np.sqrt(L1 * L2)

    print(f"Operating frequency: {f0/1e6:.2f} MHz")
    print(f"L1 = L2 = {L1*1e6:.1f} uH")
    print(f"k = {k:.2f}, M = {M*1e9:.1f} nH")

    # Build matrices
    L_matrix = np.array([[L1, M], [M, L2]])
    R_matrix = np.array([[R1, 0], [0, R2]])

    analysis = UnifiedAnalysis()
    analysis.set_multiport_model(L_matrix, R_matrix)

    # Frequency sweep around resonance
    frequencies = np.linspace(1e6, 15e6, 200)
    result = analysis.frequency_sweep_multiport(frequencies)

    # Get S-parameters (50 Ohm reference)
    S = result.get_S(Z0=50.0)

    # S11 = reflection at Tx
    # S21 = transmission to Rx
    S11 = S[:, 0, 0]
    S21 = S[:, 1, 0]

    # Convert to dB
    S11_dB = 20 * np.log10(np.abs(S11) + 1e-30)
    S21_dB = 20 * np.log10(np.abs(S21) + 1e-30)

    print(f"\nAt f = {f0/1e6:.2f} MHz:")
    idx_f0 = np.argmin(np.abs(frequencies - f0))
    print(f"  |S11| = {np.abs(S11[idx_f0]):.3f} ({S11_dB[idx_f0]:.1f} dB)")
    print(f"  |S21| = {np.abs(S21[idx_f0]):.3f} ({S21_dB[idx_f0]:.1f} dB)")

    # Find minimum S11 (best match)
    idx_min_S11 = np.argmin(np.abs(S11))
    f_match = frequencies[idx_min_S11]
    print(f"\nBest match at f = {f_match/1e6:.2f} MHz:")
    print(f"  |S11| = {np.abs(S11[idx_min_S11]):.4f} ({S11_dB[idx_min_S11]:.1f} dB)")

    # Plot S-parameters
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.plot(frequencies / 1e6, S11_dB, label='|S11| (Reflection)')
    ax.axvline(x=f0 / 1e6, color='r', linestyle='--', alpha=0.5, label=f'f0 = {f0/1e6:.2f} MHz')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|S11| [dB]')
    ax.set_title('Input Reflection')
    ax.legend()
    ax.grid(True)
    ax.set_ylim([-30, 5])

    ax = axes[1]
    ax.plot(frequencies / 1e6, S21_dB, color='orange', label='|S21| (Transmission)')
    ax.axvline(x=f0 / 1e6, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|S21| [dB]')
    ax.set_title('Forward Transmission')
    ax.legend()
    ax.grid(True)
    ax.set_ylim([-30, 5])

    plt.tight_layout()
    plt.savefig('demo_wpt_s_params.png', dpi=150)
    print("\nSaved: demo_wpt_s_params.png")

    return result


def demo_multi_turn_coils():
    """
    Multi-turn coil WPT with more realistic model.

    Each coil has multiple turns with inter-turn coupling.
    """
    print("\n" + "=" * 70)
    print("4. Multi-Turn Coil WPT")
    print("=" * 70)

    # Coil parameters
    n_turns_tx = 5   # Tx turns
    n_turns_rx = 5   # Rx turns
    L_turn = 1e-6    # Inductance per turn
    R_turn = 0.02    # Resistance per turn
    k_inter = 0.9    # Inter-turn coupling (within same coil)
    k_mutual = 0.2   # Mutual coupling between coils

    print(f"Tx coil: {n_turns_tx} turns")
    print(f"Rx coil: {n_turns_rx} turns")
    print(f"L per turn: {L_turn*1e6:.1f} uH")
    print(f"Inter-turn k: {k_inter:.2f}")
    print(f"Coil-to-coil k: {k_mutual:.2f}")

    # Build full inductance matrix (all turns)
    n_total = n_turns_tx + n_turns_rx
    L_matrix = np.zeros((n_total, n_total))
    R_matrix = np.zeros((n_total, n_total))

    # Fill Tx self-inductance block
    for i in range(n_turns_tx):
        for j in range(n_turns_tx):
            if i == j:
                L_matrix[i, j] = L_turn
            else:
                L_matrix[i, j] = k_inter * L_turn

    # Fill Rx self-inductance block
    for i in range(n_turns_rx):
        for j in range(n_turns_rx):
            ii = n_turns_tx + i
            jj = n_turns_tx + j
            if i == j:
                L_matrix[ii, jj] = L_turn
            else:
                L_matrix[ii, jj] = k_inter * L_turn

    # Fill Tx-Rx mutual inductance blocks
    for i in range(n_turns_tx):
        for j in range(n_turns_rx):
            M_ij = k_mutual * L_turn
            L_matrix[i, n_turns_tx + j] = M_ij
            L_matrix[n_turns_tx + j, i] = M_ij

    # Resistance (diagonal only)
    for i in range(n_total):
        R_matrix[i, i] = R_turn

    # Effective total inductance
    L1_eff = np.sum(L_matrix[:n_turns_tx, :n_turns_tx])
    L2_eff = np.sum(L_matrix[n_turns_tx:, n_turns_tx:])
    M_eff = np.sum(L_matrix[:n_turns_tx, n_turns_tx:])
    k_eff = M_eff / np.sqrt(L1_eff * L2_eff)

    print(f"\nEffective values:")
    print(f"  L1_eff = {L1_eff*1e6:.2f} uH")
    print(f"  L2_eff = {L2_eff*1e6:.2f} uH")
    print(f"  M_eff = {M_eff*1e6:.2f} uH")
    print(f"  k_eff = {k_eff:.3f}")

    # Port grouping: all Tx turns = port 0, all Rx turns = port 1
    port_indices = [
        list(range(n_turns_tx)),           # Port 0: Tx turns
        list(range(n_turns_tx, n_total))   # Port 1: Rx turns
    ]

    # Create analysis
    analysis = UnifiedAnalysis()
    analysis.set_multiport_model(L_matrix, R_matrix, port_indices)

    # Frequency sweep
    frequencies = np.logspace(4, 6, 100)
    result = analysis.frequency_sweep_multiport(frequencies)

    print(f"\nMulti-port analysis: {result.n_ports} ports")

    # Verify effective coupling
    idx_mid = len(frequencies) // 2
    k_computed = result.coupling_coefficient[idx_mid]
    print(f"Computed k_eff at {frequencies[idx_mid]/1e3:.0f} kHz: {k_computed:.4f}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.semilogx(frequencies / 1e3, result.self_inductance[:, 0] * 1e6, label='L1_eff')
    ax.semilogx(frequencies / 1e3, result.self_inductance[:, 1] * 1e6, label='L2_eff')
    ax.axhline(y=L1_eff * 1e6, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Effective L [uH]')
    ax.set_title('Effective Self-Inductance')
    ax.legend()
    ax.grid(True)

    ax = axes[1]
    ax.semilogx(frequencies / 1e3, result.mutual_inductance * 1e6)
    ax.axhline(y=M_eff * 1e6, color='r', linestyle='--', alpha=0.5, label=f'M_eff = {M_eff*1e6:.1f} uH')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Mutual inductance [uH]')
    ax.set_title('Effective Mutual Inductance')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('demo_wpt_multi_turn.png', dpi=150)
    print("\nSaved: demo_wpt_multi_turn.png")

    return result


def demo_skin_effect():
    """
    Skin effect demonstration for conductor losses.

    Shows how AC resistance increases with frequency due to skin effect.
    Uses UnifiedSurfaceSolver with Dowell and SIBC models.
    """
    print("\n" + "=" * 70)
    print("6. Skin Effect (Conductor Losses)")
    print("=" * 70)

    from dielectric_solver import UnifiedSurfaceSolver

    # Coil parameters
    L1 = 10e-6   # 10 uH
    R_dc = 0.1   # DC resistance 0.1 Ohm

    # Conductor geometry (rectangular wire, e.g., Litz wire bundle equivalent)
    conductor_height = 2e-3   # 2 mm
    conductor_width = 3e-3    # 3 mm
    sigma = 5.8e7             # Copper conductivity

    print(f"Inductor: L = {L1*1e6:.1f} uH, R_dc = {R_dc:.3f} Ohm")
    print(f"Conductor: h = {conductor_height*1e3:.1f} mm, w = {conductor_width*1e3:.1f} mm")
    print(f"Conductivity: sigma = {sigma:.2e} S/m (copper)")

    # Frequency sweep: 10 kHz to 1 MHz
    frequencies = np.logspace(4, 6, 100)

    # Test 1: Dowell model
    print("\n--- Dowell Model ---")
    solver_dowell = UnifiedSurfaceSolver()
    solver_dowell.set_conductor_loop(
        L=np.array([[L1]]),
        R=np.array([[R_dc]])
    )
    solver_dowell.set_skin_effect(
        sigma=sigma,
        model='dowell',
        conductor_height=conductor_height,
        conductor_width=conductor_width
    )

    result_dowell = solver_dowell.solve_frequency(frequencies)
    R_dowell = result_dowell['resistance']
    L_dowell = result_dowell['inductance']

    # Test 2: SIBC model
    print("--- SIBC Model ---")
    solver_sibc = UnifiedSurfaceSolver()
    solver_sibc.set_conductor_loop(
        L=np.array([[L1]]),
        R=np.array([[R_dc]])
    )
    solver_sibc.set_skin_effect(
        sigma=sigma,
        model='sibc',
        conductor_height=conductor_height
    )

    result_sibc = solver_sibc.solve_frequency(frequencies)
    R_sibc = result_sibc['resistance']

    # Test 3: No skin effect (reference)
    solver_noskin = UnifiedSurfaceSolver()
    solver_noskin.set_conductor_loop(
        L=np.array([[L1]]),
        R=np.array([[R_dc]])
    )
    result_noskin = solver_noskin.solve_frequency(frequencies)
    R_noskin = result_noskin['resistance']
    L_noskin = result_noskin['inductance']

    # Compute skin depth at key frequencies
    mu0 = 4 * np.pi * 1e-7
    for f in [10e3, 100e3, 1e6]:
        delta = np.sqrt(2 / (2 * np.pi * f * mu0 * sigma))
        print(f"  Skin depth at {f/1e3:.0f} kHz: delta = {delta*1e3:.3f} mm")

    # Report R_ac/R_dc ratio
    print("\nR_ac / R_dc ratio:")
    for i, f in enumerate([10e3, 100e3, 1e6]):
        idx = np.argmin(np.abs(frequencies - f))
        print(f"  {f/1e3:.0f} kHz: Dowell = {R_dowell[idx]/R_dc:.2f}, SIBC = {np.real(R_sibc[idx])/R_dc:.2f}")

    # Plot results
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # AC resistance
    ax = axes[0]
    ax.loglog(frequencies / 1e3, R_noskin, 'b--', label='No skin effect', alpha=0.5)
    ax.loglog(frequencies / 1e3, R_dowell, 'r-', label='Dowell model', linewidth=2)
    ax.loglog(frequencies / 1e3, np.real(R_sibc), 'g-', label='SIBC model', linewidth=2)
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('AC Resistance [Ohm]')
    ax.set_title('Skin Effect: AC Resistance vs Frequency')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Inductance (Dowell only - SIBC doesn't change L)
    ax = axes[1]
    ax.semilogx(frequencies / 1e3, L_noskin * 1e6, 'b--', label='No skin effect', alpha=0.5)
    ax.semilogx(frequencies / 1e3, L_dowell * 1e6, 'r-', label='Dowell model', linewidth=2)
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Inductance [uH]')
    ax.set_title('Internal Inductance Reduction (Dowell)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_wpt_skin_effect.png', dpi=150)
    print("\nSaved: demo_wpt_skin_effect.png")

    return {
        'frequencies': frequencies,
        'R_dowell': R_dowell,
        'R_sibc': R_sibc,
        'L_dowell': L_dowell
    }


def demo_compensation_design():
    """
    Compensation capacitor design for resonant WPT.
    """
    print("\n" + "=" * 70)
    print("5. Resonant Compensation Design")
    print("=" * 70)

    # Coil parameters
    L1 = 10e-6   # 10 uH
    L2 = 10e-6   # 10 uH
    R1 = 0.1
    R2 = 0.1
    k = 0.3

    M = k * np.sqrt(L1 * L2)

    # Target frequency
    f0 = 85e3  # 85 kHz (EV charging)

    # Design compensation capacitors
    analysis = UnifiedAnalysis()
    C1 = analysis.design_compensation_capacitor(L1, f0)
    C2 = analysis.design_compensation_capacitor(L2, f0)

    print(f"Target resonant frequency: {f0/1e3:.0f} kHz")
    print(f"L1 = {L1*1e6:.1f} uH -> C1 = {C1*1e9:.2f} nF")
    print(f"L2 = {L2*1e6:.1f} uH -> C2 = {C2*1e9:.2f} nF")

    # Verify resonance
    f1_res = analysis.find_resonant_frequency(L1, C1)
    f2_res = analysis.find_resonant_frequency(L2, C2)
    print(f"\nVerification:")
    print(f"  Tx resonance: {f1_res/1e3:.2f} kHz")
    print(f"  Rx resonance: {f2_res/1e3:.2f} kHz")

    # Q factors
    omega0 = 2 * np.pi * f0
    Q1 = omega0 * L1 / R1
    Q2 = omega0 * L2 / R2
    print(f"\nQ factors at resonance:")
    print(f"  Q1 = {Q1:.0f}")
    print(f"  Q2 = {Q2:.0f}")

    # Maximum efficiency at resonance
    # eta_max = k^2 * Q1 * Q2 / (1 + sqrt(1 + k^2 * Q1 * Q2))^2
    kQ_sq = k ** 2 * Q1 * Q2
    eta_max = kQ_sq / (1 + np.sqrt(1 + kQ_sq)) ** 2
    print(f"\nMaximum efficiency (at optimal load):")
    print(f"  eta_max = {eta_max*100:.1f}%")

    return {
        'f0': f0,
        'C1': C1,
        'C2': C2,
        'Q1': Q1,
        'Q2': Q2,
        'eta_max': eta_max
    }


def main():
    """Main demonstration."""
    print("=" * 70)
    print("Multi-Port WPT Analysis Demonstration")
    print("=" * 70)

    # Test 1: Basic WPT
    demo_basic_wpt()

    # Test 2: Variable coupling
    demo_variable_coupling()

    # Test 3: S-parameters
    demo_s_parameters()

    # Test 4: Multi-turn coils
    demo_multi_turn_coils()

    # Test 5: Compensation design
    demo_compensation_design()

    # Test 6: Skin effect (conductor losses)
    demo_skin_effect()

    print("\n" + "=" * 70)
    print("Demo completed!")
    print("=" * 70)

    plt.show()


if __name__ == "__main__":
    main()
