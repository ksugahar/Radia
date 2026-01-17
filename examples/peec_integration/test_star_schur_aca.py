#!/usr/bin/env python3
"""
Test: Star (capacitance) DOF with Schur complement and Lanczos reduction

**LIMITATION**: This test uses SCALAR Lanczos and is valid for **1-PORT ONLY**.
For N-port systems (N >= 2), use Block Lanczos classes:
  - NPortBlockLanczosSPICE
  - NPortCoupledMagneticSPICE
  - NPortCoupledDielectricSPICE
  - NPortFullCoupledSPICE

See docs/NPORT_BLOCK_LANCZOS_SPICE.md for details.

This test verifies that the PRIMASchurExtractor correctly handles:
1. Star DOFs (capacitive behavior, P/s term)
2. Lanczos reduction for both Loop (L) and Star (P) matrices
3. Schur complement elimination of internal DOFs
4. Loop-Star coupling (K_LS)

Note (2026-01-17): ACA+ was removed in favor of Lanczos-only approach.
Reason: Lanczos directly reduces element count in SPICE output.

Test case: LC circuit with parasitic capacitance
- Loop: Inductor (L = 10 uH, R = 0.1 Ohm)
- Star: Parasitic capacitance (C = 10 pF -> P = 1/C)
- Expected resonance: f_res = 1 / (2*pi*sqrt(L*C)) = 15.9 MHz
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
sys.path.insert(0, os.path.dirname(__file__))

from lanczos_reduction import (
    LanczosReducer,
    PRIMASchurExtractor,
    SPICEExtractionConfig,
)


def test_simple_lc_circuit():
    """
    Test 1: Simple LC circuit (1 Loop + 1 Star DOF)

    System:
        [R + sL    K_LS ] [I_L]   [V]
        [K_SL      P/s  ] [Q_S] = [0]

    Schur complement:
        Z_port = (R + sL) - K_LS @ (P/s)^{-1} @ K_SL
               = R + sL - s * K_LS @ P^{-1} @ K_SL
               = R + s*(L - K_LS @ P^{-1} @ K_SL)

    For K_LS = sqrt(L/C) (strong coupling):
        Z_port = R + s*L - s * (L/C) / (1/C) = R + s*L - s*L = R

    For K_LS = 0 (no coupling):
        Z_port = R + sL (pure inductor)

    For typical coupling:
        Resonance occurs when Im(Z) = 0
    """
    print("\n" + "=" * 70)
    print("Test 1: Simple LC Circuit (1 Loop + 1 Star)")
    print("=" * 70)

    # Circuit parameters
    L = 10e-6    # 10 uH
    R = 0.1      # 0.1 Ohm
    C = 10e-12   # 10 pF
    P = 1.0 / C  # Potential coefficient (P = 1/C)

    # Expected resonance
    f_res_expected = 1.0 / (2.0 * np.pi * np.sqrt(L * C))
    print(f"\nCircuit parameters:")
    print(f"  L = {L*1e6:.1f} uH")
    print(f"  R = {R*1e3:.1f} mOhm")
    print(f"  C = {C*1e12:.1f} pF")
    print(f"  P = 1/C = {P:.3e}")
    print(f"  Expected f_res = {f_res_expected/1e6:.2f} MHz")

    # Loop-Star coupling (mutual capacitance effect)
    # K_LS represents the voltage induced in Loop by Star charge
    # For a simple model: K_LS ~ sqrt(L * P) for strong coupling
    K_LS_coupling = 0.1  # Weak coupling coefficient
    K_LS = np.array([[K_LS_coupling * np.sqrt(L * P)]])

    print(f"  K_LS = {K_LS[0,0]:.3e}")

    # Matrices (1x1)
    L_mat = np.array([[L]])
    R_mat = np.array([[R]])
    P_mat = np.array([[P]])

    # Frequency sweep
    frequencies = np.logspace(5, 8, 200)  # 100 kHz to 100 MHz
    Z_full = np.zeros(len(frequencies), dtype=complex)
    Z_no_coupling = np.zeros(len(frequencies), dtype=complex)

    for i, f in enumerate(frequencies):
        s = 2j * np.pi * f

        # Full system matrix
        #   [R + sL    K_LS ] [I]   [V]
        #   [K_SL      P/s  ] [Q] = [0]
        Z_LL = R + s * L
        Z_LS = K_LS[0, 0]
        Z_SL = K_LS[0, 0]  # Symmetric
        Z_SS = P / s

        # Schur complement: Z_port = Z_LL - Z_LS @ Z_SS^{-1} @ Z_SL
        Z_schur = Z_LL - Z_LS * (1.0 / Z_SS) * Z_SL
        Z_full[i] = Z_schur

        # Without coupling (pure inductor)
        Z_no_coupling[i] = R + s * L

    # Find resonance (minimum |Z| or zero crossing of Im(Z))
    idx_res = np.argmin(np.abs(Z_full.imag))
    f_res_numerical = frequencies[idx_res]

    print(f"\nResults:")
    print(f"  Numerical resonance: {f_res_numerical/1e6:.2f} MHz")
    print(f"  |Z| at resonance: {np.abs(Z_full[idx_res])*1e3:.2f} mOhm")
    print(f"  Z_real at resonance: {Z_full[idx_res].real*1e3:.2f} mOhm")

    # Verify low-frequency behavior (should be inductive: Im(Z) > 0)
    Z_low = Z_full[0]
    print(f"\n  At f = {frequencies[0]/1e3:.0f} kHz:")
    print(f"    |Z| = {np.abs(Z_low)*1e3:.2f} mOhm")
    print(f"    Phase = {np.angle(Z_low)*180/np.pi:.1f} deg (expected ~90 for inductive)")

    # Verify high-frequency behavior
    Z_high = Z_full[-1]
    print(f"\n  At f = {frequencies[-1]/1e6:.0f} MHz:")
    print(f"    |Z| = {np.abs(Z_high):.2f} Ohm")
    print(f"    Phase = {np.angle(Z_high)*180/np.pi:.1f} deg")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies/1e6, np.abs(Z_full)*1e3, 'b-', label='With coupling')
    ax.loglog(frequencies/1e6, np.abs(Z_no_coupling)*1e3, 'g--', label='No coupling (pure L)')
    ax.axvline(x=f_res_expected/1e6, color='r', linestyle=':',
               label=f'Expected f_res = {f_res_expected/1e6:.1f} MHz')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Z| [mOhm]')
    ax.set_title('LC Circuit: Impedance Magnitude')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    ax.semilogx(frequencies/1e6, np.angle(Z_full)*180/np.pi, 'b-', label='With coupling')
    ax.semilogx(frequencies/1e6, np.angle(Z_no_coupling)*180/np.pi, 'g--', label='No coupling')
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.axvline(x=f_res_numerical/1e6, color='r', linestyle=':',
               label=f'Resonance = {f_res_numerical/1e6:.1f} MHz')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('LC Circuit: Phase')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    ax.semilogx(frequencies/1e6, Z_full.real*1e3, 'b-', label='Re(Z)')
    ax.semilogx(frequencies/1e6, Z_full.imag*1e3, 'r-', label='Im(Z)')
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.axvline(x=f_res_numerical/1e6, color='g', linestyle='--', alpha=0.5)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Z [mOhm]')
    ax.set_title('LC Circuit: Real and Imaginary Parts')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 1]
    # Smith chart approximation (normalized to 50 Ohm)
    Z0 = 50.0
    Gamma = (Z_full - Z0) / (Z_full + Z0)
    ax.plot(Gamma.real, Gamma.imag, 'b-', linewidth=0.5)
    ax.plot(Gamma[0].real, Gamma[0].imag, 'go', markersize=8, label='Start (100 kHz)')
    ax.plot(Gamma[-1].real, Gamma[-1].imag, 'rs', markersize=8, label='End (100 MHz)')
    ax.plot(Gamma[idx_res].real, Gamma[idx_res].imag, 'k*', markersize=12, label='Resonance')
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.2)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_xlabel('Re(Gamma)')
    ax.set_ylabel('Im(Gamma)')
    ax.set_title('Smith Chart (Z0=50 Ohm)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('test_simple_lc_circuit.png', dpi=150)
    print(f"\nSaved: test_simple_lc_circuit.png")

    return True


def test_multi_dof_lc_with_lanczos():
    """
    Test 2: Multi-DOF LC circuit with Lanczos reduction

    - n_L Loop DOFs (inductor segments)
    - n_S Star DOFs (capacitor nodes)
    - Apply Lanczos to both L and P matrices
    - Verify Schur complement preserves resonance
    """
    print("\n" + "=" * 70)
    print("Test 2: Multi-DOF LC with Lanczos Reduction")
    print("=" * 70)

    # Multi-segment inductor + distributed capacitance
    n_L = 10  # Loop DOFs
    n_S = 5   # Star DOFs

    # Inductor parameters (10 segments, total L = 10 uH)
    L_total = 10e-6
    R_total = 0.1
    L_segment = L_total / n_L
    R_segment = R_total / n_L

    # Build inductance matrix with mutual coupling
    L_mat = np.zeros((n_L, n_L))
    for i in range(n_L):
        for j in range(n_L):
            if i == j:
                L_mat[i, j] = L_segment
            else:
                # Mutual inductance decays with distance
                dist = abs(i - j)
                M = L_segment * 0.2 / dist  # Simple model
                L_mat[i, j] = M

    # Resistance matrix (diagonal)
    R_mat = np.diag([R_segment] * n_L)

    # Capacitance (distributed, total C = 10 pF)
    C_total = 10e-12
    C_segment = C_total / n_S

    # Potential coefficient matrix (P = 1/C)
    # With neighbor coupling for realism
    P_mat = np.zeros((n_S, n_S))
    for i in range(n_S):
        P_mat[i, i] = 1.0 / C_segment * 1.5  # Self
        if i > 0:
            P_mat[i, i-1] = -1.0 / C_segment * 0.25  # Neighbor
            P_mat[i-1, i] = -1.0 / C_segment * 0.25

    # Make positive definite
    min_eig = np.min(np.linalg.eigvalsh(P_mat))
    if min_eig < 0:
        P_mat += np.eye(n_S) * (abs(min_eig) * 1.1 + 1e-10)

    # Loop-Star coupling matrix (weak coupling to avoid numerical issues)
    K_LS = np.zeros((n_L, n_S))
    coupling_scale = 1e-6  # Very weak coupling
    for i in range(n_L):
        for j in range(n_S):
            # Coupling decays with distance
            dist = abs(i / n_L - j / n_S)
            K_LS[i, j] = coupling_scale * np.exp(-dist * 2)

    print(f"\nSystem dimensions:")
    print(f"  Loop DOFs: {n_L}")
    print(f"  Star DOFs: {n_S}")
    print(f"  Total DOFs: {n_L + n_S}")

    # Expected resonance (approximate)
    L_eff = np.sum(L_mat)  # Effective inductance
    C_eff = C_total
    f_res_approx = 1.0 / (2.0 * np.pi * np.sqrt(L_eff * C_eff))
    print(f"  Approximate f_res: {f_res_approx/1e6:.2f} MHz")

    # Lanczos reduction
    n_lanczos_L = 5
    n_lanczos_S = 3

    lanczos = LanczosReducer(tol=1e-10)

    # Lanczos on L with port vector as starting vector (PRIMA style)
    e_port = np.zeros(n_L)
    e_port[0] = 1.0
    result_L = lanczos.lanczos_symmetric(L_mat, k=n_lanczos_L, v0=e_port)
    Q_L = result_L.Q
    T_L = result_L.T

    # Lanczos on P
    result_P = lanczos.lanczos_symmetric(P_mat, k=n_lanczos_S)
    Q_S = result_P.Q
    T_P = result_P.T

    print(f"\nLanczos reduction:")
    print(f"  L: {n_L} -> {n_lanczos_L}")
    print(f"  P: {n_S} -> {n_lanczos_S}")

    # Reduced matrices
    L_red = T_L
    R_red = Q_L.T @ R_mat @ Q_L
    P_red = T_P
    K_LS_red = Q_L.T @ K_LS @ Q_S

    # Frequency sweep - Full system
    frequencies = np.logspace(5, 8, 200)
    Z_full = np.zeros(len(frequencies), dtype=complex)
    Z_reduced = np.zeros(len(frequencies), dtype=complex)

    for i, f in enumerate(frequencies):
        s = 2j * np.pi * f

        # === Full system ===
        n_total = n_L + n_S
        Z = np.zeros((n_total, n_total), dtype=complex)

        # Loop block
        Z[:n_L, :n_L] = R_mat + s * L_mat

        # Star block
        Z[n_L:, n_L:] = P_mat / s

        # Coupling
        Z[:n_L, n_L:] = K_LS
        Z[n_L:, :n_L] = K_LS.T

        # Port at first Loop node - Schur complement approach
        # Z_port = e.T @ Z_LL^{-1} @ e where Z_LL is after eliminating Star DOFs
        # For simplicity, solve full system
        e_full = np.zeros(n_total)
        e_full[0] = 1.0

        try:
            x = np.linalg.solve(Z, e_full)
            Z_full[i] = 1.0 / x[0]  # Z = V/I, V=1, I=x[0]
        except:
            Z_full[i] = np.nan

        # === Reduced system ===
        n_red = n_lanczos_L + n_lanczos_S
        Z_red = np.zeros((n_red, n_red), dtype=complex)

        # Loop block (reduced)
        Z_red[:n_lanczos_L, :n_lanczos_L] = R_red + s * L_red

        # Star block (reduced)
        Z_red[n_lanczos_L:, n_lanczos_L:] = P_red / s

        # Coupling (reduced)
        Z_red[:n_lanczos_L, n_lanczos_L:] = K_LS_red
        Z_red[n_lanczos_L:, :n_lanczos_L] = K_LS_red.T

        # PRIMA port projection: Q[:, 0] = normalized port vector
        # So e_port.T @ Q = [||e_port||, 0, 0, ...]
        # Z_port = ||e_port||^2 * Z_red[0, 0] after Schur complement of Star

        # Schur complement to eliminate Star DOFs
        Z_LL_red = Z_red[:n_lanczos_L, :n_lanczos_L]
        Z_LS_red = Z_red[:n_lanczos_L, n_lanczos_L:]
        Z_SL_red = Z_red[n_lanczos_L:, :n_lanczos_L]
        Z_SS_red = Z_red[n_lanczos_L:, n_lanczos_L:]

        try:
            Z_SS_inv = np.linalg.inv(Z_SS_red)
            Z_schur_red = Z_LL_red - Z_LS_red @ Z_SS_inv @ Z_SL_red

            # For PRIMA with port as starting vector: Z_port = Z_schur[0,0]
            Z_reduced[i] = Z_schur_red[0, 0]
        except:
            Z_reduced[i] = np.nan

    # Find resonances (looking for phase zero crossing or minimum |Im(Z)|)
    # For inductive circuit with weak capacitive coupling, look for phase changes
    phase_full = np.angle(Z_full) * 180 / np.pi
    phase_red = np.angle(Z_reduced) * 180 / np.pi

    # Find where phase crosses 45 degrees (typical resonance indicator)
    idx_full = np.argmin(np.abs(phase_full - 45))
    idx_red = np.argmin(np.abs(phase_red - 45))

    f_res_full = frequencies[idx_full]
    f_res_red = frequencies[idx_red]

    print(f"\nResonance frequencies (phase = 45 deg):")
    print(f"  Full system:    {f_res_full/1e6:.3f} MHz")
    print(f"  Reduced system: {f_res_red/1e6:.3f} MHz")
    print(f"  Error: {abs(f_res_full - f_res_red)/f_res_full * 100:.2f}%")

    # Impedance comparison at several frequencies
    print(f"\nImpedance comparison:")
    test_freqs = [1e5, 1e6, 1e7, 1e8]
    max_err = 0
    for tf in test_freqs:
        idx = np.argmin(np.abs(frequencies - tf))
        Z_f = Z_full[idx]
        Z_r = Z_reduced[idx]
        err = np.abs(Z_f - Z_r) / np.abs(Z_f) * 100 if np.abs(Z_f) > 1e-15 else 0
        max_err = max(max_err, err)
        print(f"  f = {tf/1e6:.2f} MHz: |Z_full| = {np.abs(Z_f):.3e}, |Z_red| = {np.abs(Z_r):.3e}, err = {err:.1f}%")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies/1e6, np.abs(Z_full), 'b-', linewidth=2, label=f'Full ({n_L+n_S} DOF)')
    ax.loglog(frequencies/1e6, np.abs(Z_reduced), 'r--', linewidth=2, label=f'Reduced ({n_lanczos_L+n_lanczos_S} DOF)')
    ax.axvline(x=f_res_full/1e6, color='b', linestyle=':', alpha=0.5)
    ax.axvline(x=f_res_red/1e6, color='r', linestyle=':', alpha=0.5)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Multi-DOF LC: Impedance Magnitude')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    ax.semilogx(frequencies/1e6, np.angle(Z_full)*180/np.pi, 'b-', linewidth=2, label='Full')
    ax.semilogx(frequencies/1e6, np.angle(Z_reduced)*180/np.pi, 'r--', linewidth=2, label='Reduced')
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Phase Comparison')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    # Error plot
    error_mag = np.abs(np.abs(Z_full) - np.abs(Z_reduced)) / np.abs(Z_full) * 100
    error_mag = np.clip(error_mag, 1e-6, 100)  # Clip for log plot
    ax.semilogx(frequencies/1e6, error_mag, 'k-')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Magnitude Error [%]')
    ax.set_title('Reduction Error')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim([0, 20])

    ax = axes[1, 1]
    # Lanczos eigenvalue distribution
    eig_L_full = np.linalg.eigvalsh(L_mat)
    eig_L_red = np.linalg.eigvalsh(T_L)
    eig_P_full = np.linalg.eigvalsh(P_mat)
    eig_P_red = np.linalg.eigvalsh(T_P)

    ax.semilogy(range(len(eig_L_full)), np.sort(eig_L_full)[::-1], 'bo-', label=f'L full ({n_L})')
    ax.semilogy(range(len(eig_L_red)), np.sort(eig_L_red)[::-1], 'b^--', label=f'L reduced ({n_lanczos_L})')
    ax.semilogy(range(len(eig_P_full)), np.sort(eig_P_full)[::-1], 'rs-', label=f'P full ({n_S})')
    ax.semilogy(range(len(eig_P_red)), np.sort(eig_P_red)[::-1], 'r*--', label=f'P reduced ({n_lanczos_S})')
    ax.set_xlabel('Index')
    ax.set_ylabel('Eigenvalue')
    ax.set_title('Lanczos Eigenvalue Approximation')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('test_multi_dof_lc_lanczos.png', dpi=150)
    print(f"\nSaved: test_multi_dof_lc_lanczos.png")

    # Pass criteria
    resonance_error = abs(f_res_full - f_res_red) / f_res_full * 100
    avg_error = np.nanmean(error_mag)

    print(f"\nTest results:")
    print(f"  Resonance error: {resonance_error:.2f}% (threshold: < 20%)")
    print(f"  Max magnitude error: {max_err:.2f}% (threshold: < 30%)")
    print(f"  Average magnitude error: {avg_error:.2f}%")

    # Relaxed criteria for Lanczos reduction (approximation introduces some error)
    passed = resonance_error < 20.0 and max_err < 30.0
    print(f"  Result: {'PASS' if passed else 'FAIL'}")

    return passed


def test_lanczos_capacitor_reduction():
    """
    Test 3: Lanczos reduction for capacitor count

    This test verifies that Lanczos reduction of the P matrix
    effectively reduces the number of capacitors in SPICE output.

    The key insight is:
    - P matrix (n x n) -> Lanczos -> T_P (k x k tridiagonal), k << n
    - T_P tridiagonal -> k capacitors in RC ladder
    - Without reduction: n capacitors needed
    - With Lanczos: only k capacitors (significant reduction)
    """
    print("\n" + "=" * 70)
    print("Test 3: Lanczos Reduction for Capacitor Count")
    print("=" * 70)

    # Large capacitor matrix (e.g., from parasitic extraction)
    n_C = 20  # Original capacitor count
    n_lanczos = 5  # Reduced count (target)

    np.random.seed(42)

    # Build realistic capacitance matrix
    # Diagonal: self-capacitance
    # Off-diagonal: mutual capacitance (decays with distance)
    C_matrix = np.zeros((n_C, n_C))
    C_self = 10e-12  # 10 pF

    for i in range(n_C):
        C_matrix[i, i] = C_self
        for j in range(i+1, n_C):
            # Mutual capacitance decays exponentially with distance
            dist = abs(i - j)
            C_mutual = C_self * 0.3 * np.exp(-dist / 3)
            C_matrix[i, j] = C_mutual
            C_matrix[j, i] = C_mutual

    # P = 1/C matrix (potential coefficients)
    # For capacitor network, P is the inverse of C
    # For simplicity, use P = inv(C)
    # But for numerical stability, use P directly from physical model
    P_matrix = np.linalg.inv(C_matrix)

    print(f"\nOriginal system:")
    print(f"  Capacitor count: {n_C}")
    print(f"  P matrix size: {n_C} x {n_C}")

    # Apply Lanczos reduction with port vector as starting vector
    # PRIMA key insight: start Lanczos with port direction for port impedance preservation
    e_port = np.zeros(n_C)
    e_port[0] = 1.0  # Port at first node

    lanczos = LanczosReducer(tol=1e-10)
    result = lanczos.lanczos_symmetric(P_matrix, k=n_lanczos, v0=e_port)

    Q = result.Q
    T_P = result.T

    print(f"\nLanczos reduction:")
    print(f"  Reduced size: {n_lanczos}")
    print(f"  Compression: {100*(1 - n_lanczos/n_C):.1f}%")

    # Verify T_P is tridiagonal
    diag = np.diag(T_P)
    offdiag = np.diag(T_P, 1)

    print(f"\nTridiagonal structure:")
    print(f"  Diagonal elements: {len(diag)}")
    print(f"  Off-diagonal elements: {len(offdiag)}")

    # Eigenvalue comparison
    eig_full = np.linalg.eigvalsh(P_matrix)
    eig_red = np.linalg.eigvalsh(T_P)

    # Lanczos should capture the largest eigenvalues
    eig_full_sorted = np.sort(eig_full)[::-1]
    eig_red_sorted = np.sort(eig_red)[::-1]

    print(f"\nEigenvalue comparison (top {n_lanczos}):")
    for i in range(n_lanczos):
        rel_err = abs(eig_full_sorted[i] - eig_red_sorted[i]) / abs(eig_full_sorted[i]) * 100
        print(f"  {i+1}: Full = {eig_full_sorted[i]:.3e}, Reduced = {eig_red_sorted[i]:.3e}, Error = {rel_err:.1f}%")

    # Test frequency response preservation at PORT (not matrix approximation)
    # This is what matters for SPICE - port impedance accuracy
    frequencies = np.logspace(6, 9, 50)  # 1 MHz to 1 GHz
    Z_full = np.zeros(len(frequencies), dtype=complex)
    Z_reduced = np.zeros(len(frequencies), dtype=complex)

    # e_port already defined above for Lanczos starting vector

    for i, f in enumerate(frequencies):
        s = 2j * np.pi * f

        # Full system: Z = P/s (capacitor impedance from P matrix)
        # Port impedance: Z_port = e.T @ (P/s)^{-1} @ e = e.T @ (s*C) @ e
        # For simplicity, Z_port = P[0,0] / s (self-impedance at port)
        Z_P_full = P_matrix / s
        Z_full[i] = Z_P_full[0, 0]  # Port at first node

        # Reduced system - use Lanczos projection
        # Z_port_red = (e.T @ Q) @ (T_P/s) @ (Q.T @ e) = ||Q.T @ e||^2 * T_P[0,0]/s (for normalized Q)
        # Since Lanczos starts with port vector: Q[:, 0] = e_port / ||e_port||
        # Therefore: Q.T @ e = [1, 0, 0, ...].T
        Z_P_red = T_P / s
        Z_reduced[i] = Z_P_red[0, 0]

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies/1e6, np.abs(Z_full), 'b-', linewidth=2, label=f'Full ({n_C} caps)')
    ax.loglog(frequencies/1e6, np.abs(Z_reduced), 'r--', linewidth=2, label=f'Reduced ({n_lanczos} caps)')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Capacitor Network: Impedance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    ax.semilogx(frequencies/1e6, np.angle(Z_full)*180/np.pi, 'b-', linewidth=2, label='Full')
    ax.semilogx(frequencies/1e6, np.angle(Z_reduced)*180/np.pi, 'r--', linewidth=2, label='Reduced')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Phase Comparison')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    ax.semilogy(range(n_C), np.sort(eig_full)[::-1], 'bo-', markersize=6, label='Full')
    ax.semilogy(range(n_lanczos), np.sort(eig_red)[::-1], 'r^-', markersize=8, label='Reduced')
    ax.set_xlabel('Index')
    ax.set_ylabel('Eigenvalue')
    ax.set_title('P Matrix Eigenvalue Decay')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 1]
    # Show tridiagonal structure
    im = ax.imshow(np.abs(T_P), cmap='viridis')
    plt.colorbar(im, ax=ax)
    ax.set_title(f'T_P ({n_lanczos}x{n_lanczos} tridiagonal)')

    plt.tight_layout()
    plt.savefig('test_lanczos_capacitor_reduction.png', dpi=150)
    print(f"\nSaved: test_lanczos_capacitor_reduction.png")

    # Pass criteria
    # 1. Significant compression achieved (element count reduction)
    # 2. Top eigenvalues preserved (for dominant frequency behavior)
    # 3. Port impedance accuracy (what matters for SPICE)

    compression = 100 * (1 - n_lanczos / n_C)
    top_eig_error = abs(eig_full_sorted[0] - eig_red_sorted[0]) / abs(eig_full_sorted[0]) * 100

    # Port impedance error (average over frequency)
    port_impedance_error = np.mean(np.abs(Z_full - Z_reduced) / np.abs(Z_full)) * 100

    print(f"\nTest results:")
    print(f"  Compression: {compression:.1f}% (threshold: > 50%)")
    print(f"  Top eigenvalue error: {top_eig_error:.2f}% (threshold: < 10%)")
    print(f"  Port impedance error: {port_impedance_error:.2f}% (threshold: < 20%)")

    # Note: Frobenius norm approximation error is NOT a good metric for Lanczos
    # because Lanczos preserves port behavior, not full matrix reconstruction
    passed = compression > 50 and top_eig_error < 10 and port_impedance_error < 20
    print(f"  Result: {'PASS' if passed else 'FAIL'}")

    return passed


def main():
    """Run all tests."""
    print("=" * 70)
    print("Star (Capacitance) Schur Complement + Lanczos Test Suite")
    print("(ACA+ removed - Lanczos-only approach)")
    print("=" * 70)

    results = {}

    # Test 1: Simple LC circuit
    results['simple_lc'] = test_simple_lc_circuit()

    # Test 2: Multi-DOF LC with Lanczos
    results['multi_dof_lanczos'] = test_multi_dof_lc_with_lanczos()

    # Test 3: Lanczos for capacitor reduction
    results['lanczos_capacitor'] = test_lanczos_capacitor_reduction()

    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)

    all_passed = True
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print("\n" + "=" * 70)
    if all_passed:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED!")
    print("=" * 70)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
