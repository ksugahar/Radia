"""
Verify PRIMA I frequency response matches original PEEC Loop impedance

This script verifies that the PRIMA I transformation preserves the
frequency-dependent impedance of the original PEEC circuit.

Test 1: Without surface impedance Z_s
    Original: Z(s) = R + s*L  (dense matrices)
    PRIMA I:    Z'(s) = R_diag + s*L_tridiag  (sparse matrices)

Test 2: With surface impedance Z_s (frequency-dependent)
    Original: Z(s) = R + s*L + Z_s(s)*I  (dense matrices)
    PRIMA I:    Z'(s) = R_diag + s*L_tridiag + Z_s(s)*I  (sparse matrices)

The impedance should be identical when measured from the same port.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Add PRIMA library path
sys.path.insert(0, r'W:\30_CauerLadderNetwork\2021_01_22_CauerI_to_CauerII\Python')

from prima import lanczos


def compute_surface_impedance(frequencies, sigma, mu_r=1.0):
    """
    Compute surface impedance Z_s(f) for conductor with skin effect.

    Z_s = (1 + j) / (sigma * delta) = (1 + j) * sqrt(pi * f * mu / sigma)

    Args:
        frequencies: Array of frequencies [Hz]
        sigma: Conductivity [S/m]
        mu_r: Relative permeability (default 1.0)

    Returns:
        Z_s: Complex surface impedance array (len(frequencies),)
    """
    mu_0 = 4 * np.pi * 1e-7  # H/m
    mu = mu_r * mu_0

    Z_s = np.zeros(len(frequencies), dtype=complex)
    for i, f in enumerate(frequencies):
        if f > 0:
            # Z_s = (1+j) * sqrt(pi * f * mu / sigma)
            Z_s[i] = (1 + 1j) * np.sqrt(np.pi * f * mu / sigma)
        else:
            Z_s[i] = 0.0
    return Z_s


def compute_impedance(R, L, frequencies, port_vector, Z_s=None):
    """
    Compute impedance Z(f) = v^T * (R + j*omega*L + Z_s*I)^{-1} * v

    Args:
        R: Resistance matrix (n x n)
        L: Inductance matrix (n x n)
        frequencies: Array of frequencies [Hz]
        port_vector: Port excitation vector (n,)
        Z_s: Surface impedance array (len(frequencies),) or None

    Returns:
        Z: Complex impedance array (len(frequencies),)
    """
    n = R.shape[0]
    Z = np.zeros(len(frequencies), dtype=complex)
    I_mat = np.eye(n)

    for i, f in enumerate(frequencies):
        omega = 2 * np.pi * f
        # Z_matrix = R + j*omega*L + Z_s*I
        Z_matrix = R + 1j * omega * L
        if Z_s is not None:
            Z_matrix = Z_matrix + Z_s[i] * I_mat
        # Solve Z_matrix * I = v for I, then Z = v^T * I
        try:
            I = np.linalg.solve(Z_matrix, port_vector)
            Z[i] = np.dot(port_vector, I)
        except np.linalg.LinAlgError:
            Z[i] = np.nan

    return Z


def compute_impedance_with_transformed_Zs(R, L, frequencies, port_vector, Z_s, I_transformed):
    """
    Compute impedance Z(f) = v^T * (R + j*omega*L + Z_s*I')^{-1} * v

    where I' = U^H * I * V is the transformed identity matrix.

    Args:
        R: Resistance matrix (n x n)
        L: Inductance matrix (n x n)
        frequencies: Array of frequencies [Hz]
        port_vector: Port excitation vector (n,)
        Z_s: Surface impedance array (len(frequencies),)
        I_transformed: Transformed identity matrix U^H*I*V (n x n)

    Returns:
        Z: Complex impedance array (len(frequencies),)
    """
    n = R.shape[0]
    Z = np.zeros(len(frequencies), dtype=complex)

    for i, f in enumerate(frequencies):
        omega = 2 * np.pi * f
        # Z_matrix = R + j*omega*L + Z_s*I'
        Z_matrix = R + 1j * omega * L + Z_s[i] * I_transformed
        # Solve Z_matrix * I = v for I, then Z = v^T * I
        try:
            I = np.linalg.solve(Z_matrix, port_vector)
            Z[i] = np.dot(port_vector, I)
        except np.linalg.LinAlgError:
            Z[i] = np.nan

    return Z


def test_prima_frequency_response():
    """Test that PRIMA I transformation preserves frequency response."""
    print("=" * 70)
    print("Verify PRIMA I Frequency Response Matches Original PEEC")
    print("=" * 70)

    # Create test PEEC matrices (dense, symmetric)
    n = 10  # Number of loops

    # Dense inductance matrix with mutual coupling
    # L_ij = L0 * exp(-|i-j| / decay)
    L0 = 1e-6  # 1 uH base inductance
    decay = 3.0
    L_dense = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            L_dense[i, j] = L0 * np.exp(-abs(i - j) / decay)

    # Diagonal resistance matrix
    R0 = 0.01  # 10 mOhm
    R_diag_orig = np.diag(np.full(n, R0))

    print(f"\nOriginal PEEC matrices:")
    print(f"  L_dense: {n}x{n} dense symmetric (mutual inductance)")
    print(f"  R_diag:  {n}x{n} diagonal")
    print(f"  L0 = {L0*1e6:.2f} uH, R0 = {R0*1e3:.1f} mOhm")

    # Apply PRIMA I transformation
    # lanczos(K=R, N=L) -> R_diag (diagonal), L_tridiag (tridiagonal)
    print("\n--- Applying PRIMA I Transformation ---")
    print("lanczos(K=R, N=L) -> R_diag, L_tridiag")

    result = lanczos(R_diag_orig, L_dense)

    R_prima = result.R_diag
    L_prima = result.L_tridiag
    U = result.U
    V = result.V

    print(f"\nPRIMA I transformed matrices:")
    print(f"  R_diag:    {R_prima.shape[0]}x{R_prima.shape[1]} diagonal")
    print(f"  L_tridiag: {L_prima.shape[0]}x{L_prima.shape[1]} tridiagonal")
    print(f"  n_output = {result.n_output}")

    # Verify structure (use relative tolerance)
    R_norm = np.linalg.norm(R_prima, 'fro')
    L_norm = np.linalg.norm(L_prima, 'fro')
    is_diag_R = np.allclose(R_prima, np.diag(np.diag(R_prima)), atol=1e-10 * R_norm)
    is_tridiag_L = True
    for i in range(L_prima.shape[0]):
        for j in range(L_prima.shape[1]):
            if abs(i - j) > 1 and abs(L_prima[i, j]) > 1e-10 * L_norm:
                is_tridiag_L = False
                break

    print(f"\n  R_prima is diagonal: {is_diag_R}")
    print(f"  L_prima is tridiagonal: {is_tridiag_L}")

    # Frequency sweep
    frequencies = np.logspace(1, 7, 100)  # 10 Hz to 10 MHz

    # Port vector (excitation at first node)
    port_orig = np.zeros(n)
    port_orig[0] = 1.0

    # Transform port vector for PRIMA: port_prima = U^T * port_orig
    port_prima = U.T @ port_orig

    print(f"\n--- Computing Frequency Response ---")
    print(f"  Frequency range: {frequencies[0]:.0f} Hz to {frequencies[-1]/1e6:.1f} MHz")
    print(f"  Port: node 0")

    # Compute impedance for original PEEC
    Z_orig = compute_impedance(R_diag_orig, L_dense, frequencies, port_orig)

    # Compute impedance for PRIMA I
    Z_prima = compute_impedance(R_prima, L_prima, frequencies, port_prima)

    # Compare
    rel_error = np.abs(Z_prima - Z_orig) / np.abs(Z_orig)
    max_rel_error = np.nanmax(rel_error)
    mean_rel_error = np.nanmean(rel_error)

    print(f"\n--- Comparison Results ---")
    print(f"  Max relative error:  {max_rel_error:.2e}")
    print(f"  Mean relative error: {mean_rel_error:.2e}")

    if max_rel_error < 1e-10:
        print("\n=> PASS: PRIMA I frequency response matches original PEEC!")
    else:
        print(f"\n=> WARNING: Relative error {max_rel_error:.2e} > 1e-10")

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Magnitude
    ax1 = axes[0, 0]
    ax1.loglog(frequencies, np.abs(Z_orig), 'b-', label='Original PEEC', linewidth=2)
    ax1.loglog(frequencies, np.abs(Z_prima), 'r--', label='PRIMA I', linewidth=2)
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Phase
    ax2 = axes[0, 1]
    ax2.semilogx(frequencies, np.angle(Z_orig, deg=True), 'b-', label='Original PEEC', linewidth=2)
    ax2.semilogx(frequencies, np.angle(Z_prima, deg=True), 'r--', label='PRIMA I', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Phase [deg]')
    ax2.set_title('Impedance Phase')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Real part
    ax3 = axes[1, 0]
    ax3.semilogx(frequencies, np.real(Z_orig), 'b-', label='Original PEEC', linewidth=2)
    ax3.semilogx(frequencies, np.real(Z_prima), 'r--', label='PRIMA I', linewidth=2)
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z) [Ohm]')
    ax3.set_title('Resistance (Real Part)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Relative error
    ax4 = axes[1, 1]
    ax4.loglog(frequencies, rel_error, 'k-', linewidth=1)
    ax4.axhline(y=1e-10, color='g', linestyle='--', label='1e-10 threshold')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Relative Error')
    ax4.set_title('Relative Error |Z_prima - Z_orig| / |Z_orig|')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('prima_frequency_response_verification.png', dpi=150)
    print(f"\nPlot saved: prima_frequency_response_verification.png")
    plt.close()

    return max_rel_error < 1e-10


def test_prima_with_surface_impedance():
    """Test that PRIMA I transformation preserves frequency response with Z_s."""
    print("\n" + "=" * 70)
    print("Test 2: Verify PRIMA I Frequency Response with Surface Impedance Z_s")
    print("=" * 70)

    # Create test PEEC matrices (dense, symmetric)
    n = 10  # Number of loops

    # Dense inductance matrix with mutual coupling
    L0 = 1e-6  # 1 uH base inductance
    decay = 3.0
    L_dense = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            L_dense[i, j] = L0 * np.exp(-abs(i - j) / decay)

    # Diagonal resistance matrix
    R0 = 0.01  # 10 mOhm
    R_diag_orig = np.diag(np.full(n, R0))

    print(f"\nOriginal PEEC matrices:")
    print(f"  L_dense: {n}x{n} dense symmetric")
    print(f"  R_diag:  {n}x{n} diagonal")
    print(f"  L0 = {L0*1e6:.2f} uH, R0 = {R0*1e3:.1f} mOhm")

    # Surface impedance parameters (copper conductor)
    sigma = 5.8e7  # S/m (copper conductivity)
    mu_r = 1.0     # relative permeability

    print(f"\nSurface impedance parameters:")
    print(f"  sigma = {sigma:.2e} S/m (copper)")
    print(f"  mu_r = {mu_r}")

    # Apply PRIMA I transformation
    print("\n--- Applying PRIMA I Transformation ---")
    result = lanczos(R_diag_orig, L_dense)

    R_prima = result.R_diag
    L_prima = result.L_tridiag
    U = result.U
    V = result.V

    print(f"  R_diag:    {R_prima.shape[0]}x{R_prima.shape[1]} diagonal")
    print(f"  L_tridiag: {L_prima.shape[0]}x{L_prima.shape[1]} tridiagonal")

    # Frequency sweep
    frequencies = np.logspace(1, 7, 100)  # 10 Hz to 10 MHz

    # Compute surface impedance Z_s(f)
    Z_s = compute_surface_impedance(frequencies, sigma, mu_r)

    # Show Z_s at a few frequencies
    print(f"\n--- Surface Impedance Z_s(f) ---")
    for f_idx in [0, 50, 99]:
        f = frequencies[f_idx]
        delta = np.sqrt(2 / (2*np.pi*f * 4*np.pi*1e-7 * mu_r * sigma)) if f > 0 else np.inf
        print(f"  f={f:.0e} Hz: Z_s={Z_s[f_idx]:.4e}, delta={delta*1e3:.3f} mm")

    # Port vectors
    port_orig = np.zeros(n)
    port_orig[0] = 1.0
    port_prima = U.T @ port_orig

    print(f"\n--- Computing Frequency Response with Z_s ---")
    print(f"  Frequency range: {frequencies[0]:.0f} Hz to {frequencies[-1]/1e6:.1f} MHz")

    # Compute impedance for original PEEC with Z_s
    Z_orig = compute_impedance(R_diag_orig, L_dense, frequencies, port_orig, Z_s)

    # For PRIMA I with Z_s, we need to transform the identity matrix:
    # In original: Z_matrix = R + sL + Z_s*I
    # In PRIMA I:    Z_matrix = R' + sL' + Z_s*I'
    # where I' = U^H * I * V (NOT identity in general)
    I_transformed = U.T @ np.eye(n) @ V
    print(f"\n  I_transformed = U^T * I * V:")
    print(f"    Is I' = I? {np.allclose(I_transformed, np.eye(n))}")
    print(f"    ||I' - I||_F = {np.linalg.norm(I_transformed - np.eye(n), 'fro'):.6e}")

    # Compute impedance for PRIMA I with Z_s (using transformed identity)
    Z_prima = compute_impedance_with_transformed_Zs(R_prima, L_prima, frequencies, port_prima, Z_s, I_transformed)

    # Compare
    rel_error = np.abs(Z_prima - Z_orig) / np.abs(Z_orig)
    max_rel_error = np.nanmax(rel_error)
    mean_rel_error = np.nanmean(rel_error)

    print(f"\n--- Comparison Results (with Z_s) ---")
    print(f"  Max relative error:  {max_rel_error:.2e}")
    print(f"  Mean relative error: {mean_rel_error:.2e}")

    if max_rel_error < 1e-10:
        print("\n=> PASS: PRIMA I with Z_s matches original PEEC + Z_s!")
    else:
        print(f"\n=> WARNING: Relative error {max_rel_error:.2e} > 1e-10")

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('PRIMA I Frequency Response with Surface Impedance Z_s', fontsize=14)

    # Magnitude
    ax1 = axes[0, 0]
    ax1.loglog(frequencies, np.abs(Z_orig), 'b-', label='Original PEEC + Z_s', linewidth=2)
    ax1.loglog(frequencies, np.abs(Z_prima), 'r--', label='PRIMA I + Z_s', linewidth=2)
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Phase
    ax2 = axes[0, 1]
    ax2.semilogx(frequencies, np.angle(Z_orig, deg=True), 'b-', label='Original PEEC + Z_s', linewidth=2)
    ax2.semilogx(frequencies, np.angle(Z_prima, deg=True), 'r--', label='PRIMA I + Z_s', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Phase [deg]')
    ax2.set_title('Impedance Phase')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Real part (shows skin effect resistance)
    ax3 = axes[1, 0]
    ax3.loglog(frequencies, np.real(Z_orig), 'b-', label='Original PEEC + Z_s', linewidth=2)
    ax3.loglog(frequencies, np.real(Z_prima), 'r--', label='PRIMA I + Z_s', linewidth=2)
    ax3.loglog(frequencies, np.real(Z_s), 'g:', label='Z_s (skin effect)', linewidth=1)
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z) [Ohm]')
    ax3.set_title('Resistance (Real Part)')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Relative error
    ax4 = axes[1, 1]
    ax4.loglog(frequencies, rel_error, 'k-', linewidth=1)
    ax4.axhline(y=1e-10, color='g', linestyle='--', label='1e-10 threshold')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Relative Error')
    ax4.set_title('Relative Error |Z_prima - Z_orig| / |Z_orig|')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('prima_frequency_response_with_zs.png', dpi=150)
    print(f"\nPlot saved: prima_frequency_response_with_zs.png")
    plt.close()

    return max_rel_error < 1e-10


if __name__ == '__main__':
    # Test 1: Without surface impedance
    success1 = test_prima_frequency_response()

    # Test 2: With surface impedance Z_s
    success2 = test_prima_with_surface_impedance()

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Test 1 (without Z_s): {'PASS' if success1 else 'FAIL'}")
    print(f"  Test 2 (with Z_s):    {'PASS' if success2 else 'FAIL'}")
    print("=" * 70)

    sys.exit(0 if (success1 and success2) else 1)
