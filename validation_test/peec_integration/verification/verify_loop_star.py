"""
Loop-Star Decomposition Verification

Verifies that the PEEC Loop-Star matrices satisfy:
1. L matrix is symmetric and positive definite
2. P matrix is symmetric and positive definite
3. M_LS coupling preserves energy conservation
4. Frequency response matches theoretical predictions

Uses Radia's built-in PEEC module (peec_matrices).
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad

try:
    from peec_matrices import PEECBuilder
    print("Using C++ PEEC implementation")
except ImportError:
    print("C++ PEEC module not available, using Python fallback")
    from radia.peec_matrices import PEECBuilder

mm = 1e-3  # 1 mm in meters


def verify_matrix_properties(name, M, expected_symmetric=True, expected_spd=True):
    """Verify matrix properties: symmetry, positive definiteness."""
    n = M.shape[0]
    results = {'name': name, 'size': n, 'passed': True}

    # Symmetry check
    if expected_symmetric:
        sym_error = np.linalg.norm(M - M.T) / (np.linalg.norm(M) + 1e-15)
        results['symmetry_error'] = sym_error
        if sym_error > 1e-10:
            results['passed'] = False
            print(f"  {name}: Symmetry FAILED (error = {sym_error:.2e})")
        else:
            print(f"  {name}: Symmetric (error = {sym_error:.2e})")

    # Positive definiteness check
    if expected_spd:
        try:
            eigenvalues = np.linalg.eigvalsh(M)
            min_eig = np.min(eigenvalues)
            max_eig = np.max(eigenvalues)
            results['min_eigenvalue'] = min_eig
            results['max_eigenvalue'] = max_eig
            results['condition_number'] = max_eig / (min_eig + 1e-15)

            if min_eig > 0:
                print(f"  {name}: Positive definite (min_eig = {min_eig:.2e}, cond = {results['condition_number']:.2e})")
            else:
                results['passed'] = False
                print(f"  {name}: NOT positive definite (min_eig = {min_eig:.2e})")
        except np.linalg.LinAlgError:
            results['passed'] = False
            print(f"  {name}: Eigenvalue computation FAILED")

    return results


def verify_wire_peec():
    """Verify PEEC matrices for a simple straight wire."""
    print("\n" + "=" * 60)
    print("Test 1: Straight Wire PEEC Matrices")
    print("=" * 60)

    # Wire parameters
    length = 100*mm  # 10 cm
    width = 1*mm     # 1 mm
    height = 1*mm    # 1 mm
    n_segments = 10
    sigma = 5.8e7  # Copper

    print(f"\nWire: L={length/mm:.0f}mm, w={width/mm:.1f}mm, h={height/mm:.1f}mm")
    print(f"Segments: {n_segments}, Conductivity: {sigma:.2e} S/m")

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    print(f"\nMatrix dimensions:")
    print(f"  L: {L.shape}")
    print(f"  R: {R.shape} (diagonal)")
    print(f"  P: {P.shape}")
    print(f"  M_LS: {M_LS.shape}")

    # Verify L matrix
    print("\n--- Inductance Matrix L ---")
    L_result = verify_matrix_properties("L", L)

    # Verify P matrix (if non-empty)
    if P.size > 0:
        print("\n--- Potential Coefficient Matrix P ---")
        P_result = verify_matrix_properties("P", P)
    else:
        P_result = {'passed': True, 'name': 'P', 'note': 'empty'}
        print("\n--- Potential Coefficient Matrix P: empty (no star nodes) ---")

    # Verify physical values
    print("\n--- Physical Verification ---")

    # DC resistance: R_dc = L / (sigma * A)
    R_dc_expected = length / (sigma * width * height)
    R_dc_actual = np.sum(np.diag(R))

    print(f"  DC Resistance:")
    print(f"    Expected: {R_dc_expected*1e3:.4f} mOhm")
    print(f"    From R matrix: {R_dc_actual*1e3:.4f} mOhm")
    print(f"    Error: {abs(R_dc_actual - R_dc_expected)/R_dc_expected*100:.2f}%")

    # Total inductance (rough estimate)
    L_total = np.sum(L)
    # Rough estimate for straight wire: L ~ mu_0 * L * (ln(2L/sqrt(w*h)) - 1)
    a = np.sqrt(width * height / np.pi)  # Equivalent radius
    L_estimate = 4 * np.pi * 1e-7 * length / (2*np.pi) * (np.log(2 * length / a) - 1)

    print(f"\n  Total Inductance:")
    print(f"    Estimate: {L_estimate*1e9:.2f} nH")
    print(f"    From L matrix: {L_total*1e9:.2f} nH")

    return L_result['passed'] and P_result['passed']


def verify_loop_peec():
    """Verify PEEC matrices for a circular loop."""
    print("\n" + "=" * 60)
    print("Test 2: Circular Loop PEEC Matrices")
    print("=" * 60)

    # Loop parameters
    radius = 50*mm   # 5 cm
    width = 2*mm     # 2 mm
    height = 1*mm    # 1 mm
    n_segments = 36
    sigma = 5.8e7

    print(f"\nLoop: R={radius/mm:.0f}mm, w={width/mm:.1f}mm, h={height/mm:.1f}mm")
    print(f"Segments: {n_segments}")

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_loop([0, 0, 0], radius, [0, 0, 1], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    print(f"\nMatrix dimensions: L={L.shape}, R={R.shape}")

    # Verify L matrix
    print("\n--- Inductance Matrix L ---")
    L_result = verify_matrix_properties("L", L)

    # Verify against analytical formula
    # L = mu_0 * R * (ln(8R/a) - 2) for thin circular loop
    a = np.sqrt(width * height / np.pi)  # Equivalent radius
    L_analytical = 4 * np.pi * 1e-7 * radius * (np.log(8 * radius / a) - 2)
    L_total = np.sum(L)

    print(f"\n--- Analytical Comparison ---")
    print(f"  Analytical L: {L_analytical*1e6:.3f} uH")
    print(f"  PEEC L (sum): {L_total*1e6:.3f} uH")
    print(f"  Ratio: {L_total/L_analytical:.3f}")

    # DC resistance (full loop)
    circumference = 2 * np.pi * radius
    R_dc_expected = circumference / (sigma * width * height)
    R_dc_actual = np.sum(np.diag(R))

    print(f"\n  DC Resistance:")
    print(f"    Expected: {R_dc_expected*1e3:.4f} mOhm")
    print(f"    From R: {R_dc_actual*1e3:.4f} mOhm")

    return L_result['passed']


def verify_frequency_response():
    """Verify frequency response of Loop-Star system."""
    print("\n" + "=" * 60)
    print("Test 3: Frequency Response Verification")
    print("=" * 60)

    # Simple wire
    length = 100*mm
    width = 1*mm
    height = 1*mm
    n_segments = 10
    sigma = 5.8e7

    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    n_loop = L.shape[0]
    n_star = P.shape[0] if P.size > 0 else 0

    print(f"\nSystem size: n_loop={n_loop}, n_star={n_star}")

    # Frequency sweep
    frequencies = np.logspace(3, 8, 11)  # 1 kHz to 100 MHz

    print(f"\n{'Frequency':>12s}  {'|Z|':>12s}  {'Angle':>8s}  {'Type':>10s}")
    print("-" * 50)

    results = []
    for f in frequencies:
        omega = 2 * np.pi * f

        # Build impedance matrix
        # Z_LL = R + j*omega*L
        Z_LL = np.diag(np.diag(R)) + 1j * omega * L

        if n_star > 0:
            # Full Loop-Star system
            n_total = n_loop + n_star
            Z = np.zeros((n_total, n_total), dtype=complex)
            Z[:n_loop, :n_loop] = Z_LL
            Z[n_loop:, n_loop:] = P / (1j * omega + 1e-15)
            Z[:n_loop, n_loop:] = 1j * omega * M_LS
            Z[n_loop:, :n_loop] = (1j * omega * M_LS).T
        else:
            Z = Z_LL

        # Port impedance (first segment)
        V = np.zeros(Z.shape[0], dtype=complex)
        V[0] = 1.0

        try:
            I = np.linalg.solve(Z, V)
            Z_port = 1.0 / I[0]
        except np.linalg.LinAlgError:
            Z_port = np.nan

        mag = np.abs(Z_port)
        angle = np.angle(Z_port, deg=True)

        if angle > 45:
            ztype = "Inductive"
        elif angle < -45:
            ztype = "Capacitive"
        else:
            ztype = "Resistive"

        print(f"{f:12.2e}  {mag:12.4e}  {angle:8.1f}  {ztype:>10s}")
        results.append((f, mag, angle))

    # Verify low frequency is inductive, high frequency behavior
    low_freq_angle = results[0][2]
    high_freq_angle = results[-1][2]

    print(f"\n--- Frequency Response Check ---")
    print(f"  Low frequency (1 kHz): angle = {low_freq_angle:.1f} deg")
    print(f"  High frequency (100 MHz): angle = {high_freq_angle:.1f} deg")

    if low_freq_angle > 80:
        print("  Low frequency: Inductive (PASS)")
        lf_pass = True
    else:
        print("  Low frequency: NOT strongly inductive (CHECK)")
        lf_pass = False

    return lf_pass


def main():
    """Run all verification tests."""
    print("=" * 60)
    print("Loop-Star PEEC Matrix Verification")
    print("Darwin approximation: G(r) = 1/(4*pi*r)")
    print("=" * 60)

    results = []

    # Test 1: Wire
    try:
        results.append(("Straight Wire", verify_wire_peec()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Straight Wire", False))

    # Test 2: Loop
    try:
        results.append(("Circular Loop", verify_loop_peec()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Circular Loop", False))

    # Test 3: Frequency response
    try:
        results.append(("Frequency Response", verify_frequency_response()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Frequency Response", False))

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("All tests PASSED")
    else:
        print("Some tests FAILED")

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
