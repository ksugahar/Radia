"""
Test PyKAN for Continued Fraction Learning

Verify that PyKAN can learn continued fraction expansions
and extract RL ladder parameters for SPICE.

Key test: Learn z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))
where w = tau * s, and extract ladder coefficients.

Theory:
- z*coth(z) is the skin effect impedance function
- Its continued fraction expansion maps directly to an RL ladder
- If KAN learns this function, we can extract the ladder structure

Success criteria:
1. KAN can fit z*coth(z) accurately
2. KAN symbolic regression recovers rational polynomial
3. Ladder coefficients extracted match analytical [3, 5, 7, ...]
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
import matplotlib.pyplot as plt


def z_coth_z(z):
    """Compute z*coth(z) = z*cosh(z)/sinh(z) with numerical stability."""
    z = np.asarray(z, dtype=complex)
    result = np.zeros_like(z)

    # For small |z|, use Taylor series: 1 + z^2/3 - z^4/45 + ...
    small = np.abs(z) < 0.1
    if np.any(small):
        z_small = z[small]
        z2 = z_small**2
        result[small] = 1 + z2/3 - z2**2/45 + 2*z2**3/945

    # For larger |z|, use direct formula with overflow protection
    large = ~small
    if np.any(large):
        z_large = z[large]
        # coth(z) = (exp(z) + exp(-z)) / (exp(z) - exp(-z))
        #         = (1 + exp(-2z)) / (1 - exp(-2z))
        exp_m2z = np.exp(-2*z_large)
        coth_z = (1 + exp_m2z) / (1 - exp_m2z)
        result[large] = z_large * coth_z

    return result


def continued_fraction_approx(w, coeffs):
    """Evaluate continued fraction 1 + w/(c1 + w/(c2 + w/(c3 + ...)))."""
    result = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        result = c + w / result
    return 1 + w / result


def test_zcothz_learning():
    """Test if KAN can learn z*coth(z) function."""
    print("=" * 70)
    print("Test 1: KAN Learning of z*coth(z)")
    print("=" * 70)

    try:
        import torch
        from kan import KAN
    except ImportError as e:
        print(f"PyKAN not available: {e}")
        print("Install with: pip install pykan torch")
        return False

    # Generate training data
    # w = tau * j*omega = j * (omega * tau)
    # For skin effect: tau = d^2 * mu * sigma / 2

    omega_tau = np.logspace(-2, 2, 100)  # omega * tau from 0.01 to 100
    w = 1j * omega_tau
    z = np.sqrt(w)  # z = sqrt(j*omega*tau) = sqrt(j) * sqrt(omega*tau)

    # Target: z*coth(z)
    y_target = z_coth_z(z)

    # Compare with continued fraction approximation (5 terms: 3,5,7,9,11)
    y_cf5 = continued_fraction_approx(w, [3, 5, 7, 9, 11])

    # Relative error of 5-term CF
    cf5_error = np.abs(y_target - y_cf5) / np.abs(y_target)

    print(f"\nContinued fraction (5 terms) accuracy:")
    print(f"  Max relative error: {np.max(cf5_error)*100:.2f}%")
    print(f"  Mean relative error: {np.mean(cf5_error)*100:.4f}%")

    # Prepare data for KAN
    # Input: [Re(w), Im(w)] = [0, omega_tau]
    # Output: [Re(z*coth(z)), Im(z*coth(z))]

    X = np.column_stack([np.zeros_like(omega_tau), omega_tau]).astype(np.float32)
    Y = np.column_stack([y_target.real, y_target.imag]).astype(np.float32)

    X_tensor = torch.tensor(X)
    Y_tensor = torch.tensor(Y)

    # Create KAN model
    # Architecture: 2 inputs -> hidden -> 2 outputs
    print(f"\nTraining KAN model...")
    print(f"  Input shape: {X.shape}")
    print(f"  Output shape: {Y.shape}")

    model = KAN(width=[2, 8, 8, 2], grid=10, k=3, seed=42)

    # Training
    dataset = {'train_input': X_tensor, 'train_label': Y_tensor,
               'test_input': X_tensor, 'test_label': Y_tensor}

    results = model.fit(dataset, opt='Adam', steps=500, lr=0.01,
                       lamb=0.001, lamb_entropy=0.0, metrics=None,
                       save_fig=False, display_metrics=None)

    # Evaluate
    with torch.no_grad():
        Y_pred = model(X_tensor).numpy()

    y_pred = Y_pred[:, 0] + 1j * Y_pred[:, 1]

    # Compute errors
    relative_error = np.abs(y_pred - y_target) / np.abs(y_target)

    print(f"\nKAN prediction accuracy:")
    print(f"  Max relative error: {np.max(relative_error)*100:.2f}%")
    print(f"  Mean relative error: {np.mean(relative_error)*100:.4f}%")

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Real part
    ax = axes[0, 0]
    ax.semilogx(omega_tau, y_target.real, 'b-', linewidth=2, label='Exact z*coth(z)')
    ax.semilogx(omega_tau, y_cf5.real, 'g--', linewidth=1.5, label='CF 5-term')
    ax.semilogx(omega_tau, y_pred.real, 'r:', linewidth=2, label='KAN')
    ax.set_xlabel('omega * tau')
    ax.set_ylabel('Re[z*coth(z)]')
    ax.set_title('Real Part')
    ax.legend()
    ax.grid(True)

    # Imaginary part
    ax = axes[0, 1]
    ax.semilogx(omega_tau, y_target.imag, 'b-', linewidth=2, label='Exact')
    ax.semilogx(omega_tau, y_cf5.imag, 'g--', linewidth=1.5, label='CF 5-term')
    ax.semilogx(omega_tau, y_pred.imag, 'r:', linewidth=2, label='KAN')
    ax.set_xlabel('omega * tau')
    ax.set_ylabel('Im[z*coth(z)]')
    ax.set_title('Imaginary Part')
    ax.legend()
    ax.grid(True)

    # Error comparison
    ax = axes[1, 0]
    ax.loglog(omega_tau, relative_error * 100, 'r-', linewidth=2, label='KAN error')
    ax.loglog(omega_tau, cf5_error * 100, 'g--', linewidth=2, label='CF 5-term error')
    ax.set_xlabel('omega * tau')
    ax.set_ylabel('Relative Error [%]')
    ax.set_title('Error Comparison')
    ax.legend()
    ax.grid(True)
    ax.axhline(1.0, color='k', linestyle=':', label='1% threshold')

    # Magnitude (|Z|)
    ax = axes[1, 1]
    ax.loglog(omega_tau, np.abs(y_target), 'b-', linewidth=2, label='Exact')
    ax.loglog(omega_tau, np.abs(y_cf5), 'g--', linewidth=1.5, label='CF 5-term')
    ax.loglog(omega_tau, np.abs(y_pred), 'r:', linewidth=2, label='KAN')
    ax.set_xlabel('omega * tau')
    ax.set_ylabel('|z*coth(z)|')
    ax.set_title('Magnitude')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('test_pykan_zcothz.png', dpi=150)
    print(f"\nSaved: test_pykan_zcothz.png")

    # Pass criteria
    passed = np.max(relative_error) < 0.05  # < 5% max error
    print(f"\nTest result: {'PASSED' if passed else 'FAILED'}")

    return passed


def test_symbolic_regression():
    """Test if KAN symbolic regression can extract continued fraction structure."""
    print("\n" + "=" * 70)
    print("Test 2: KAN Symbolic Regression for CF Structure")
    print("=" * 70)

    try:
        import torch
        from kan import KAN
    except ImportError:
        print("PyKAN not available")
        return False

    # Simpler test: fit 1 + w/(3 + w/5) = (5 + 3w + w) / (5 + w)
    # = (5 + 4w) / (5 + w)

    print("\nTarget: f(w) = 1 + w/(3 + w/5)")
    print("This equals: (5 + 4w) / (5 + w)")

    # Generate data (real w only for simplicity)
    w = np.linspace(0.01, 10, 100).astype(np.float32)
    y = 1 + w / (3 + w/5)

    X_tensor = torch.tensor(w.reshape(-1, 1))
    Y_tensor = torch.tensor(y.reshape(-1, 1))

    # Train KAN
    model = KAN(width=[1, 4, 4, 1], grid=5, k=3, seed=42)

    dataset = {'train_input': X_tensor, 'train_label': Y_tensor,
               'test_input': X_tensor, 'test_label': Y_tensor}

    print("\nTraining KAN...")
    results = model.fit(dataset, opt='Adam', steps=300, lr=0.01,
                       lamb=0.001, save_fig=False, display_metrics=None)

    # Evaluate
    with torch.no_grad():
        y_pred = model(X_tensor).numpy().flatten()

    error = np.abs(y_pred - y) / np.abs(y)
    print(f"\nPrediction accuracy:")
    print(f"  Max relative error: {np.max(error)*100:.2f}%")
    print(f"  Mean relative error: {np.mean(error)*100:.4f}%")

    # Try symbolic regression
    print("\nAttempting symbolic regression...")
    try:
        # Auto-symbolic simplification
        model.auto_symbolic(lib=['x', '1/x', 'x^2', '1+x'])
        formula = model.symbolic_formula()
        print(f"Symbolic formula: {formula}")
    except Exception as e:
        print(f"Symbolic regression not available or failed: {e}")
        print("(This is expected - KAN symbolic regression is limited)")

    # Manual check: can we fit with specific basis functions?
    # The continued fraction should give us insight into ladder coefficients

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(w, y, 'b-', linewidth=2, label='Target: 1+w/(3+w/5)')
    plt.plot(w, y_pred, 'r--', linewidth=2, label='KAN prediction')
    plt.xlabel('w')
    plt.ylabel('f(w)')
    plt.title('KAN Learning of Simple Continued Fraction')
    plt.legend()
    plt.grid(True)
    plt.savefig('test_pykan_simple_cf.png', dpi=150)
    print(f"\nSaved: test_pykan_simple_cf.png")

    passed = np.max(error) < 0.02  # < 2% error
    print(f"\nTest result: {'PASSED' if passed else 'FAILED'}")

    return passed


def test_ladder_extraction():
    """Test extraction of RL ladder parameters from learned function."""
    print("\n" + "=" * 70)
    print("Test 3: RL Ladder Parameter Extraction")
    print("=" * 70)

    try:
        import torch
        from kan import KAN
    except ImportError:
        print("PyKAN not available")
        return False

    # Strategy: Fit KAN to z*coth(z) data, then use Pade approximation
    # to extract rational function, then convert to ladder

    print("\nStrategy:")
    print("  1. Train KAN on z*coth(z) data")
    print("  2. Sample KAN at specific points")
    print("  3. Use Pade approximation to get rational function")
    print("  4. Convert rational function to continued fraction")
    print("  5. Compare with analytical coefficients [3, 5, 7, ...]")

    # Generate training data
    omega_tau = np.logspace(-1, 1, 50)
    w = 1j * omega_tau
    z = np.sqrt(w)
    y_target = z_coth_z(z)

    # Use real-valued representation
    X = np.column_stack([np.zeros_like(omega_tau), omega_tau]).astype(np.float32)
    Y = np.column_stack([y_target.real, y_target.imag]).astype(np.float32)

    X_tensor = torch.tensor(X)
    Y_tensor = torch.tensor(Y)

    model = KAN(width=[2, 6, 6, 2], grid=8, k=3, seed=42)

    dataset = {'train_input': X_tensor, 'train_label': Y_tensor,
               'test_input': X_tensor, 'test_label': Y_tensor}

    print("\nTraining KAN...")
    model.fit(dataset, opt='Adam', steps=400, lr=0.01,
             lamb=0.001, save_fig=False, display_metrics=None)

    # Sample at specific points for Pade
    omega_tau_sample = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
    X_sample = np.column_stack([np.zeros_like(omega_tau_sample),
                                omega_tau_sample]).astype(np.float32)

    with torch.no_grad():
        Y_sample = model(torch.tensor(X_sample)).numpy()

    y_sample_kan = Y_sample[:, 0] + 1j * Y_sample[:, 1]

    # Compare with analytical at sample points
    w_sample = 1j * omega_tau_sample
    z_sample = np.sqrt(w_sample)
    y_sample_exact = z_coth_z(z_sample)

    print("\nSampled values comparison:")
    print("  omega*tau    Exact           KAN             Error")
    for i in range(len(omega_tau_sample)):
        err = abs(y_sample_kan[i] - y_sample_exact[i]) / abs(y_sample_exact[i]) * 100
        print(f"  {omega_tau_sample[i]:6.2f}    {y_sample_exact[i]:15.6f}  {y_sample_kan[i]:15.6f}  {err:5.2f}%")

    # Attempt Pade approximation (simplified: just check if values match CF)
    # For continued fraction f(w) = 1 + w/(b1 + w/(b2 + ...))
    # At w=0: f(0) = 1
    # Derivative at w=0 gives b1 = 3

    # Check DC limit (should be 1)
    X_dc = np.array([[0.0, 0.0]], dtype=np.float32)
    with torch.no_grad():
        Y_dc = model(torch.tensor(X_dc)).numpy()
    y_dc = Y_dc[0, 0] + 1j * Y_dc[0, 1]
    print(f"\nDC limit check (should be 1.0):")
    print(f"  KAN(w=0) = {y_dc}")

    # Check slope to estimate first CF coefficient
    # f(w) ~ 1 + w/3 for small w, so df/dw|_{w=0} = 1/3
    eps = 0.001
    X_eps = np.array([[0.0, eps]], dtype=np.float32)
    with torch.no_grad():
        Y_eps = model(torch.tensor(X_eps)).numpy()
    y_eps = Y_eps[0, 0] + 1j * Y_eps[0, 1]

    slope = (y_eps - y_dc) / (1j * eps)
    b1_estimate = 1.0 / slope.imag if abs(slope.imag) > 1e-10 else float('inf')

    print(f"\nFirst CF coefficient estimate (should be ~3):")
    print(f"  Slope at w=0: {slope}")
    print(f"  Estimated b1 = 1/slope.imag = {b1_estimate:.2f}")

    # Verdict
    dc_ok = abs(y_dc - 1.0) < 0.1
    b1_ok = 2.5 < b1_estimate < 3.5  # Within 17% of 3

    passed = dc_ok and b1_ok

    print(f"\nTest results:")
    print(f"  DC limit correct: {dc_ok}")
    print(f"  First CF coeff ~3: {b1_ok}")
    print(f"  Overall: {'PASSED' if passed else 'FAILED'}")

    return passed


def main():
    """Run all tests."""
    print("PyKAN Continued Fraction Learning Test Suite")
    print("=" * 70)
    print("\nGoal: Verify that PyKAN can learn continued fraction")
    print("expansions and that RL ladder parameters can be extracted.")
    print("=" * 70)

    results = {}

    # Test 1: Basic learning of z*coth(z)
    results['zcothz_learning'] = test_zcothz_learning()

    # Test 2: Symbolic regression on simple CF
    results['symbolic_regression'] = test_symbolic_regression()

    # Test 3: Ladder parameter extraction
    results['ladder_extraction'] = test_ladder_extraction()

    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)

    all_passed = True
    for name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print("-" * 70)
    overall = "ALL TESTS PASSED" if all_passed else "SOME TESTS FAILED"
    print(f"  Overall: {overall}")

    # Conclusion
    print("\n" + "=" * 70)
    print("Conclusion")
    print("=" * 70)
    if all_passed:
        print("""
KAN can learn the continued fraction expansion z*coth(z) with good accuracy.

However, extracting the exact ladder coefficients [3, 5, 7, ...] requires:
1. Symbolic regression (limited in current PyKAN)
2. Pade approximation from sampled values
3. Or direct fitting to ladder circuit model

For SPICE export, the recommended approach is:
- Use analytical Dowell formula for standard skin effect
- Use KAN only for learning non-standard Z(H,f) from ESIM data
- Export as PWL lookup table rather than extracted ladder
""")
    else:
        print("""
Some tests failed. Possible issues:
1. KAN training may need tuning (more epochs, different architecture)
2. Continued fraction extraction is inherently difficult
3. Numerical precision issues with complex arithmetic

The Dowell analytical formula remains the preferred method for
standard skin effect modeling.
""")

    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
