"""
CF-KAN: Continued Fraction KAN for Skin Effect

Tests whether a neural network with continued fraction structure
can learn the Dowell coefficients [3, 5, 7, 9, 11] from z*coth(z) data.

Key insight: Instead of using standard KAN (sum of compositions),
we use a network that directly outputs continued fraction coefficients.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim


class CFKAN(nn.Module):
    """
    Continued Fraction KAN

    Learns coefficients a_1, a_2, ..., a_n such that:

    f(z) = 1 + z^2 / (a_1 + z^2 / (a_2 + z^2 / (a_3 + ...)))

    For z*coth(z), the exact coefficients are [3, 5, 7, 9, 11, ...]
    """

    def __init__(self, n_stages=5, init_mode='random'):
        super().__init__()
        self.n_stages = n_stages

        # Initialize coefficients
        if init_mode == 'random':
            # Random initialization around expected values
            init_vals = torch.rand(n_stages) * 10 + 1  # 1 to 11
        elif init_mode == 'ones':
            init_vals = torch.ones(n_stages)
        elif init_mode == 'near_dowell':
            # Initialize near true values (for debugging)
            dowell = torch.tensor([3.0, 5.0, 7.0, 9.0, 11.0])[:n_stages]
            init_vals = dowell + torch.randn(n_stages) * 0.5
        else:
            init_vals = torch.ones(n_stages) * 5.0

        self.a = nn.Parameter(init_vals)

    def forward(self, z):
        """
        Evaluate continued fraction: 1 + z^2 / (a_1 + z^2 / (a_2 + ...))

        z can be complex (torch.complex64 or complex128)
        """
        z_sq = z ** 2

        # Evaluate from the innermost term outward
        # Start with the last coefficient
        result = self.a[-1].to(z.dtype)

        # Work backwards
        for i in range(self.n_stages - 2, -1, -1):
            result = self.a[i].to(z.dtype) + z_sq / result

        # Final: 1 + z^2 / result
        return 1.0 + z_sq / result

    def get_coefficients(self):
        """Return learned coefficients as numpy array."""
        return self.a.detach().cpu().numpy()


def zcothz_exact(z):
    """Exact z*coth(z) using numpy."""
    # Handle z=0 case
    z = np.asarray(z, dtype=complex)
    result = np.ones_like(z, dtype=complex)
    nonzero = np.abs(z) > 1e-10
    result[nonzero] = z[nonzero] / np.tanh(z[nonzero])
    return result


def test_cfkan_dowell():
    """Test if CF-KAN can learn Dowell coefficients."""

    print("=" * 60)
    print("CF-KAN: Learning Dowell Coefficients from z*coth(z)")
    print("=" * 60)

    # Target: Dowell coefficients
    dowell_exact = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
    n_stages = 5

    # Generate training data
    # z = (1+j) * d/delta, typical range for skin effect
    # Use real z for simplicity first (imaginary part similar behavior)
    np.random.seed(42)

    # Training points: z values where skin effect is significant
    z_train_real = np.linspace(0.1, 5.0, 50)
    z_train = z_train_real * (1 + 1j) / np.sqrt(2)  # Complex z

    # Target values
    y_train = zcothz_exact(z_train)

    # Convert to torch tensors (complex)
    z_torch = torch.tensor(z_train, dtype=torch.complex128)
    y_torch = torch.tensor(y_train, dtype=torch.complex128)

    # Create model
    model = CFKAN(n_stages=n_stages, init_mode='random')

    print(f"\nInitial coefficients: {model.get_coefficients()}")
    print(f"Target (Dowell):      {dowell_exact}")

    # Optimizer
    optimizer = optim.Adam(model.parameters(), lr=0.1)
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=500, gamma=0.5)

    # Training loop
    n_epochs = 2000
    losses = []

    print(f"\nTraining for {n_epochs} epochs...")

    for epoch in range(n_epochs):
        optimizer.zero_grad()

        # Forward pass
        y_pred = model(z_torch)

        # Loss: MSE on complex values (real + imag)
        diff = y_pred - y_torch
        loss = torch.mean(torch.abs(diff) ** 2)

        # Backward pass
        loss.backward()

        # Gradient clipping (important for stability)
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)

        optimizer.step()
        scheduler.step()

        # Ensure coefficients stay positive (physical constraint)
        with torch.no_grad():
            model.a.data = torch.clamp(model.a.data, min=0.1)

        losses.append(loss.item())

        if (epoch + 1) % 200 == 0:
            coeffs = model.get_coefficients()
            print(f"  Epoch {epoch+1:4d}: loss={loss.item():.6f}, "
                  f"coeffs={np.array2string(coeffs, precision=2)}")

    # Final results
    final_coeffs = model.get_coefficients()

    print("\n" + "=" * 60)
    print("Results")
    print("=" * 60)
    print(f"Learned coefficients: {final_coeffs}")
    print(f"Target (Dowell):      {dowell_exact}")
    print(f"Absolute errors:      {np.abs(final_coeffs - dowell_exact)}")
    print(f"Relative errors (%):  {100 * np.abs(final_coeffs - dowell_exact) / dowell_exact}")

    # Verify on test points
    z_test = np.linspace(0.5, 4.0, 10) * (1 + 1j) / np.sqrt(2)
    y_test_exact = zcothz_exact(z_test)

    with torch.no_grad():
        y_test_pred = model(torch.tensor(z_test, dtype=torch.complex128)).numpy()

    test_error = np.max(np.abs(y_test_pred - y_test_exact) / np.abs(y_test_exact)) * 100

    print(f"\nTest set max relative error: {test_error:.2f}%")

    # Check if coefficients converged to Dowell values
    coeff_error = np.max(np.abs(final_coeffs - dowell_exact) / dowell_exact) * 100

    if coeff_error < 5.0:
        print(f"\nSUCCESS: Coefficients within 5% of Dowell values!")
        success = True
    else:
        print(f"\nPartial success: Coefficients have {coeff_error:.1f}% max error")
        success = False

    return {
        'success': success,
        'learned_coeffs': final_coeffs,
        'dowell_coeffs': dowell_exact,
        'coeff_error_pct': coeff_error,
        'test_error_pct': test_error,
        'final_loss': losses[-1],
    }


def test_cfkan_convergence_study():
    """Study how initialization affects convergence."""

    print("\n" + "=" * 60)
    print("Convergence Study: Effect of Initialization")
    print("=" * 60)

    init_modes = ['random', 'ones', 'near_dowell']
    results = {}

    for mode in init_modes:
        print(f"\n--- Init mode: {mode} ---")

        torch.manual_seed(123)  # Reproducibility

        model = CFKAN(n_stages=5, init_mode=mode)

        # Quick training
        z_train = np.linspace(0.1, 5.0, 50) * (1 + 1j) / np.sqrt(2)
        y_train = zcothz_exact(z_train)
        z_torch = torch.tensor(z_train, dtype=torch.complex128)
        y_torch = torch.tensor(y_train, dtype=torch.complex128)

        optimizer = optim.Adam(model.parameters(), lr=0.1)

        for epoch in range(1000):
            optimizer.zero_grad()
            y_pred = model(z_torch)
            loss = torch.mean(torch.abs(y_pred - y_torch) ** 2)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            with torch.no_grad():
                model.a.data = torch.clamp(model.a.data, min=0.1)

        final_coeffs = model.get_coefficients()
        dowell = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
        error = np.max(np.abs(final_coeffs - dowell) / dowell) * 100

        print(f"  Final coeffs: {np.array2string(final_coeffs, precision=2)}")
        print(f"  Max error: {error:.1f}%")

        results[mode] = {
            'coeffs': final_coeffs,
            'error_pct': error,
        }

    return results


if __name__ == '__main__':
    print("CF-KAN (Continued Fraction KAN) Test Suite")
    print("=" * 60)

    # Main test
    result = test_cfkan_dowell()

    # Convergence study
    conv_results = test_cfkan_convergence_study()

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    if result['success']:
        print("CF-KAN successfully learned Dowell coefficients [3,5,7,9,11]!")
        print("This proves the concept: CF structure enables coefficient extraction.")
    else:
        print("CF-KAN partially converged. Further tuning may be needed.")

    print(f"\nBest coefficient error: {result['coeff_error_pct']:.2f}%")
    print(f"Function approximation error: {result['test_error_pct']:.2f}%")
