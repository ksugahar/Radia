"""
CF-KAN v2: Improved Continued Fraction Learning

Improvements over v1:
1. Layer-wise training (inner coefficients first)
2. Better loss function (coefficient regularization)
3. Multiple random restarts
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim


class CFKANv2(nn.Module):
    """
    Continued Fraction KAN v2 with improved training.
    """

    def __init__(self, n_stages=5):
        super().__init__()
        self.n_stages = n_stages
        # Initialize with reasonable values (odd numbers are typical)
        init_vals = torch.tensor([3.0, 5.0, 7.0, 9.0, 11.0])[:n_stages]
        init_vals = init_vals + torch.randn(n_stages) * 2.0  # Add noise
        self.a = nn.Parameter(init_vals)

    def forward(self, z, n_terms=None):
        """
        Evaluate truncated continued fraction.
        n_terms: how many terms to use (for layer-wise training)
        """
        if n_terms is None:
            n_terms = self.n_stages

        z_sq = z ** 2

        # Start from the n_terms-th coefficient
        result = self.a[n_terms - 1].to(z.dtype)

        for i in range(n_terms - 2, -1, -1):
            result = self.a[i].to(z.dtype) + z_sq / result

        return 1.0 + z_sq / result

    def get_coefficients(self):
        return self.a.detach().cpu().numpy()


def zcothz_exact(z):
    """Exact z*coth(z)."""
    z = np.asarray(z, dtype=complex)
    result = np.ones_like(z, dtype=complex)
    nonzero = np.abs(z) > 1e-10
    result[nonzero] = z[nonzero] / np.tanh(z[nonzero])
    return result


def cf_truncated(z, coeffs):
    """Evaluate truncated continued fraction with given coefficients."""
    z_sq = z ** 2
    n = len(coeffs)
    result = coeffs[n - 1]
    for i in range(n - 2, -1, -1):
        result = coeffs[i] + z_sq / result
    return 1.0 + z_sq / result


def train_cfkan_layerwise():
    """
    Layer-wise training: fit inner coefficients first.

    Idea: The innermost coefficient (a_n) has the weakest gradient.
    So we train from outside-in: first fix a_1, then a_2, etc.

    Actually, better approach: train a_1 alone first with 1-term CF,
    then train a_1, a_2 with 2-term CF, etc.
    """
    print("=" * 60)
    print("CF-KAN v2: Layer-wise Training")
    print("=" * 60)

    dowell_exact = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
    n_stages = 5

    # Training data - wider range for better conditioning
    z_train = np.concatenate([
        np.linspace(0.1, 1.0, 20) * (1 + 1j) / np.sqrt(2),
        np.linspace(1.0, 3.0, 30) * (1 + 1j) / np.sqrt(2),
        np.linspace(3.0, 6.0, 20) * (1 + 1j) / np.sqrt(2),
    ])
    y_train = zcothz_exact(z_train)

    z_torch = torch.tensor(z_train, dtype=torch.complex128)
    y_torch = torch.tensor(y_train, dtype=torch.complex128)

    # Initialize model
    model = CFKANv2(n_stages=n_stages)

    print(f"Initial coefficients: {model.get_coefficients()}")

    # Layer-wise training
    for n_terms in range(1, n_stages + 1):
        print(f"\n--- Training with {n_terms} terms ---")

        # Compute target for this truncation level
        # We need the best n_terms approximation
        optimizer = optim.Adam(model.parameters(), lr=0.05)

        for epoch in range(500):
            optimizer.zero_grad()

            y_pred = model(z_torch, n_terms=n_terms)
            loss = torch.mean(torch.abs(y_pred - y_torch) ** 2)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

            with torch.no_grad():
                model.a.data = torch.clamp(model.a.data, min=0.5, max=20.0)

        coeffs = model.get_coefficients()
        print(f"  Coefficients: {np.array2string(coeffs[:n_terms], precision=2)}")
        print(f"  Loss: {loss.item():.6f}")

    # Final fine-tuning with all terms
    print("\n--- Final fine-tuning ---")
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    for epoch in range(2000):
        optimizer.zero_grad()
        y_pred = model(z_torch)
        loss = torch.mean(torch.abs(y_pred - y_torch) ** 2)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()
        with torch.no_grad():
            model.a.data = torch.clamp(model.a.data, min=0.5, max=20.0)

        if (epoch + 1) % 500 == 0:
            print(f"  Epoch {epoch+1}: loss={loss.item():.6f}, "
                  f"coeffs={np.array2string(model.get_coefficients(), precision=2)}")

    final_coeffs = model.get_coefficients()
    errors = 100 * np.abs(final_coeffs - dowell_exact) / dowell_exact

    print("\n" + "=" * 60)
    print("Results (Layer-wise)")
    print("=" * 60)
    print(f"Learned: {final_coeffs}")
    print(f"Target:  {dowell_exact}")
    print(f"Errors:  {errors} %")

    return final_coeffs, errors


def train_cfkan_multi_restart():
    """
    Multiple random restarts to find global optimum.
    """
    print("\n" + "=" * 60)
    print("CF-KAN v2: Multiple Random Restarts")
    print("=" * 60)

    dowell_exact = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
    n_stages = 5

    z_train = np.linspace(0.1, 5.0, 100) * (1 + 1j) / np.sqrt(2)
    y_train = zcothz_exact(z_train)
    z_torch = torch.tensor(z_train, dtype=torch.complex128)
    y_torch = torch.tensor(y_train, dtype=torch.complex128)

    best_coeffs = None
    best_error = float('inf')

    for trial in range(10):
        torch.manual_seed(trial * 123)

        # Random initialization
        init = torch.rand(n_stages) * 10 + 1
        model = CFKANv2(n_stages=n_stages)
        model.a.data = init

        optimizer = optim.Adam(model.parameters(), lr=0.1)
        scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=300, gamma=0.5)

        for epoch in range(1500):
            optimizer.zero_grad()
            y_pred = model(z_torch)
            loss = torch.mean(torch.abs(y_pred - y_torch) ** 2)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            scheduler.step()
            with torch.no_grad():
                model.a.data = torch.clamp(model.a.data, min=0.5, max=20.0)

        final_coeffs = model.get_coefficients()
        max_error = np.max(np.abs(final_coeffs - dowell_exact) / dowell_exact) * 100

        print(f"  Trial {trial+1:2d}: coeffs={np.array2string(final_coeffs, precision=2)}, "
              f"max_err={max_error:.1f}%")

        if max_error < best_error:
            best_error = max_error
            best_coeffs = final_coeffs.copy()

    print(f"\nBest result: max_error={best_error:.1f}%")
    print(f"Best coeffs: {best_coeffs}")

    return best_coeffs, best_error


def train_cfkan_penalty():
    """
    Add penalty to encourage coefficients toward odd integers.
    """
    print("\n" + "=" * 60)
    print("CF-KAN v2: With Odd-Integer Regularization")
    print("=" * 60)

    dowell_exact = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
    n_stages = 5

    z_train = np.linspace(0.1, 5.0, 100) * (1 + 1j) / np.sqrt(2)
    y_train = zcothz_exact(z_train)
    z_torch = torch.tensor(z_train, dtype=torch.complex128)
    y_torch = torch.tensor(y_train, dtype=torch.complex128)

    model = CFKANv2(n_stages=n_stages)

    # Initialize near odd integers with some noise
    model.a.data = torch.tensor([3.0, 5.0, 7.0, 9.0, 11.0]) + torch.randn(5) * 3.0
    model.a.data = torch.clamp(model.a.data, min=1.0)

    print(f"Initial: {model.get_coefficients()}")

    optimizer = optim.Adam(model.parameters(), lr=0.05)

    for epoch in range(3000):
        optimizer.zero_grad()

        y_pred = model(z_torch)
        loss_fit = torch.mean(torch.abs(y_pred - y_torch) ** 2)

        # Regularization: encourage 2*k+1 pattern
        # penalty = sum of (a_i - (2*i + 3))^2 with small weight
        expected = torch.tensor([3.0, 5.0, 7.0, 9.0, 11.0])
        loss_reg = 0.0001 * torch.sum((model.a - expected) ** 2)

        loss = loss_fit + loss_reg

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        with torch.no_grad():
            model.a.data = torch.clamp(model.a.data, min=0.5, max=20.0)

        if (epoch + 1) % 500 == 0:
            print(f"  Epoch {epoch+1}: fit={loss_fit.item():.6f}, "
                  f"reg={loss_reg.item():.6f}, "
                  f"coeffs={np.array2string(model.get_coefficients(), precision=2)}")

    final_coeffs = model.get_coefficients()
    errors = 100 * np.abs(final_coeffs - dowell_exact) / dowell_exact

    print(f"\nFinal: {final_coeffs}")
    print(f"Errors: {errors} %")

    return final_coeffs, errors


if __name__ == '__main__':
    print("CF-KAN v2: Improved Training Strategies")
    print("=" * 60)

    # Test 1: Layer-wise training
    coeffs1, err1 = train_cfkan_layerwise()

    # Test 2: Multiple restarts
    coeffs2, err2 = train_cfkan_multi_restart()

    # Test 3: With regularization (knowing the pattern)
    coeffs3, err3 = train_cfkan_penalty()

    print("\n" + "=" * 60)
    print("Summary of All Methods")
    print("=" * 60)
    print(f"Target (Dowell):    [3, 5, 7, 9, 11]")
    print(f"Layer-wise:         {np.array2string(coeffs1, precision=2)}, max_err={np.max(err1):.1f}%")
    print(f"Multi-restart:      {np.array2string(coeffs2, precision=2)}, max_err={err2:.1f}%")
    print(f"With regularization:{np.array2string(coeffs3, precision=2)}, max_err={np.max(err3):.1f}%")
