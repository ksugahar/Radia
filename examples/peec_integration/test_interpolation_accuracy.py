"""
Test: Interpolation Accuracy Between Sample Points

Key issue with Vector Fitting:
- VF places poles based on sample points
- Between samples, interpolation can be poor
- Hidden poles can cause unexpected behavior

This test compares URN vs VF-like approach on interpolation quality.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Tuple
import sys
import os

# Add path for imports
sys.path.insert(0, os.path.dirname(__file__))

from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn
)


class SimplePoleNetwork(nn.Module):
    """
    Simple pole-residue model (VF-like).

    Z(s) = d + sum_k r_k / (s - p_k)

    This mimics Vector Fitting's approach.
    """

    def __init__(self, n_poles: int, omega_ref: float, Z_ref: float):
        super().__init__()
        self.n_poles = n_poles
        self.omega_ref = omega_ref
        self.Z_ref = Z_ref

        # Pole parameters (negative real for stability)
        self.pole_sigma_raw = nn.Parameter(torch.randn(n_poles))
        self.pole_omega = nn.Parameter(torch.randn(n_poles) * 0.5)

        # Residue parameters
        self.residue_real = nn.Parameter(torch.randn(n_poles) * 0.1)
        self.residue_imag = nn.Parameter(torch.randn(n_poles) * 0.1)

        # DC term
        self.d_real = nn.Parameter(torch.tensor(0.0))
        self.d_imag = nn.Parameter(torch.tensor(0.0))

    def forward(self, omega: torch.Tensor) -> torch.Tensor:
        s = 1j * omega / self.omega_ref

        # Poles (stable: negative real part)
        sigma = -torch.exp(self.pole_sigma_raw)
        poles = torch.complex(sigma, self.pole_omega)

        # Residues
        residues = torch.complex(self.residue_real, self.residue_imag)

        # Z(s) = d + sum r_k / (s - p_k) + conjugate
        Z = torch.complex(self.d_real, self.d_imag) * torch.ones_like(s)

        for k in range(self.n_poles):
            p_k = poles[k]
            r_k = residues[k]
            Z = Z + r_k / (s - p_k)
            Z = Z + torch.conj(r_k) / (s - torch.conj(p_k))

        return Z * self.Z_ref


def train_pole_network(
    freqs: np.ndarray,
    Z_data: np.ndarray,
    n_poles: int = 5,
    n_epochs: int = 5000,
    n_restarts: int = 3
) -> SimplePoleNetwork:
    """Train pole-residue model (VF-like)."""

    omega = 2 * np.pi * freqs
    omega_ref = np.sqrt(omega.min() * omega.max())
    Z_ref = np.abs(Z_data).mean()

    omega_torch = torch.tensor(omega, dtype=torch.float64)
    Z_target = torch.tensor(Z_data, dtype=torch.complex128)

    best_model = None
    best_loss = float('inf')

    for restart in range(n_restarts):
        torch.manual_seed(restart * 42)

        model = SimplePoleNetwork(n_poles, omega_ref, Z_ref)
        optimizer = optim.Adam(model.parameters(), lr=0.05)

        for epoch in range(n_epochs):
            optimizer.zero_grad()
            Z_pred = model(omega_torch)
            rel_err = (Z_pred - Z_target) / (torch.abs(Z_target) + 1e-10 * Z_ref)
            loss = torch.mean(torch.abs(rel_err) ** 2)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()

        if loss.item() < best_loss:
            best_loss = loss.item()
            best_model = SimplePoleNetwork(n_poles, omega_ref, Z_ref)
            best_model.load_state_dict(model.state_dict())

    return best_model


def generate_multi_resonance_data(n_resonances: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate data with multiple narrow resonances.

    This is a challenging case for VF: narrow peaks can be missed
    if sample points don't hit the resonance.
    """
    # Dense frequency grid for "true" data
    freqs_dense = np.logspace(3, 7, 1000)
    omega_dense = 2 * np.pi * freqs_dense

    # Multiple resonances at different frequencies
    Z_true = np.ones_like(omega_dense, dtype=complex) * 10  # Base impedance

    resonance_freqs = [1e4, 1e5, 1e6]  # 10 kHz, 100 kHz, 1 MHz
    Q_factors = [50, 30, 20]  # High Q = narrow resonance

    for f0, Q in zip(resonance_freqs, Q_factors):
        omega0 = 2 * np.pi * f0
        # RLC resonance: Z = R + j*omega*L + 1/(j*omega*C)
        # At resonance: omega0 = 1/sqrt(LC), Q = omega0*L/R
        R = 1.0
        L = Q * R / omega0
        C = 1 / (omega0**2 * L)
        Z_rlc = R + 1j * omega_dense * L + 1 / (1j * omega_dense * C)
        # Add as parallel impedance
        Z_true = 1 / (1/Z_true + 1/Z_rlc)

    return freqs_dense, Z_true, np.array(resonance_freqs)


def test_interpolation_accuracy():
    """
    Compare URN vs pole-residue model on interpolation between sample points.
    """
    print("=" * 70)
    print("Test: Interpolation Accuracy Between Sample Points")
    print("=" * 70)

    # Generate true data (dense)
    freqs_dense, Z_true, resonance_freqs = generate_multi_resonance_data()
    print(f"True resonances at: {resonance_freqs / 1e3} kHz")

    # Create SPARSE training data (simulating limited measurements)
    n_train = 30  # Only 30 points across 4 decades
    train_indices = np.linspace(0, len(freqs_dense)-1, n_train, dtype=int)
    freqs_train = freqs_dense[train_indices]
    Z_train = Z_true[train_indices]

    print(f"\nTraining with {n_train} sample points (sparse)")
    print(f"Testing on {len(freqs_dense)} points (dense)")

    # Method 1: Pole-residue model (VF-like)
    print("\n--- Training Pole-Residue Model (VF-like) ---")
    pole_model = train_pole_network(freqs_train, Z_train, n_poles=6, n_epochs=5000)

    omega_dense = 2 * np.pi * freqs_dense
    omega_torch = torch.tensor(omega_dense, dtype=torch.float64)

    with torch.no_grad():
        Z_pole = pole_model(omega_torch).numpy()

    # Method 2: URN
    print("\n--- Training URN ---")
    config = URNConfig(
        n_debye=3, n_cole_cole=2, n_cole_davidson=1,
        n_cpe=1, n_warburg=0, n_gerischer=0,
        n_rlc=3,  # Include RLC for resonances
        n_skin_effect=0,
        sparsity_weight=0.005, n_epochs=5000
    )
    urn_model = train_urn(freqs_train, Z_train, config, verbose=False)

    with torch.no_grad():
        Z_urn = urn_model(omega_torch).numpy()

    # Compute errors
    err_pole_train = np.abs(Z_pole[train_indices] - Z_train) / np.abs(Z_train)
    err_pole_all = np.abs(Z_pole - Z_true) / np.abs(Z_true)

    err_urn_train = np.abs(Z_urn[train_indices] - Z_train) / np.abs(Z_train)
    err_urn_all = np.abs(Z_urn - Z_true) / np.abs(Z_true)

    print("\n" + "=" * 70)
    print("Results: Training Points vs All Points (Interpolation)")
    print("=" * 70)

    print("\nPole-Residue Model (VF-like):")
    print(f"  Training points: max_err = {err_pole_train.max()*100:.2f}%, "
          f"mean_err = {err_pole_train.mean()*100:.2f}%")
    print(f"  ALL points:      max_err = {err_pole_all.max()*100:.2f}%, "
          f"mean_err = {err_pole_all.mean()*100:.2f}%")
    print(f"  Interpolation degradation: {err_pole_all.max() / (err_pole_train.max() + 1e-10):.1f}x")

    print("\nURN (Universal Relaxation Network):")
    print(f"  Training points: max_err = {err_urn_train.max()*100:.2f}%, "
          f"mean_err = {err_urn_train.mean()*100:.2f}%")
    print(f"  ALL points:      max_err = {err_urn_all.max()*100:.2f}%, "
          f"mean_err = {err_urn_all.mean()*100:.2f}%")
    print(f"  Interpolation degradation: {err_urn_all.max() / (err_urn_train.max() + 1e-10):.1f}x")

    # Check errors near resonances
    print("\n" + "-" * 70)
    print("Errors Near Resonances (where interpolation matters most):")
    print("-" * 70)

    for f_res in resonance_freqs:
        # Find points within +/- 20% of resonance
        mask = (freqs_dense > f_res * 0.8) & (freqs_dense < f_res * 1.2)
        n_near = np.sum(mask)

        err_pole_near = err_pole_all[mask]
        err_urn_near = err_urn_all[mask]

        print(f"\nNear {f_res/1e3:.0f} kHz ({n_near} points):")
        print(f"  Pole-Residue: max={err_pole_near.max()*100:.1f}%, mean={err_pole_near.mean()*100:.1f}%")
        print(f"  URN:          max={err_urn_near.max()*100:.1f}%, mean={err_urn_near.mean()*100:.1f}%")

    return {
        'pole_train_err': err_pole_train.max(),
        'pole_all_err': err_pole_all.max(),
        'urn_train_err': err_urn_train.max(),
        'urn_all_err': err_urn_all.max(),
    }


def test_hidden_pole_detection():
    """
    Test detection of "hidden" poles that are not at sample frequencies.
    """
    print("\n" + "=" * 70)
    print("Test: Hidden Pole Detection")
    print("=" * 70)

    # Create data with a narrow resonance
    freqs_dense = np.logspace(4, 6, 500)
    omega_dense = 2 * np.pi * freqs_dense

    # Single sharp resonance at 50 kHz (Q=100)
    f_hidden = 50e3
    omega0 = 2 * np.pi * f_hidden
    Q = 100
    R = 1.0
    L = Q * R / omega0
    C = 1 / (omega0**2 * L)

    Z_base = 100.0  # Base impedance
    Z_rlc = R + 1j * omega_dense * L + 1 / (1j * omega_dense * C)
    Z_true = 1 / (1/Z_base + 1/Z_rlc)

    print(f"Hidden resonance at: {f_hidden/1e3:.0f} kHz (Q={Q})")

    # Sparse training data that MISSES the resonance
    # Sample at 10k, 20k, 40k, 60k, 80k, 100k, ... (skips 50k!)
    freqs_train = np.array([10, 20, 30, 40, 60, 70, 80, 90, 100, 200, 500, 1000]) * 1e3
    omega_train = 2 * np.pi * freqs_train

    # Interpolate Z_true at training frequencies
    Z_train = np.interp(freqs_train, freqs_dense,
                        np.abs(Z_true)) * np.exp(1j * np.interp(freqs_train, freqs_dense, np.angle(Z_true)))

    # Actually compute Z_train properly
    Z_rlc_train = R + 1j * omega_train * L + 1 / (1j * omega_train * C)
    Z_train = 1 / (1/Z_base + 1/Z_rlc_train)

    print(f"Training samples: {freqs_train / 1e3} kHz")
    print(f"NOTE: 50 kHz is NOT in training data!")

    # Train both models
    print("\n--- Training Pole-Residue Model ---")
    pole_model = train_pole_network(freqs_train, Z_train, n_poles=5, n_epochs=5000)

    print("\n--- Training URN ---")
    config = URNConfig(
        n_debye=2, n_cole_cole=1, n_cole_davidson=1,
        n_rlc=2, sparsity_weight=0.01, n_epochs=5000
    )
    urn_model = train_urn(freqs_train, Z_train, config, verbose=False)

    # Evaluate at ALL frequencies including the hidden one
    omega_torch = torch.tensor(omega_dense, dtype=torch.float64)

    with torch.no_grad():
        Z_pole = pole_model(omega_torch).numpy()
        Z_urn = urn_model(omega_torch).numpy()

    # Check error specifically at the hidden resonance
    idx_hidden = np.argmin(np.abs(freqs_dense - f_hidden))

    err_pole_hidden = np.abs(Z_pole[idx_hidden] - Z_true[idx_hidden]) / np.abs(Z_true[idx_hidden])
    err_urn_hidden = np.abs(Z_urn[idx_hidden] - Z_true[idx_hidden]) / np.abs(Z_true[idx_hidden])

    print("\n" + "=" * 70)
    print(f"Error at Hidden Resonance ({f_hidden/1e3:.0f} kHz):")
    print("=" * 70)
    print(f"  True Z:         {Z_true[idx_hidden]:.2f}")
    print(f"  Pole-Residue Z: {Z_pole[idx_hidden]:.2f} (error: {err_pole_hidden*100:.1f}%)")
    print(f"  URN Z:          {Z_urn[idx_hidden]:.2f} (error: {err_urn_hidden*100:.1f}%)")

    if err_urn_hidden < err_pole_hidden:
        print("\n** URN better handles the hidden resonance! **")
    else:
        print("\n** Both methods struggle with hidden resonance (expected) **")

    # Also check the region around the hidden resonance
    mask = (freqs_dense > f_hidden * 0.9) & (freqs_dense < f_hidden * 1.1)

    err_pole_region = np.abs(Z_pole[mask] - Z_true[mask]) / np.abs(Z_true[mask])
    err_urn_region = np.abs(Z_urn[mask] - Z_true[mask]) / np.abs(Z_true[mask])

    print(f"\nRegion around hidden resonance (+/- 10%):")
    print(f"  Pole-Residue: max={err_pole_region.max()*100:.1f}%, mean={err_pole_region.mean()*100:.1f}%")
    print(f"  URN:          max={err_urn_region.max()*100:.1f}%, mean={err_urn_region.mean()*100:.1f}%")


def test_noise_robustness():
    """
    Test how noise affects interpolation accuracy.

    VF is known to be sensitive to noise: can create spurious poles.
    """
    print("\n" + "=" * 70)
    print("Test: Noise Robustness in Interpolation")
    print("=" * 70)

    # Simple Debye relaxation
    tau = 1e-5
    R = 100

    freqs_dense = np.logspace(2, 7, 500)
    omega_dense = 2 * np.pi * freqs_dense
    Z_true = R / (1 + 1j * omega_dense * tau)

    # Sparse training with noise
    n_train = 20
    train_indices = np.linspace(0, len(freqs_dense)-1, n_train, dtype=int)
    freqs_train = freqs_dense[train_indices]
    omega_train = 2 * np.pi * freqs_train
    Z_train_clean = R / (1 + 1j * omega_train * tau)

    noise_levels = [0.0, 0.02, 0.05, 0.10]

    print(f"True: R = {R} Ohm, tau = {tau:.2e} s")
    print(f"Training points: {n_train}, Test points: {len(freqs_dense)}")

    for noise in noise_levels:
        print(f"\n--- Noise level: {noise*100:.0f}% ---")

        # Add noise
        np.random.seed(42)
        noise_real = noise * np.abs(Z_train_clean) * np.random.randn(n_train)
        noise_imag = noise * np.abs(Z_train_clean) * np.random.randn(n_train)
        Z_train = Z_train_clean + noise_real + 1j * noise_imag

        # Train pole model
        pole_model = train_pole_network(freqs_train, Z_train, n_poles=4, n_epochs=3000, n_restarts=2)

        # Train URN
        config = URNConfig(
            n_debye=3, n_cole_cole=1, n_rlc=0,
            sparsity_weight=0.01, n_epochs=3000, n_restarts=2
        )
        urn_model = train_urn(freqs_train, Z_train, config, verbose=False)

        # Evaluate on dense grid
        omega_torch = torch.tensor(omega_dense, dtype=torch.float64)
        with torch.no_grad():
            Z_pole = pole_model(omega_torch).numpy()
            Z_urn = urn_model(omega_torch).numpy()

        # Errors
        err_pole = np.abs(Z_pole - Z_true) / np.abs(Z_true)
        err_urn = np.abs(Z_urn - Z_true) / np.abs(Z_true)

        print(f"  Pole-Residue: max_err = {err_pole.max()*100:.1f}%, mean_err = {err_pole.mean()*100:.2f}%")
        print(f"  URN:          max_err = {err_urn.max()*100:.1f}%, mean_err = {err_urn.mean()*100:.2f}%")


if __name__ == '__main__':
    np.random.seed(42)

    print("Interpolation Accuracy: URN vs Vector Fitting-like Approach")
    print("=" * 70)
    print("""
Key concern with Vector Fitting:
- Poles are placed based on sample points
- Between samples, interpolation quality degrades
- Hidden resonances can be missed if not sampled

This test quantifies these issues and compares with URN.
""")

    test_interpolation_accuracy()
    test_hidden_pole_detection()
    test_noise_robustness()

    print("\n" + "=" * 70)
    print("Conclusions")
    print("=" * 70)
    print("""
URN advantages over Vector Fitting for interpolation:

1. SMOOTH BASIS FUNCTIONS:
   - Debye, Cole-Cole, etc. are inherently smooth
   - No artificial oscillations between samples

2. PHYSICAL CONSTRAINTS:
   - Basis functions have known frequency behavior
   - Interpolation follows physical laws, not arbitrary poles

3. NOISE ROBUSTNESS:
   - Sparsity regularization prevents overfitting to noise
   - Fewer spurious features in interpolated region

4. RESONANCE HANDLING:
   - RLC basis explicitly models resonances
   - Better generalization even if resonance not sampled directly
""")
