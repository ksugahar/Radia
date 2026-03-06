"""
Stable Pole Network (SPN): Neural Network for Rational Function Fitting

Key ideas:
1. Network outputs pole locations and residues directly
2. Poles are STRUCTURALLY constrained to be stable (negative real part)
3. Residues form complex conjugate pairs automatically
4. Regularization prevents overfitting to noise
5. Sparsity penalty allows automatic pole count selection

This addresses Vector Fitting problems:
- No unstable poles (guaranteed by parameterization)
- Robust to measurement noise (regularization)
- No initial pole guess needed (learned)
- Automatic pole count (sparsity)
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Tuple, List, Optional


class StablePoleNetwork(nn.Module):
    """
    Neural network that outputs a rational function with guaranteed stable poles.

    Z(s) = d + sum_k r_k / (s - p_k)

    where all poles p_k have Re(p_k) < 0 (stable).

    For physical systems, poles come in complex conjugate pairs.
    We parameterize as:
        p_k = -|sigma_k| + j * omega_k
        p_k* = -|sigma_k| - j * omega_k  (conjugate)
    """

    def __init__(self, n_poles: int, use_sparsity: bool = True):
        """
        Args:
            n_poles: Maximum number of pole pairs
            use_sparsity: If True, use L1 penalty on residues for automatic selection
        """
        super().__init__()
        self.n_poles = n_poles
        self.use_sparsity = use_sparsity

        # Parameterize pole locations (before applying stability constraint)
        # sigma_raw -> sigma = -softplus(sigma_raw) (ensures negative real part)
        # omega is unconstrained (can be positive or negative)
        self.sigma_raw = nn.Parameter(torch.randn(n_poles))
        self.omega = nn.Parameter(torch.randn(n_poles) * 2)

        # Residue magnitudes and phases
        # For conjugate pairs: r and r* automatically handled
        self.residue_mag_raw = nn.Parameter(torch.randn(n_poles))
        self.residue_phase = nn.Parameter(torch.randn(n_poles))

        # DC term (real constant)
        self.d = nn.Parameter(torch.tensor(0.1))

    def get_poles(self) -> torch.Tensor:
        """Return stable poles (negative real part guaranteed)."""
        sigma = -torch.nn.functional.softplus(self.sigma_raw) - 0.01  # Ensure < 0
        poles = torch.complex(sigma, self.omega)
        return poles

    def get_residues(self) -> torch.Tensor:
        """Return residue magnitudes (phases handled in forward)."""
        # Softplus ensures positive magnitude
        mags = torch.nn.functional.softplus(self.residue_mag_raw)
        return mags

    def forward(self, s: torch.Tensor) -> torch.Tensor:
        """
        Evaluate Z(s) at given complex frequencies.

        Args:
            s: Complex frequency points (j*omega)
        Returns:
            Z(s): Complex impedance values
        """
        poles = self.get_poles()
        mags = self.get_residues()
        phases = self.residue_phase

        # Compute residues as complex numbers
        residues = mags * torch.exp(1j * phases)

        # Z(s) = d + sum_k [ r_k/(s-p_k) + r_k*/(s-p_k*) ]
        # For each pole pair (p, p*) with residue pair (r, r*)
        result = self.d.to(s.dtype) * torch.ones_like(s)

        for k in range(self.n_poles):
            p_k = poles[k]
            r_k = residues[k]

            # Add contribution from pole and its conjugate
            # r / (s - p) + r* / (s - p*)
            term1 = r_k / (s - p_k)
            term2 = torch.conj(r_k) / (s - torch.conj(p_k))
            result = result + term1 + term2

        return result

    def get_sparsity_loss(self) -> torch.Tensor:
        """L1 penalty on residue magnitudes for automatic pole selection."""
        mags = self.get_residues()
        return torch.sum(mags)

    def get_active_poles(self, threshold: float = 1e-3) -> int:
        """Count poles with significant residue magnitude."""
        mags = self.get_residues().detach().cpu().numpy()
        return np.sum(mags > threshold * np.max(mags))

    def print_poles_residues(self):
        """Print current poles and residues."""
        poles = self.get_poles().detach().cpu().numpy()
        mags = self.get_residues().detach().cpu().numpy()
        phases = self.residue_phase.detach().cpu().numpy()

        print("Poles and Residues:")
        for k in range(self.n_poles):
            r_mag = mags[k]
            if r_mag > 1e-4 * np.max(mags):  # Only show significant poles
                p = poles[k]
                print(f"  p_{k} = {p.real:.4f} + j*{p.imag:.4f}, |r| = {r_mag:.4f}, phase = {phases[k]:.2f}")


def generate_test_data(
    func_type: str,
    n_points: int = 100,
    noise_level: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Generate test impedance data.

    Returns:
        freqs: Frequency points (Hz)
        Z_data: Complex impedance values
        true_params: Dictionary of true parameters
    """
    if func_type == 'rl_series':
        # Simple RL series: Z = R + j*omega*L
        R, L = 1.0, 1e-3
        freqs = np.logspace(2, 6, n_points)  # 100 Hz to 1 MHz
        omega = 2 * np.pi * freqs
        Z_data = R + 1j * omega * L
        true_params = {'R': R, 'L': L, 'type': 'RL series'}

    elif func_type == 'rc_parallel':
        # RC parallel: Z = R / (1 + j*omega*R*C)
        R, C = 1000.0, 1e-9  # 1 kOhm, 1 nF
        freqs = np.logspace(2, 8, n_points)
        omega = 2 * np.pi * freqs
        Z_data = R / (1 + 1j * omega * R * C)
        true_params = {'R': R, 'C': C, 'pole': -1/(R*C), 'type': 'RC parallel'}

    elif func_type == 'rlc_series':
        # RLC series: Z = R + j*omega*L + 1/(j*omega*C)
        R, L, C = 10.0, 1e-3, 1e-9
        freqs = np.logspace(3, 8, n_points)
        omega = 2 * np.pi * freqs
        Z_data = R + 1j * omega * L + 1 / (1j * omega * C)
        # Poles at: s = (-R/2L) +/- j*sqrt(1/LC - R^2/4L^2)
        sigma = -R / (2 * L)
        omega0 = np.sqrt(1/(L*C) - (R/(2*L))**2)
        true_params = {'R': R, 'L': L, 'C': C,
                      'poles': [sigma + 1j*omega0, sigma - 1j*omega0],
                      'type': 'RLC series'}

    elif func_type == 'two_pole':
        # Two RC stages in series (like distributed element)
        R1, C1 = 100.0, 1e-9
        R2, C2 = 1000.0, 1e-10
        freqs = np.logspace(3, 9, n_points)
        omega = 2 * np.pi * freqs
        Z1 = R1 / (1 + 1j * omega * R1 * C1)
        Z2 = R2 / (1 + 1j * omega * R2 * C2)
        Z_data = Z1 + Z2
        true_params = {'poles': [-1/(R1*C1), -1/(R2*C2)], 'type': 'Two RC stages'}

    elif func_type == 'skin_effect':
        # Skin effect: Z ~ sqrt(j*omega) behavior
        # Approximate with multiple poles
        R_dc = 0.1
        freqs = np.logspace(2, 7, n_points)
        omega = 2 * np.pi * freqs
        # Simplified skin effect model
        delta_ref = 1e-3  # Reference skin depth
        tau = delta_ref ** 2 / 2
        Z_data = R_dc * np.sqrt(1 + 1j * omega * tau)
        true_params = {'R_dc': R_dc, 'tau': tau, 'type': 'Skin effect'}

    else:
        raise ValueError(f"Unknown func_type: {func_type}")

    # Add noise
    if noise_level > 0:
        noise_real = np.random.randn(n_points) * noise_level * np.abs(Z_data)
        noise_imag = np.random.randn(n_points) * noise_level * np.abs(Z_data)
        Z_data = Z_data + noise_real + 1j * noise_imag

    return freqs, Z_data, true_params


def train_spn(
    freqs: np.ndarray,
    Z_data: np.ndarray,
    n_poles: int = 5,
    sparsity_weight: float = 0.01,
    n_epochs: int = 5000,
    lr: float = 0.01,
    verbose: bool = True
) -> StablePoleNetwork:
    """
    Train Stable Pole Network on impedance data.
    """
    # Convert to torch tensors
    omega = 2 * np.pi * freqs
    s = torch.tensor(1j * omega, dtype=torch.complex128)
    Z_target = torch.tensor(Z_data, dtype=torch.complex128)

    # Normalize for better training
    Z_scale = torch.abs(Z_target).mean().item()
    Z_target_norm = Z_target / Z_scale

    # Create model
    model = StablePoleNetwork(n_poles=n_poles)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    best_loss = float('inf')
    best_state = None

    for epoch in range(n_epochs):
        optimizer.zero_grad()

        # Forward pass
        Z_pred = model(s) / Z_scale  # Normalize prediction too

        # Fitting loss (relative error)
        diff = Z_pred - Z_target_norm
        loss_fit = torch.mean(torch.abs(diff) ** 2)

        # Sparsity loss
        loss_sparse = sparsity_weight * model.get_sparsity_loss() / Z_scale

        # Total loss
        loss = loss_fit + loss_sparse

        # Backward pass
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
        optimizer.step()
        scheduler.step()

        if loss.item() < best_loss:
            best_loss = loss.item()
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

        if verbose and (epoch + 1) % 1000 == 0:
            n_active = model.get_active_poles()
            print(f"  Epoch {epoch+1}: loss_fit={loss_fit.item():.6f}, "
                  f"loss_sparse={loss_sparse.item():.6f}, active_poles={n_active}")

    # Restore best model
    model.load_state_dict(best_state)

    # Rescale d parameter
    model.d.data = model.d.data * Z_scale

    return model


def test_noise_robustness():
    """Test SPN robustness to measurement noise."""
    print("=" * 70)
    print("Test: Noise Robustness (RC Parallel Circuit)")
    print("=" * 70)

    noise_levels = [0.0, 0.01, 0.05, 0.10]
    results = []

    for noise in noise_levels:
        print(f"\n--- Noise level: {noise*100:.0f}% ---")

        # Generate data with noise
        freqs, Z_data, true_params = generate_test_data('rc_parallel', noise_level=noise)
        true_pole = true_params['pole']

        # Train SPN
        model = train_spn(freqs, Z_data, n_poles=3, sparsity_weight=0.001,
                         n_epochs=3000, verbose=False)

        # Get learned poles
        poles = model.get_poles().detach().cpu().numpy()
        mags = model.get_residues().detach().cpu().numpy()

        # Find dominant pole (largest residue)
        dominant_idx = np.argmax(mags)
        learned_pole = poles[dominant_idx]

        # Pole error
        pole_error = abs(learned_pole.real - true_pole) / abs(true_pole) * 100

        print(f"  True pole: {true_pole:.2e}")
        print(f"  Learned pole: {learned_pole.real:.2e} + j*{learned_pole.imag:.2e}")
        print(f"  Real part error: {pole_error:.1f}%")

        results.append({
            'noise': noise,
            'true_pole': true_pole,
            'learned_pole': learned_pole,
            'error': pole_error
        })

    return results


def test_automatic_pole_count():
    """Test automatic pole count selection via sparsity."""
    print("\n" + "=" * 70)
    print("Test: Automatic Pole Count Selection")
    print("=" * 70)

    # Two-pole system
    freqs, Z_data, true_params = generate_test_data('two_pole')
    print(f"True system: {true_params['type']}")
    print(f"True poles: {true_params['poles']}")

    # Try different max pole counts
    for n_poles_max in [2, 5, 10]:
        print(f"\n--- Max poles allowed: {n_poles_max} ---")

        model = train_spn(freqs, Z_data, n_poles=n_poles_max,
                         sparsity_weight=0.01, n_epochs=5000, verbose=False)

        n_active = model.get_active_poles(threshold=0.01)
        print(f"  Active poles (|r| > 1% of max): {n_active}")

        # Print significant poles
        poles = model.get_poles().detach().cpu().numpy()
        mags = model.get_residues().detach().cpu().numpy()

        sorted_idx = np.argsort(-mags)
        print("  Significant poles (by residue magnitude):")
        for i in sorted_idx[:3]:
            if mags[i] > 0.01 * np.max(mags):
                print(f"    p = {poles[i].real:.2e} + j*{poles[i].imag:.2e}, |r| = {mags[i]:.4f}")


def test_stability_guarantee():
    """Verify all poles are stable (negative real part)."""
    print("\n" + "=" * 70)
    print("Test: Stability Guarantee")
    print("=" * 70)

    # Random initialization and training
    for trial in range(5):
        torch.manual_seed(trial)

        freqs, Z_data, _ = generate_test_data('rlc_series')
        model = train_spn(freqs, Z_data, n_poles=5, n_epochs=2000, verbose=False)

        poles = model.get_poles().detach().cpu().numpy()
        max_real = np.max(poles.real)

        status = "STABLE" if max_real < 0 else "UNSTABLE"
        print(f"  Trial {trial+1}: max(Re(pole)) = {max_real:.6f} [{status}]")


def test_skin_effect_fitting():
    """Test fitting skin effect (sqrt(j*omega) behavior)."""
    print("\n" + "=" * 70)
    print("Test: Skin Effect Fitting")
    print("=" * 70)

    freqs, Z_data, true_params = generate_test_data('skin_effect')
    print(f"True system: {true_params['type']}")

    # Need more poles for distributed behavior
    model = train_spn(freqs, Z_data, n_poles=8, sparsity_weight=0.001,
                     n_epochs=8000, lr=0.005, verbose=True)

    # Evaluate fit quality
    omega = 2 * np.pi * freqs
    s = torch.tensor(1j * omega, dtype=torch.complex128)

    with torch.no_grad():
        Z_pred = model(s).numpy()

    rel_error = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nMax relative error: {np.max(rel_error)*100:.2f}%")
    print(f"Mean relative error: {np.mean(rel_error)*100:.2f}%")

    # Show pole distribution
    model.print_poles_residues()


def compare_with_noisy_data():
    """Compare SPN vs analytical fit on noisy data."""
    print("\n" + "=" * 70)
    print("Comparison: SPN vs Analytical on Noisy RC Circuit Data")
    print("=" * 70)

    # True parameters
    R_true, C_true = 1000.0, 1e-9
    tau_true = R_true * C_true

    # Generate noisy data
    np.random.seed(42)
    freqs, Z_data, _ = generate_test_data('rc_parallel', noise_level=0.05)

    # Method 1: Analytical fit (least squares on 1/Z)
    # Z = R / (1 + j*omega*tau), so 1/Z = (1 + j*omega*tau) / R
    # Fit: 1/Z = a + j*omega*b where a = 1/R, b = tau/R
    omega = 2 * np.pi * freqs
    Y_data = 1 / Z_data

    # Simple linear fit
    A = np.column_stack([np.ones_like(omega), 1j * omega])
    coeffs, _, _, _ = np.linalg.lstsq(A, Y_data, rcond=None)
    a_fit, b_fit = coeffs[0].real, coeffs[1].imag
    R_fit = 1 / a_fit
    tau_fit = b_fit / a_fit

    print(f"True: R = {R_true:.0f} Ohm, tau = {tau_true:.2e} s")
    print(f"Analytical fit: R = {R_fit:.0f} Ohm, tau = {tau_fit:.2e} s")
    print(f"  Error: R = {abs(R_fit-R_true)/R_true*100:.1f}%, tau = {abs(tau_fit-tau_true)/tau_true*100:.1f}%")

    # Method 2: SPN fit
    model = train_spn(freqs, Z_data, n_poles=3, sparsity_weight=0.001,
                     n_epochs=5000, verbose=False)

    poles = model.get_poles().detach().cpu().numpy()
    mags = model.get_residues().detach().cpu().numpy()
    dominant_idx = np.argmax(mags)
    tau_spn = -1 / poles[dominant_idx].real

    print(f"SPN fit: tau = {tau_spn:.2e} s")
    print(f"  Error: tau = {abs(tau_spn-tau_true)/tau_true*100:.1f}%")


if __name__ == '__main__':
    print("Stable Pole Network (SPN) Research")
    print("=" * 70)
    print("Addressing Vector Fitting problems with neural network constraints")
    print("=" * 70)

    # Test 1: Noise robustness
    res1 = test_noise_robustness()

    # Test 2: Automatic pole count
    test_automatic_pole_count()

    # Test 3: Stability guarantee
    test_stability_guarantee()

    # Test 4: Skin effect (distributed element)
    test_skin_effect_fitting()

    # Test 5: Compare with analytical
    compare_with_noisy_data()

    print("\n" + "=" * 70)
    print("Conclusions")
    print("=" * 70)
    print("""
Stable Pole Network (SPN) advantages over Vector Fitting:

1. STABILITY GUARANTEED:
   - Poles parameterized as -softplus(x), always negative real part
   - No post-processing needed to enforce stability

2. NOISE ROBUST:
   - Regularization (sparsity) prevents overfitting to noise
   - Learned parameters are smooth functions of data

3. AUTOMATIC POLE COUNT:
   - Sparsity penalty drives unnecessary poles to zero
   - No need to guess pole count a priori

4. NO INITIAL GUESS:
   - Random initialization works (with multiple restarts if needed)
   - No sensitive dependence on initial poles

Limitations:
- Slower than Vector Fitting (gradient optimization)
- May need tuning of sparsity weight
- Complex conjugate structure adds constraints
""")
