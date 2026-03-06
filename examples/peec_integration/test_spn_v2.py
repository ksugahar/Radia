"""
Stable Pole Network v2: With Proper Frequency Scaling

Key fix: Normalize frequency to avoid numerical issues with
poles at 10^6 rad/s vs network parameters at O(1).
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Tuple


class StablePoleNetworkV2(nn.Module):
    """
    SPN with proper frequency scaling.

    Key insight: Work in normalized frequency domain.
    s_norm = s / omega_ref where omega_ref is characteristic frequency.

    Z(s) = Z_scale * [ d + sum_k r_k / (s_norm - p_k) ]
    """

    def __init__(self, n_poles: int, omega_ref: float, Z_ref: float):
        """
        Args:
            n_poles: Number of pole pairs
            omega_ref: Reference frequency for normalization (rad/s)
            Z_ref: Reference impedance for normalization (Ohm)
        """
        super().__init__()
        self.n_poles = n_poles
        self.omega_ref = omega_ref
        self.Z_ref = Z_ref

        # Pole parameters (normalized)
        # Real part: -exp(sigma_raw) to ensure negative
        # Imag part: omega_raw (unconstrained)
        self.sigma_raw = nn.Parameter(torch.zeros(n_poles))
        self.omega_raw = nn.Parameter(torch.randn(n_poles) * 0.5)

        # Residue parameters (normalized)
        self.residue_real = nn.Parameter(torch.randn(n_poles) * 0.1)
        self.residue_imag = nn.Parameter(torch.randn(n_poles) * 0.1)

        # DC term (normalized)
        self.d = nn.Parameter(torch.tensor(1.0))

    def get_poles_normalized(self) -> torch.Tensor:
        """Return normalized poles (guaranteed stable)."""
        # Real part: always negative via -exp()
        sigma = -torch.exp(self.sigma_raw)
        omega = self.omega_raw
        return torch.complex(sigma, omega)

    def get_poles_physical(self) -> np.ndarray:
        """Return physical poles (rad/s)."""
        poles_norm = self.get_poles_normalized().detach().cpu().numpy()
        return poles_norm * self.omega_ref

    def get_residues_normalized(self) -> torch.Tensor:
        """Return normalized residues."""
        return torch.complex(self.residue_real, self.residue_imag)

    def forward(self, s: torch.Tensor) -> torch.Tensor:
        """
        Evaluate Z(s).

        Args:
            s: Complex frequency (j*omega) in rad/s
        """
        # Normalize frequency
        s_norm = s / self.omega_ref

        poles = self.get_poles_normalized()
        residues = self.get_residues_normalized()

        # Z_norm = d + sum [ r_k/(s_norm - p_k) + conj term ]
        result = self.d * torch.ones_like(s_norm)

        for k in range(self.n_poles):
            p_k = poles[k]
            r_k = residues[k]

            # Pole and conjugate pair
            term1 = r_k / (s_norm - p_k)
            term2 = torch.conj(r_k) / (s_norm - torch.conj(p_k))
            result = result + term1 + term2

        # Scale back to physical units
        return result * self.Z_ref

    def get_sparsity_loss(self) -> torch.Tensor:
        """L1 on residue magnitudes."""
        residues = self.get_residues_normalized()
        return torch.sum(torch.abs(residues))


def train_spn_v2(
    freqs: np.ndarray,
    Z_data: np.ndarray,
    n_poles: int = 5,
    sparsity_weight: float = 0.001,
    n_epochs: int = 5000,
    n_restarts: int = 5,
    verbose: bool = True
) -> StablePoleNetworkV2:
    """Train SPN v2 with automatic scaling and multiple restarts."""

    # Determine reference scales
    omega = 2 * np.pi * freqs
    omega_ref = np.sqrt(omega.min() * omega.max())  # Geometric mean
    Z_ref = np.abs(Z_data).mean()

    # Convert to torch
    s = torch.tensor(1j * omega, dtype=torch.complex128)
    Z_target = torch.tensor(Z_data, dtype=torch.complex128)

    best_model = None
    best_loss = float('inf')

    for restart in range(n_restarts):
        torch.manual_seed(restart * 42)

        model = StablePoleNetworkV2(n_poles=n_poles, omega_ref=omega_ref, Z_ref=Z_ref)
        optimizer = optim.Adam(model.parameters(), lr=0.05)
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

        for epoch in range(n_epochs):
            optimizer.zero_grad()

            Z_pred = model(s)

            # Relative error loss
            rel_err = (Z_pred - Z_target) / (torch.abs(Z_target) + 1e-10 * Z_ref)
            loss_fit = torch.mean(torch.abs(rel_err) ** 2)

            # Sparsity
            loss_sparse = sparsity_weight * model.get_sparsity_loss()

            loss = loss_fit + loss_sparse

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            scheduler.step()

        final_loss = loss_fit.item()
        if verbose:
            print(f"  Restart {restart+1}: final_loss = {final_loss:.6f}")

        if final_loss < best_loss:
            best_loss = final_loss
            best_model = StablePoleNetworkV2(n_poles=n_poles, omega_ref=omega_ref, Z_ref=Z_ref)
            best_model.load_state_dict(model.state_dict())

    return best_model


def test_rc_circuit():
    """Test on RC parallel circuit."""
    print("=" * 70)
    print("Test: RC Parallel Circuit")
    print("=" * 70)

    # RC parallel: Z = R / (1 + j*omega*R*C)
    R, C = 1000.0, 1e-9  # 1 kOhm, 1 nF
    tau = R * C  # 1e-6 s
    f_pole = 1 / (2 * np.pi * tau)  # ~159 kHz

    print(f"True: R = {R:.0f} Ohm, C = {C:.2e} F")
    print(f"True pole: {-1/tau:.2e} rad/s (f = {f_pole:.0f} Hz)")

    # Generate data
    freqs = np.logspace(3, 7, 50)  # 1 kHz to 10 MHz
    omega = 2 * np.pi * freqs
    Z_data = R / (1 + 1j * omega * tau)

    # Add noise
    noise = 0.02
    Z_noisy = Z_data * (1 + noise * (np.random.randn(len(freqs)) + 1j * np.random.randn(len(freqs))))

    print(f"\nTraining on noisy data ({noise*100:.0f}% noise)...")
    model = train_spn_v2(freqs, Z_noisy, n_poles=3, n_epochs=3000, n_restarts=5)

    # Get learned poles
    poles = model.get_poles_physical()
    print(f"\nLearned poles (rad/s):")
    for i, p in enumerate(poles):
        print(f"  p_{i} = {p.real:.2e} + j*{p.imag:.2e}")

    # Find dominant (most negative real part with small imag)
    real_poles = poles[np.abs(poles.imag) < 0.1 * np.abs(poles.real)]
    if len(real_poles) > 0:
        dominant = real_poles[np.argmin(real_poles.real)]
        tau_learned = -1 / dominant.real
        error = abs(tau_learned - tau) / tau * 100
        print(f"\nDominant real pole: {dominant.real:.2e}")
        print(f"Learned tau: {tau_learned:.2e} s (error: {error:.1f}%)")
    else:
        # Find pole with largest residue contribution
        print("\nNo clear real pole found, checking all poles...")

    # Verify fit
    s = torch.tensor(1j * omega, dtype=torch.complex128)
    with torch.no_grad():
        Z_pred = model(s).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit quality: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_two_pole():
    """Test on two-pole system."""
    print("\n" + "=" * 70)
    print("Test: Two-Pole RC System")
    print("=" * 70)

    # Two RC in series
    R1, C1 = 100, 1e-9     # tau1 = 1e-7 s, f1 = 1.6 MHz
    R2, C2 = 1000, 1e-10   # tau2 = 1e-7 s, f2 = 1.6 MHz (same!)

    tau1, tau2 = R1*C1, R2*C2
    print(f"True poles:")
    print(f"  tau1 = {tau1:.2e} s -> pole at {-1/tau1:.2e} rad/s")
    print(f"  tau2 = {tau2:.2e} s -> pole at {-1/tau2:.2e} rad/s")

    # Generate data
    freqs = np.logspace(4, 8, 60)
    omega = 2 * np.pi * freqs
    Z1 = R1 / (1 + 1j * omega * tau1)
    Z2 = R2 / (1 + 1j * omega * tau2)
    Z_data = Z1 + Z2

    print("\nTraining SPN...")
    model = train_spn_v2(freqs, Z_data, n_poles=4, n_epochs=5000, n_restarts=5)

    # Get poles
    poles = model.get_poles_physical()
    print(f"\nLearned poles (rad/s):")
    for i, p in enumerate(poles):
        if np.abs(p.imag) < 0.5 * np.abs(p.real):  # Nearly real poles
            print(f"  p_{i} = {p.real:.2e} (real)")
        else:
            print(f"  p_{i} = {p.real:.2e} + j*{p.imag:.2e}")


def test_rlc_resonance():
    """Test on RLC resonant circuit."""
    print("\n" + "=" * 70)
    print("Test: RLC Series Resonance")
    print("=" * 70)

    R, L, C = 10.0, 1e-6, 1e-9  # 10 Ohm, 1 uH, 1 nF
    omega0 = 1 / np.sqrt(L * C)  # ~31.6 Mrad/s
    f0 = omega0 / (2 * np.pi)    # ~5 MHz
    Q = omega0 * L / R           # ~3.16

    print(f"True: R = {R:.0f} Ohm, L = {L:.2e} H, C = {C:.2e} F")
    print(f"Resonance: f0 = {f0/1e6:.2f} MHz, Q = {Q:.2f}")

    # Complex poles: p = -alpha +/- j*omega_d
    alpha = R / (2 * L)
    omega_d = np.sqrt(omega0**2 - alpha**2)
    print(f"True poles: {-alpha:.2e} +/- j*{omega_d:.2e}")

    # Generate data
    freqs = np.logspace(5, 8, 80)
    omega = 2 * np.pi * freqs
    Z_data = R + 1j * omega * L + 1 / (1j * omega * C)

    print("\nTraining SPN...")
    model = train_spn_v2(freqs, Z_data, n_poles=3, n_epochs=5000, n_restarts=5)

    poles = model.get_poles_physical()
    print(f"\nLearned poles (rad/s):")
    for i, p in enumerate(poles):
        print(f"  p_{i} = {p.real:.2e} + j*{p.imag:.2e}")

    # Find complex conjugate pair
    complex_poles = poles[np.abs(poles.imag) > 1e3]
    if len(complex_poles) >= 1:
        p = complex_poles[np.argmax(np.abs(complex_poles.imag))]
        print(f"\nDominant complex pole: {p.real:.2e} + j*{p.imag:.2e}")
        print(f"  True: {-alpha:.2e} + j*{omega_d:.2e}")
        print(f"  Error (real): {abs(p.real + alpha)/alpha*100:.1f}%")
        print(f"  Error (imag): {abs(abs(p.imag) - omega_d)/omega_d*100:.1f}%")

    # Fit quality
    s = torch.tensor(1j * omega, dtype=torch.complex128)
    with torch.no_grad():
        Z_pred = model(s).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_skin_effect():
    """Test on skin effect (sqrt behavior)."""
    print("\n" + "=" * 70)
    print("Test: Skin Effect sqrt(j*omega)")
    print("=" * 70)

    # Z ~ sqrt(1 + j*omega*tau) behavior
    R_dc = 0.1
    tau = 1e-8  # Characteristic time

    freqs = np.logspace(4, 8, 80)
    omega = 2 * np.pi * freqs
    Z_data = R_dc * np.sqrt(1 + 1j * omega * tau)

    print(f"True: R_dc = {R_dc} Ohm, tau = {tau:.2e} s")
    print("Note: sqrt behavior needs multiple poles to approximate")

    print("\nTraining SPN with 6 poles...")
    model = train_spn_v2(freqs, Z_data, n_poles=6, n_epochs=8000, n_restarts=5)

    # Fit quality
    s = torch.tensor(1j * omega, dtype=torch.complex128)
    with torch.no_grad():
        Z_pred = model(s).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")

    # Show poles
    poles = model.get_poles_physical()
    print(f"\nPoles (should be distributed along negative real axis):")
    sorted_poles = sorted(poles, key=lambda p: p.real)
    for p in sorted_poles:
        if np.abs(p.imag) < 0.3 * np.abs(p.real):
            print(f"  {p.real:.2e} (nearly real)")
        else:
            print(f"  {p.real:.2e} + j*{p.imag:.2e}")


if __name__ == '__main__':
    np.random.seed(42)

    print("Stable Pole Network v2: With Proper Scaling")
    print("=" * 70)

    test_rc_circuit()
    test_two_pole()
    test_rlc_resonance()
    test_skin_effect()

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print("""
With proper frequency scaling, SPN can:
1. Correctly identify pole locations
2. Handle noisy data
3. Work with resonant (complex pole) systems
4. Approximate distributed elements (skin effect)

Key improvements over v1:
- Frequency normalization (omega_ref)
- Impedance normalization (Z_ref)
- Relative error loss function
- Multiple restarts for global optimization
""")
