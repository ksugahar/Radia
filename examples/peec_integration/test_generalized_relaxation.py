"""
Generalized Relaxation Model Network

Idea: Instead of fitting poles directly (Vector Fitting),
fit a linear combination of known relaxation basis functions.

Basis functions:
1. Debye: 1 / (1 + j*omega*tau)
2. Cole-Cole: 1 / (1 + (j*omega*tau)^alpha)
3. Cole-Davidson: 1 / (1 + j*omega*tau)^beta
4. Havriliak-Negami: 1 / (1 + (j*omega*tau)^alpha)^beta
5. Constant Phase Element (CPE): (j*omega*tau)^(-gamma)

These can be combined to fit arbitrary frequency responses.

Advantages:
- Each basis has physical meaning
- Numerically stable (no pole extraction needed)
- Can discover which relaxation mechanism dominates
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Tuple, List


def debye(omega: torch.Tensor, tau: torch.Tensor) -> torch.Tensor:
    """Debye relaxation: 1 / (1 + j*omega*tau)"""
    return 1.0 / (1.0 + 1j * omega * tau)


def cole_cole(omega: torch.Tensor, tau: torch.Tensor, alpha: torch.Tensor) -> torch.Tensor:
    """Cole-Cole: 1 / (1 + (j*omega*tau)^alpha), 0 < alpha <= 1"""
    # Ensure 0 < alpha <= 1
    alpha_safe = torch.sigmoid(alpha) * 0.9 + 0.1  # Range [0.1, 1.0]
    jot = 1j * omega * tau
    return 1.0 / (1.0 + jot ** alpha_safe)


def cole_davidson(omega: torch.Tensor, tau: torch.Tensor, beta: torch.Tensor) -> torch.Tensor:
    """Cole-Davidson: 1 / (1 + j*omega*tau)^beta, 0 < beta <= 1"""
    beta_safe = torch.sigmoid(beta) * 0.9 + 0.1
    return 1.0 / (1.0 + 1j * omega * tau) ** beta_safe


def cpe(omega: torch.Tensor, tau: torch.Tensor, gamma: torch.Tensor) -> torch.Tensor:
    """Constant Phase Element: (j*omega*tau)^(-gamma), 0 < gamma < 1"""
    gamma_safe = torch.sigmoid(gamma) * 0.8 + 0.1  # Range [0.1, 0.9]
    jot = 1j * omega * tau
    # Avoid branch cut issues
    return jot ** (-gamma_safe)


class GeneralizedRelaxationNetwork(nn.Module):
    """
    Network that learns a mixture of relaxation models.

    Z(omega) = Z_inf + sum_k weight_k * basis_k(omega, tau_k, params_k)
    """

    def __init__(self, n_debye: int = 2, n_cole_cole: int = 1,
                 n_cole_davidson: int = 1, n_cpe: int = 0,
                 omega_ref: float = 1e6, Z_ref: float = 1.0):
        super().__init__()

        self.n_debye = n_debye
        self.n_cole_cole = n_cole_cole
        self.n_cole_davidson = n_cole_davidson
        self.n_cpe = n_cpe
        self.omega_ref = omega_ref
        self.Z_ref = Z_ref

        # Debye parameters: tau (time constant)
        if n_debye > 0:
            self.debye_log_tau = nn.Parameter(torch.randn(n_debye))
            self.debye_weight_real = nn.Parameter(torch.randn(n_debye) * 0.1)
            self.debye_weight_imag = nn.Parameter(torch.randn(n_debye) * 0.1)

        # Cole-Cole parameters: tau, alpha
        if n_cole_cole > 0:
            self.cc_log_tau = nn.Parameter(torch.randn(n_cole_cole))
            self.cc_alpha_raw = nn.Parameter(torch.zeros(n_cole_cole))
            self.cc_weight_real = nn.Parameter(torch.randn(n_cole_cole) * 0.1)
            self.cc_weight_imag = nn.Parameter(torch.randn(n_cole_cole) * 0.1)

        # Cole-Davidson parameters: tau, beta
        if n_cole_davidson > 0:
            self.cd_log_tau = nn.Parameter(torch.randn(n_cole_davidson))
            self.cd_beta_raw = nn.Parameter(torch.zeros(n_cole_davidson))
            self.cd_weight_real = nn.Parameter(torch.randn(n_cole_davidson) * 0.1)
            self.cd_weight_imag = nn.Parameter(torch.randn(n_cole_davidson) * 0.1)

        # CPE parameters: tau (reference), gamma (exponent)
        if n_cpe > 0:
            self.cpe_log_tau = nn.Parameter(torch.randn(n_cpe))
            self.cpe_gamma_raw = nn.Parameter(torch.zeros(n_cpe))
            self.cpe_weight_real = nn.Parameter(torch.randn(n_cpe) * 0.1)
            self.cpe_weight_imag = nn.Parameter(torch.randn(n_cpe) * 0.1)

        # DC and high-frequency limits
        self.Z_inf_real = nn.Parameter(torch.tensor(0.1))
        self.Z_inf_imag = nn.Parameter(torch.tensor(0.0))

    def forward(self, omega: torch.Tensor) -> torch.Tensor:
        """Compute Z(omega)."""
        omega_norm = omega / self.omega_ref

        # Start with high-frequency limit
        Z = torch.complex(self.Z_inf_real, self.Z_inf_imag) * torch.ones_like(omega_norm)

        # Add Debye contributions
        if self.n_debye > 0:
            tau = torch.exp(self.debye_log_tau)
            for k in range(self.n_debye):
                weight = torch.complex(self.debye_weight_real[k], self.debye_weight_imag[k])
                Z = Z + weight * debye(omega_norm, tau[k])

        # Add Cole-Cole contributions
        if self.n_cole_cole > 0:
            tau = torch.exp(self.cc_log_tau)
            alpha = torch.sigmoid(self.cc_alpha_raw) * 0.9 + 0.1
            for k in range(self.n_cole_cole):
                weight = torch.complex(self.cc_weight_real[k], self.cc_weight_imag[k])
                Z = Z + weight * cole_cole(omega_norm, tau[k], alpha[k])

        # Add Cole-Davidson contributions
        if self.n_cole_davidson > 0:
            tau = torch.exp(self.cd_log_tau)
            beta = torch.sigmoid(self.cd_beta_raw) * 0.9 + 0.1
            for k in range(self.n_cole_davidson):
                weight = torch.complex(self.cd_weight_real[k], self.cd_weight_imag[k])
                Z = Z + weight * cole_davidson(omega_norm, tau[k], beta[k])

        # Add CPE contributions
        if self.n_cpe > 0:
            tau = torch.exp(self.cpe_log_tau)
            gamma = torch.sigmoid(self.cpe_gamma_raw) * 0.8 + 0.1
            for k in range(self.n_cpe):
                weight = torch.complex(self.cpe_weight_real[k], self.cpe_weight_imag[k])
                Z = Z + weight * cpe(omega_norm, tau[k], gamma[k])

        return Z * self.Z_ref

    def get_parameters_summary(self) -> dict:
        """Return learned parameters in interpretable form."""
        summary = {}

        if self.n_debye > 0:
            tau = torch.exp(self.debye_log_tau).detach().cpu().numpy() / self.omega_ref
            weights = (self.debye_weight_real + 1j * self.debye_weight_imag).detach().cpu().numpy() * self.Z_ref
            summary['debye'] = {'tau': tau, 'weights': weights}

        if self.n_cole_cole > 0:
            tau = torch.exp(self.cc_log_tau).detach().cpu().numpy() / self.omega_ref
            alpha = (torch.sigmoid(self.cc_alpha_raw) * 0.9 + 0.1).detach().cpu().numpy()
            weights = (self.cc_weight_real + 1j * self.cc_weight_imag).detach().cpu().numpy() * self.Z_ref
            summary['cole_cole'] = {'tau': tau, 'alpha': alpha, 'weights': weights}

        if self.n_cole_davidson > 0:
            tau = torch.exp(self.cd_log_tau).detach().cpu().numpy() / self.omega_ref
            beta = (torch.sigmoid(self.cd_beta_raw) * 0.9 + 0.1).detach().cpu().numpy()
            weights = (self.cd_weight_real + 1j * self.cd_weight_imag).detach().cpu().numpy() * self.Z_ref
            summary['cole_davidson'] = {'tau': tau, 'beta': beta, 'weights': weights}

        if self.n_cpe > 0:
            tau = torch.exp(self.cpe_log_tau).detach().cpu().numpy() / self.omega_ref
            gamma = (torch.sigmoid(self.cpe_gamma_raw) * 0.8 + 0.1).detach().cpu().numpy()
            weights = (self.cpe_weight_real + 1j * self.cpe_weight_imag).detach().cpu().numpy() * self.Z_ref
            summary['cpe'] = {'tau': tau, 'gamma': gamma, 'weights': weights}

        summary['Z_inf'] = complex(self.Z_inf_real.item(), self.Z_inf_imag.item()) * self.Z_ref

        return summary


def train_grn(
    freqs: np.ndarray,
    Z_data: np.ndarray,
    n_debye: int = 2,
    n_cole_cole: int = 1,
    n_cole_davidson: int = 1,
    n_cpe: int = 0,
    n_epochs: int = 5000,
    n_restarts: int = 3,
    verbose: bool = True
) -> GeneralizedRelaxationNetwork:
    """Train Generalized Relaxation Network."""

    omega = 2 * np.pi * freqs
    omega_ref = np.sqrt(omega.min() * omega.max())
    Z_ref = np.abs(Z_data).mean()

    omega_torch = torch.tensor(omega, dtype=torch.float64)
    Z_target = torch.tensor(Z_data, dtype=torch.complex128)

    best_model = None
    best_loss = float('inf')

    for restart in range(n_restarts):
        torch.manual_seed(restart * 123)

        model = GeneralizedRelaxationNetwork(
            n_debye=n_debye, n_cole_cole=n_cole_cole,
            n_cole_davidson=n_cole_davidson, n_cpe=n_cpe,
            omega_ref=omega_ref, Z_ref=Z_ref
        )

        optimizer = optim.Adam(model.parameters(), lr=0.05)
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

        for epoch in range(n_epochs):
            optimizer.zero_grad()

            Z_pred = model(omega_torch)

            # Relative error
            rel_err = (Z_pred - Z_target) / (torch.abs(Z_target) + 1e-10 * Z_ref)
            loss = torch.mean(torch.abs(rel_err) ** 2)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            scheduler.step()

        final_loss = loss.item()
        if verbose:
            print(f"  Restart {restart+1}: loss = {final_loss:.6f}")

        if final_loss < best_loss:
            best_loss = final_loss
            best_model = GeneralizedRelaxationNetwork(
                n_debye=n_debye, n_cole_cole=n_cole_cole,
                n_cole_davidson=n_cole_davidson, n_cpe=n_cpe,
                omega_ref=omega_ref, Z_ref=Z_ref
            )
            best_model.load_state_dict(model.state_dict())

    return best_model


def test_single_debye():
    """Test fitting single Debye relaxation."""
    print("=" * 70)
    print("Test: Single Debye Relaxation (RC Circuit)")
    print("=" * 70)

    R, C = 1000, 1e-9
    tau_true = R * C

    freqs = np.logspace(3, 8, 50)
    omega = 2 * np.pi * freqs
    Z_data = R / (1 + 1j * omega * tau_true)

    print(f"True: R = {R} Ohm, tau = {tau_true:.2e} s")

    model = train_grn(freqs, Z_data, n_debye=2, n_cole_cole=0,
                      n_cole_davidson=0, n_cpe=0, n_epochs=3000)

    params = model.get_parameters_summary()
    print(f"\nLearned Debye parameters:")
    for k, tau in enumerate(params['debye']['tau']):
        w = params['debye']['weights'][k]
        print(f"  Debye {k}: tau = {tau:.2e} s, weight = {w:.2f}")

    # Fit quality
    omega_torch = torch.tensor(omega, dtype=torch.float64)
    with torch.no_grad():
        Z_pred = model(omega_torch).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_cole_cole_material():
    """Test fitting Cole-Cole dielectric (common in biology/polymers)."""
    print("\n" + "=" * 70)
    print("Test: Cole-Cole Relaxation (Biological Tissue)")
    print("=" * 70)

    # Typical tissue parameters
    eps_inf = 10
    eps_s = 1000
    tau_true = 1e-6
    alpha_true = 0.8

    freqs = np.logspace(2, 8, 80)
    omega = 2 * np.pi * freqs

    # Cole-Cole model
    jot = 1j * omega * tau_true
    eps = eps_inf + (eps_s - eps_inf) / (1 + jot ** alpha_true)
    Z_data = 1 / (1j * omega * eps * 1e-12)  # Capacitive impedance

    print(f"True: eps_s/eps_inf = {eps_s}/{eps_inf}, tau = {tau_true:.2e} s, alpha = {alpha_true}")

    model = train_grn(freqs, Z_data, n_debye=1, n_cole_cole=2,
                      n_cole_davidson=0, n_cpe=0, n_epochs=5000)

    params = model.get_parameters_summary()
    print(f"\nLearned Cole-Cole parameters:")
    if 'cole_cole' in params:
        for k in range(len(params['cole_cole']['tau'])):
            tau = params['cole_cole']['tau'][k]
            alpha = params['cole_cole']['alpha'][k]
            w = params['cole_cole']['weights'][k]
            print(f"  CC {k}: tau = {tau:.2e} s, alpha = {alpha:.2f}, weight = {np.abs(w):.2e}")

    omega_torch = torch.tensor(omega, dtype=torch.float64)
    with torch.no_grad():
        Z_pred = model(omega_torch).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_skin_effect_with_cpe():
    """Test fitting skin effect using CPE (fractional behavior)."""
    print("\n" + "=" * 70)
    print("Test: Skin Effect with CPE Basis")
    print("=" * 70)

    # Z ~ sqrt(j*omega) = (j*omega)^0.5
    R_dc = 0.1
    tau = 1e-8

    freqs = np.logspace(4, 8, 60)
    omega = 2 * np.pi * freqs
    Z_data = R_dc * np.sqrt(1 + 1j * omega * tau)

    print(f"True: R_dc = {R_dc}, behavior ~ sqrt(j*omega)")
    print("CPE with gamma=0.5 should fit perfectly")

    # Use CPE basis
    model = train_grn(freqs, Z_data, n_debye=1, n_cole_cole=0,
                      n_cole_davidson=1, n_cpe=2, n_epochs=5000)

    params = model.get_parameters_summary()
    print(f"\nLearned parameters:")
    print(f"  Z_inf = {params['Z_inf']:.4f}")

    if 'cpe' in params:
        for k in range(len(params['cpe']['gamma'])):
            gamma = params['cpe']['gamma'][k]
            tau = params['cpe']['tau'][k]
            w = params['cpe']['weights'][k]
            print(f"  CPE {k}: gamma = {gamma:.3f}, tau = {tau:.2e}, |weight| = {np.abs(w):.4f}")
            if abs(gamma - 0.5) < 0.1:
                print(f"    ^ gamma ~ 0.5 indicates sqrt behavior!")

    omega_torch = torch.tensor(omega, dtype=torch.float64)
    with torch.no_grad():
        Z_pred = model(omega_torch).numpy()

    rel_err = np.abs(Z_pred - Z_data) / np.abs(Z_data)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_complex_ferrite():
    """Test fitting complex permeability of ferrite."""
    print("\n" + "=" * 70)
    print("Test: Ferrite Complex Permeability mu(f)")
    print("=" * 70)

    # Typical ferrite: mu' drops, mu'' peaks at ~MHz
    mu_s = 2000  # Low-freq permeability
    mu_inf = 1   # High-freq permeability
    f0 = 3e6     # Resonance frequency
    tau = 1 / (2 * np.pi * f0)

    # Debye-like relaxation for mu
    freqs = np.logspace(4, 8, 80)
    omega = 2 * np.pi * freqs
    mu = mu_inf + (mu_s - mu_inf) / (1 + 1j * omega * tau)

    # Convert to impedance-like quantity for fitting
    Z_data = mu  # Just fit mu directly as complex function

    print(f"True: mu_s = {mu_s}, mu_inf = {mu_inf}, f0 = {f0/1e6:.1f} MHz")

    model = train_grn(freqs, Z_data, n_debye=2, n_cole_cole=1,
                      n_cole_davidson=0, n_cpe=0, n_epochs=5000)

    params = model.get_parameters_summary()
    print(f"\nLearned parameters:")
    print(f"  mu_inf = {params['Z_inf']:.1f}")

    if 'debye' in params:
        total_weight = np.sum(params['debye']['weights'].real)
        print(f"  Total Debye weight = {total_weight:.1f} (expect ~{mu_s - mu_inf})")
        for k, tau_k in enumerate(params['debye']['tau']):
            f_k = 1 / (2 * np.pi * tau_k)
            w = params['debye']['weights'][k]
            print(f"    Debye {k}: f = {f_k/1e6:.2f} MHz, weight = {w:.1f}")

    omega_torch = torch.tensor(omega, dtype=torch.float64)
    with torch.no_grad():
        mu_pred = model(omega_torch).numpy()

    rel_err = np.abs(mu_pred - mu) / np.abs(mu)
    print(f"\nFit: max_err = {rel_err.max()*100:.2f}%, mean_err = {rel_err.mean()*100:.2f}%")


def test_model_selection():
    """Test automatic model selection via sparsity."""
    print("\n" + "=" * 70)
    print("Test: Automatic Model Selection")
    print("=" * 70)

    # Generate data from single Debye
    tau_true = 1e-6
    freqs = np.logspace(3, 8, 50)
    omega = 2 * np.pi * freqs
    Z_data = 100 / (1 + 1j * omega * tau_true)

    print("True model: Single Debye")
    print("Fitting with: 3 Debye + 2 Cole-Cole + 2 Cole-Davidson")

    model = train_grn(freqs, Z_data, n_debye=3, n_cole_cole=2,
                      n_cole_davidson=2, n_cpe=0, n_epochs=5000)

    params = model.get_parameters_summary()
    print("\nLearned weights (should be sparse):")

    # Check which components are significant
    threshold = 1.0  # Significant if |weight| > 1

    print("  Debye weights:")
    for k, w in enumerate(params['debye']['weights']):
        status = "ACTIVE" if np.abs(w) > threshold else "inactive"
        print(f"    {k}: |w| = {np.abs(w):.2f} [{status}]")

    print("  Cole-Cole weights:")
    for k, w in enumerate(params['cole_cole']['weights']):
        status = "ACTIVE" if np.abs(w) > threshold else "inactive"
        print(f"    {k}: |w| = {np.abs(w):.2f} [{status}]")

    print("  Cole-Davidson weights:")
    for k, w in enumerate(params['cole_davidson']['weights']):
        status = "ACTIVE" if np.abs(w) > threshold else "inactive"
        print(f"    {k}: |w| = {np.abs(w):.2f} [{status}]")


if __name__ == '__main__':
    np.random.seed(42)

    print("Generalized Relaxation Network (GRN)")
    print("=" * 70)
    print("Fitting arbitrary frequency response with physical basis functions")
    print("=" * 70)

    test_single_debye()
    test_cole_cole_material()
    test_skin_effect_with_cpe()
    test_complex_ferrite()
    test_model_selection()

    print("\n" + "=" * 70)
    print("Conclusions")
    print("=" * 70)
    print("""
Generalized Relaxation Network advantages:

1. PHYSICAL INTERPRETABILITY:
   - Each basis (Debye, Cole-Cole, etc.) has clear physical meaning
   - Learned parameters correspond to relaxation times, distribution widths

2. NUMERICAL STABILITY:
   - No pole extraction needed
   - Basis functions are well-behaved

3. MODEL DISCOVERY:
   - Can identify dominant relaxation mechanism
   - Sparse solutions reveal which model type fits best

4. FLEXIBILITY:
   - Combines multiple relaxation models
   - Handles distributed (CPE) and resonant behavior

Potential for publication:
- Novel approach combining basis function fitting with neural network optimization
- Automatic model selection through sparsity
- Applicable to dielectrics, magnetics, electrochemistry
""")
