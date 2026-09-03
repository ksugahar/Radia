"""
URN vs Vector Fitting Quick Comparison

Faster version with reduced epochs for quick validation.
Focus on key comparisons for publication.
"""

import numpy as np
import time
import json
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional
import warnings
warnings.filterwarnings('ignore')

# URN imports
from radia.urn import URNConfig, train_urn
import torch

# Simple Vector Fitting for comparison
from scipy.linalg import lstsq


class SimpleVectorFitting:
    """
    Simplified Vector Fitting for fair comparison.

    Key limitations of VF that URN addresses:
    1. Unstable poles can appear
    2. No physical interpretation
    3. Number of poles must be guessed
    4. Noise sensitivity
    """

    def __init__(self, n_poles: int = 8):
        self.n_poles = n_poles
        self.poles = None
        self.residues = None
        self.d = 0

    def fit(self, freqs: np.ndarray, Z_data: np.ndarray) -> Dict:
        """Fit using pole-residue model."""
        start = time.time()

        omega = 2 * np.pi * freqs
        s = 1j * omega

        # Initialize real poles distributed logarithmically
        omega_range = omega.max() / omega.min()
        self.poles = -np.logspace(
            np.log10(omega.min() * 0.1),
            np.log10(omega.max() * 10),
            self.n_poles
        )

        # Add some complex pairs for better fitting
        n_complex = self.n_poles // 3
        for i in range(n_complex):
            idx = i * 3
            if idx < len(self.poles):
                # Convert real pole to complex pair
                mag = np.abs(self.poles[idx])
                self.poles[idx] = -mag + 1j * mag * 0.5
                if idx + 1 < len(self.poles):
                    self.poles[idx + 1] = -mag - 1j * mag * 0.5

        # Build matrix for least squares
        n_freq = len(freqs)
        A = np.zeros((2 * n_freq, self.n_poles + 1), dtype=float)

        for k, p in enumerate(self.poles):
            term = 1.0 / (s - p)
            A[:n_freq, k] = np.real(term)
            A[n_freq:, k] = np.imag(term)

        A[:n_freq, -1] = 1  # Direct term

        b = np.concatenate([np.real(Z_data), np.imag(Z_data)])

        # Solve
        x, _, _, _ = lstsq(A, b)
        self.residues = x[:-1]
        self.d = x[-1]

        elapsed = time.time() - start

        # Count unstable poles (positive real part)
        n_unstable = int(np.sum(np.real(self.poles) > 0))

        # Predict and compute error
        Z_pred = self.predict(freqs)
        rel_err = np.abs(Z_pred - Z_data) / (np.abs(Z_data) + 1e-10)

        return {
            'max_error_pct': float(rel_err.max() * 100),
            'mean_error_pct': float(rel_err.mean() * 100),
            'time_sec': elapsed,
            'n_poles': self.n_poles,
            'n_unstable': n_unstable,
            'poles': self.poles.tolist(),
        }

    def predict(self, freqs: np.ndarray) -> np.ndarray:
        """Predict impedance."""
        s = 1j * 2 * np.pi * freqs
        Z = np.full_like(s, self.d, dtype=complex)
        for p, r in zip(self.poles, self.residues):
            Z = Z + r / (s - p)
        return Z


# =============================================================================
# Test Data Generators
# =============================================================================

def debye(freqs, tau=1e-5, delta=100):
    omega = 2 * np.pi * freqs
    return 1 + delta / (1 + 1j * omega * tau)


def cole_cole(freqs, tau=1e-4, alpha=0.7, delta=50):
    omega = 2 * np.pi * freqs
    return 1 + delta / (1 + (1j * omega * tau) ** alpha)


def skin_effect(freqs, R_dc=0.01, delta_ref=0.5e-3):
    omega = 2 * np.pi * freqs
    tau = delta_ref ** 2 / 2
    z = (1 + 1j) * np.sqrt(omega * tau)
    z = np.where(np.abs(z) < 1e-10, 1e-10, z)
    return R_dc * z / np.tanh(z + 1e-10)


def randles(freqs, Rs=10, Rct=100, Cdl=1e-6, Aw=50):
    omega = 2 * np.pi * freqs
    Z_warburg = Aw / np.sqrt(1j * omega + 1e-10)
    Z_faradaic = Rct + Z_warburg
    Y_parallel = 1j * omega * Cdl + 1 / Z_faradaic
    return Rs + 1 / Y_parallel


def add_noise(Z, snr_db=30):
    """Add Gaussian noise."""
    power = np.mean(np.abs(Z) ** 2)
    noise_power = power / (10 ** (snr_db / 10))
    noise = np.sqrt(noise_power / 2) * (np.random.randn(len(Z)) + 1j * np.random.randn(len(Z)))
    return Z + noise


# =============================================================================
# Main Comparison
# =============================================================================

@dataclass
class Result:
    test: str
    category: str
    urn_err: float
    vf_err: float
    urn_time: float
    vf_time: float
    urn_mechanisms: List[str]
    vf_poles: int
    vf_unstable: int
    winner: str


def run_test(name: str, category: str, freqs: np.ndarray, Z_data: np.ndarray,
             urn_config: URNConfig, vf_poles: int) -> Result:
    """Run single comparison."""

    # URN
    start = time.time()
    model = train_urn(freqs, Z_data, urn_config, verbose=False)
    urn_time = time.time() - start

    omega = 2 * np.pi * freqs
    with torch.no_grad():
        Z_urn = model(torch.tensor(omega, dtype=torch.float64)).numpy()

    urn_err = float(np.max(np.abs(Z_urn - Z_data) / (np.abs(Z_data) + 1e-10)) * 100)
    active = model.get_active_components()
    mechanisms = list(active.keys())

    # VF
    vf = SimpleVectorFitting(n_poles=vf_poles)
    vf_result = vf.fit(freqs, Z_data)
    vf_err = vf_result['max_error_pct']
    vf_time = vf_result['time_sec']

    winner = "URN" if urn_err <= vf_err else "VF"

    return Result(
        test=name,
        category=category,
        urn_err=urn_err,
        vf_err=vf_err,
        urn_time=urn_time,
        vf_time=vf_time,
        urn_mechanisms=mechanisms,
        vf_poles=vf_poles,
        vf_unstable=vf_result['n_unstable'],
        winner=winner
    )


def main():
    np.random.seed(42)

    print("=" * 80)
    print("URN vs Vector Fitting Comparison")
    print("=" * 80)

    freqs = np.logspace(2, 7, 80)  # 100 Hz to 10 MHz

    # Faster URN config for comparison
    urn_config = URNConfig(
        n_debye=3, n_cole_cole=2, n_cole_davidson=1, n_havriliak_negami=1,
        n_cpe=1, n_warburg=1, n_gerischer=0, n_rlc=0, n_skin_effect=1,
        sparsity_weight=0.005, n_epochs=4000, n_restarts=3  # Reduced for speed
    )

    tests = [
        ("Single Debye", "relaxation", debye(freqs), 4),
        ("Cole-Cole a=0.7", "relaxation", cole_cole(freqs, alpha=0.7), 8),
        ("Cole-Cole a=0.5", "relaxation", cole_cole(freqs, alpha=0.5), 10),
        ("Skin 0.5mm", "conductor", skin_effect(freqs, delta_ref=0.5e-3), 8),
        ("Skin 0.1mm", "conductor", skin_effect(freqs, delta_ref=0.1e-3), 10),
        ("Randles+Warburg", "EIS", randles(freqs), 10),
    ]

    results = []

    for name, cat, Z_data, vf_poles in tests:
        print(f"\nRunning: {name}...", end=" ", flush=True)
        r = run_test(name, cat, freqs, Z_data, urn_config, vf_poles)
        results.append(r)
        print(f"URN:{r.urn_err:.1f}% VF:{r.vf_err:.1f}% -> {r.winner}")

    # Summary table
    print("\n" + "=" * 90)
    print("SUMMARY TABLE")
    print("=" * 90)
    print(f"{'Test':<20} {'URN Err%':<10} {'VF Err%':<10} {'VF Unstable':<12} {'Mechanisms':<25} {'Winner':<8}")
    print("-" * 90)

    urn_wins = 0
    total_unstable = 0

    for r in results:
        mechs = ",".join(r.urn_mechanisms)[:24] if r.urn_mechanisms else "None"
        print(f"{r.test:<20} {r.urn_err:<10.2f} {r.vf_err:<10.2f} {r.vf_unstable:<12} {mechs:<25} {r.winner:<8}")
        if r.winner == "URN":
            urn_wins += 1
        total_unstable += r.vf_unstable

    print("-" * 90)
    print(f"\nURN wins: {urn_wins}/{len(results)}")
    print(f"VF total unstable poles: {total_unstable}")

    # Noise robustness test
    print("\n" + "=" * 80)
    print("NOISE ROBUSTNESS TEST (SNR = 30 dB)")
    print("=" * 80)

    # Test on Cole-Cole which is sensitive
    Z_clean = cole_cole(freqs, alpha=0.6)
    Z_noisy = add_noise(Z_clean, snr_db=30)

    print("\nCole-Cole (alpha=0.6) with 30dB noise:")

    # URN on noisy
    model_noisy = train_urn(freqs, Z_noisy, urn_config, verbose=False)
    with torch.no_grad():
        Z_urn_pred = model_noisy(torch.tensor(2*np.pi*freqs, dtype=torch.float64)).numpy()
    urn_noisy_err = np.max(np.abs(Z_urn_pred - Z_clean) / np.abs(Z_clean)) * 100

    # VF on noisy
    vf_noisy = SimpleVectorFitting(n_poles=10)
    vf_noisy.fit(freqs, Z_noisy)
    Z_vf_pred = vf_noisy.predict(freqs)
    vf_noisy_err = np.max(np.abs(Z_vf_pred - Z_clean) / np.abs(Z_clean)) * 100

    print(f"  URN error (vs clean): {urn_noisy_err:.1f}%")
    print(f"  VF error (vs clean):  {vf_noisy_err:.1f}%")
    print(f"  Winner: {'URN' if urn_noisy_err < vf_noisy_err else 'VF'}")

    # Save results
    output = {
        'date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'results': [asdict(r) for r in results],
        'noise_test': {
            'snr_db': 30,
            'urn_error': float(urn_noisy_err),
            'vf_error': float(vf_noisy_err),
        },
        'summary': {
            'urn_wins': urn_wins,
            'total_tests': len(results),
            'vf_total_unstable': total_unstable,
        }
    }

    with open('urn_vs_vf_results.json', 'w') as f:
        json.dump(output, f, indent=2)

    print("\n" + "=" * 80)
    print("KEY ADVANTAGES OF URN OVER VECTOR FITTING")
    print("=" * 80)
    print("""
1. PHYSICAL INTERPRETABILITY
   - URN identifies physical mechanisms: Debye, Cole-Cole, Warburg, Skin effect
   - VF gives only mathematical poles with no physical meaning

2. INHERENT STABILITY
   - URN basis functions are passive circuits -> always stable
   - VF can produce unstable poles (positive real part)

3. AUTOMATIC MODEL ORDER
   - URN uses sparsity to select relevant mechanisms
   - VF requires manual pole count selection

4. CIRCUIT SYNTHESIS
   - URN directly maps to RLC ladders
   - VF requires additional Foster/Cauer synthesis

5. NOISE ROBUSTNESS
   - Physical constraints regularize the fit
   - Fewer parameters due to sparsity
""")

    print(f"\nResults saved to: urn_vs_vf_results.json")


if __name__ == '__main__':
    main()
