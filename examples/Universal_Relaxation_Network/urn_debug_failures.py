"""
Debug URN Failure Cases

Analyze why certain test cases fail and identify improvements.
"""

import numpy as np
import torch
import matplotlib.pyplot as plt
from universal_relaxation_network import URNConfig, train_urn

# =============================================================================
# FAILURE CASE 1: NiZn_ColeCole (10.11% error - just above threshold)
# =============================================================================

def debug_nizn_ferrite():
    """Debug NiZn ferrite Cole-Cole fitting."""
    print("="*70)
    print("DEBUGGING: NiZn Cole-Cole Ferrite")
    print("="*70)

    # True model
    mu_s = 500
    mu_inf = 1
    f0 = 10e6
    alpha = 0.85
    tau = 1 / (2 * np.pi * f0)

    freqs = np.logspace(3, 9, 80)
    omega = 2 * np.pi * freqs
    mu_true = mu_inf + (mu_s - mu_inf) / (1 + (1j * omega * tau)**alpha)

    print(f"True parameters: mu_s={mu_s}, mu_inf={mu_inf}, f0={f0/1e6}MHz, alpha={alpha}")
    print(f"tau = {tau:.2e} s")

    # Try different configurations
    configs = [
        ("Default", URNConfig(n_debye=3, n_cole_cole=2, n_cole_davidson=1,
                              n_havriliak_negami=1, n_cpe=1, sparsity_weight=0.005,
                              n_epochs=3000, n_restarts=3)),
        ("More epochs", URNConfig(n_debye=3, n_cole_cole=2, n_cole_davidson=1,
                                  n_havriliak_negami=1, n_cpe=1, sparsity_weight=0.005,
                                  n_epochs=6000, n_restarts=5)),
        ("Less sparsity", URNConfig(n_debye=3, n_cole_cole=2, n_cole_davidson=1,
                                    n_havriliak_negami=1, n_cpe=1, sparsity_weight=0.001,
                                    n_epochs=4000, n_restarts=3)),
        ("More Cole-Cole", URNConfig(n_debye=2, n_cole_cole=4, n_cole_davidson=1,
                                     n_havriliak_negami=1, n_cpe=1, sparsity_weight=0.003,
                                     n_epochs=4000, n_restarts=3)),
        ("More H-N", URNConfig(n_debye=2, n_cole_cole=2, n_cole_davidson=2,
                               n_havriliak_negami=2, n_cpe=1, sparsity_weight=0.003,
                               n_epochs=4000, n_restarts=3)),
    ]

    best_error = float('inf')
    best_config_name = None

    for name, config in configs:
        print(f"\n--- Config: {name} ---")
        config.n_warburg = 0
        config.n_gerischer = 0
        config.n_rlc = 0
        config.n_skin_effect = 0

        model = train_urn(freqs, mu_true, config, verbose=False)

        with torch.no_grad():
            mu_pred = model(torch.tensor(omega, dtype=torch.float64)).numpy()

        rel_err = np.abs(mu_pred - mu_true) / np.abs(mu_true)
        max_err = rel_err.max() * 100
        mean_err = rel_err.mean() * 100

        print(f"  Max error: {max_err:.2f}%, Mean: {mean_err:.2f}%")

        if max_err < best_error:
            best_error = max_err
            best_config_name = name

        # Show active components
        active = model.get_active_components()
        print(f"  Active: {list(active.keys())}")

        # Check if Cole-Cole with correct alpha was found
        if 'cole_cole' in active:
            for cc in active['cole_cole']:
                alpha_learned = cc.get('alpha', 0)
                print(f"    Cole-Cole alpha: {alpha_learned:.3f} (true: {alpha})")

    print(f"\n*** Best config: {best_config_name} with {best_error:.2f}% max error ***")
    return best_error


# =============================================================================
# FAILURE CASE 2: Cu_foil_0.5mm (23.28% error)
# =============================================================================

def debug_cu_foil():
    """Debug copper foil skin effect fitting."""
    print("\n" + "="*70)
    print("DEBUGGING: Cu Foil 0.5mm Skin Effect")
    print("="*70)

    # True model
    R_dc = 0.002  # 2 mOhm
    d = 0.5e-3   # 0.5 mm
    sigma = 5.8e7
    mu_0 = 4 * np.pi * 1e-7

    freqs = np.logspace(1, 6, 80)
    omega = 2 * np.pi * freqs
    delta = np.sqrt(2 / (omega * mu_0 * sigma + 1e-20))
    z = (1 + 1j) * d / (2 * delta)
    Z_true = R_dc * z / np.tanh(z + 1e-10)

    print(f"True parameters: R_dc={R_dc*1000:.1f} mOhm, d={d*1000:.1f} mm")
    print(f"Skin depth at 1kHz: {np.sqrt(2/(2*np.pi*1000*mu_0*sigma))*1000:.2f} mm")
    print(f"Skin depth at 1MHz: {np.sqrt(2/(2*np.pi*1e6*mu_0*sigma))*1000:.3f} mm")

    # Analyze the problem
    print(f"\nZ range: {np.abs(Z_true).min():.2e} to {np.abs(Z_true).max():.2e}")
    print(f"Z ratio (max/min): {np.abs(Z_true).max()/np.abs(Z_true).min():.1f}x")

    # Try different configurations
    configs = [
        ("Default", URNConfig(n_debye=2, n_cole_cole=1, n_cpe=2, n_warburg=1,
                              n_skin_effect=2, sparsity_weight=0.01,
                              n_epochs=4000, n_restarts=3)),
        ("More epochs", URNConfig(n_debye=2, n_cole_cole=1, n_cpe=2, n_warburg=1,
                                  n_skin_effect=2, sparsity_weight=0.01,
                                  n_epochs=8000, n_restarts=5)),
        ("Less sparsity", URNConfig(n_debye=2, n_cole_cole=1, n_cpe=2, n_warburg=1,
                                    n_skin_effect=2, sparsity_weight=0.001,
                                    n_epochs=5000, n_restarts=3)),
        ("More skin effect", URNConfig(n_debye=1, n_cole_cole=1, n_cpe=1, n_warburg=0,
                                       n_skin_effect=4, sparsity_weight=0.005,
                                       n_epochs=5000, n_restarts=3)),
        ("CPE focus", URNConfig(n_debye=2, n_cole_cole=1, n_cpe=4, n_warburg=1,
                                n_skin_effect=2, sparsity_weight=0.005,
                                n_epochs=5000, n_restarts=3)),
        ("Lower LR", URNConfig(n_debye=2, n_cole_cole=1, n_cpe=2, n_warburg=1,
                               n_skin_effect=3, sparsity_weight=0.005,
                               lr=0.01, n_epochs=6000, n_restarts=4)),
    ]

    best_error = float('inf')
    best_config_name = None

    for name, config in configs:
        print(f"\n--- Config: {name} ---")
        config.n_gerischer = 0
        config.n_rlc = 0
        config.n_cole_davidson = 0
        config.n_havriliak_negami = 0

        model = train_urn(freqs, Z_true, config, verbose=False)

        with torch.no_grad():
            Z_pred = model(torch.tensor(omega, dtype=torch.float64)).numpy()

        rel_err = np.abs(Z_pred - Z_true) / np.abs(Z_true)
        max_err = rel_err.max() * 100
        mean_err = rel_err.mean() * 100

        print(f"  Max error: {max_err:.2f}%, Mean: {mean_err:.2f}%")

        if max_err < best_error:
            best_error = max_err
            best_config_name = name

        # Show active components
        active = model.get_active_components()
        print(f"  Active: {list(active.keys())}")

        # Check skin effect parameters
        if 'skin_effect' in active:
            for se in active['skin_effect']:
                R_learned = se.get('R_dc', 0)
                delta_learned = se.get('delta', 0)
                print(f"    Skin: R_dc={R_learned*1000:.2f} mOhm, delta={delta_learned*1000:.3f} mm")

    print(f"\n*** Best config: {best_config_name} with {best_error:.2f}% max error ***")
    return best_error


# =============================================================================
# ANALYSIS: Why is thick conductor harder?
# =============================================================================

def analyze_skin_effect_difficulty():
    """Analyze why thick conductors are harder to fit."""
    print("\n" + "="*70)
    print("ANALYSIS: Skin Effect vs Conductor Thickness")
    print("="*70)

    mu_0 = 4 * np.pi * 1e-7
    sigma = 5.8e7

    thicknesses = [0.1e-3, 0.2e-3, 0.5e-3, 1e-3, 2e-3]
    R_dc_base = 0.01  # Base R_dc for 0.1mm

    print(f"\n{'Thickness':<12} {'R_dc':<10} {'delta@1kHz':<12} {'delta@1MHz':<12} {'d/delta@1MHz':<12}")
    print("-"*60)

    for d in thicknesses:
        R_dc = R_dc_base * (0.1e-3 / d)  # Scale R_dc with thickness
        delta_1k = np.sqrt(2/(2*np.pi*1000*mu_0*sigma))
        delta_1M = np.sqrt(2/(2*np.pi*1e6*mu_0*sigma))
        ratio = d / delta_1M

        print(f"{d*1000:>8.1f} mm  {R_dc*1000:>6.1f} mOhm  {delta_1k*1000:>8.2f} mm  {delta_1M*1000:>8.3f} mm  {ratio:>8.1f}")

    print("""
Observation:
- Thicker conductors have larger d/delta ratio at high frequency
- This means more nonlinear z*coth(z) behavior
- The transition from DC resistance to skin effect is sharper
- URN needs more basis functions to capture this transition

Solution approach:
- Use more skin_effect basis functions for thick conductors
- Or use transmission line basis for better approximation
""")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    err1 = debug_nizn_ferrite()
    err2 = debug_cu_foil()
    analyze_skin_effect_difficulty()

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"NiZn Cole-Cole best: {err1:.2f}%")
    print(f"Cu foil 0.5mm best: {err2:.2f}%")

    if err1 < 10 and err2 < 20:
        print("\n*** SUCCESS: Both cases improved to acceptable levels! ***")
    else:
        print("\n*** Further tuning needed ***")
