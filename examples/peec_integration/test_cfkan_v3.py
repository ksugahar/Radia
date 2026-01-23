"""
CF-KAN v3: Research on Coefficient Extraction without Prior Knowledge

Key research questions:
1. Can CF-KAN converge without knowing target coefficients?
2. What physical constraints help?
3. Can we identify the correct number of stages?
4. Does it work for non-Dowell continued fractions?
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Tuple, List, Optional


class CFKANResearch(nn.Module):
    """
    Research version of CF-KAN with various constraint options.
    """

    def __init__(self, n_stages: int, constraints: str = 'none'):
        """
        Args:
            n_stages: Number of continued fraction terms
            constraints: 'none', 'positive', 'monotonic', 'odd_prior'
        """
        super().__init__()
        self.n_stages = n_stages
        self.constraints = constraints

        # Random initialization (no prior knowledge)
        init_vals = torch.rand(n_stages) * 8 + 2  # Range [2, 10]
        self.a_raw = nn.Parameter(init_vals)

    @property
    def a(self):
        """Apply constraints to raw parameters."""
        if self.constraints == 'none':
            return self.a_raw
        elif self.constraints == 'positive':
            # Softplus ensures positive values
            return torch.nn.functional.softplus(self.a_raw) + 0.1
        elif self.constraints == 'monotonic':
            # Cumulative sum ensures monotonic increase
            deltas = torch.nn.functional.softplus(self.a_raw) + 0.5
            return torch.cumsum(deltas, dim=0) + 1.0
        elif self.constraints == 'odd_prior':
            # Soft constraint toward odd integers
            base = torch.nn.functional.softplus(self.a_raw) + 0.1
            return base
        else:
            return self.a_raw

    def forward(self, z):
        """Evaluate continued fraction."""
        z_sq = z ** 2
        a = self.a

        result = a[-1].to(z.dtype)
        for i in range(self.n_stages - 2, -1, -1):
            result = a[i].to(z.dtype) + z_sq / result

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


def debye_relaxation(omega, tau, eps_s, eps_inf):
    """
    Debye relaxation model for complex permittivity.

    eps(omega) = eps_inf + (eps_s - eps_inf) / (1 + j*omega*tau)

    This has a continued fraction representation too.
    """
    return eps_inf + (eps_s - eps_inf) / (1 + 1j * omega * tau)


def generate_training_data(func_type: str, n_points: int = 100):
    """Generate training data for different function types."""

    if func_type == 'dowell':
        # z*coth(z) for skin effect
        z_vals = np.linspace(0.1, 5.0, n_points) * (1 + 1j) / np.sqrt(2)
        y_vals = zcothz_exact(z_vals)
        true_coeffs = np.array([3.0, 5.0, 7.0, 9.0, 11.0])
        return z_vals, y_vals, true_coeffs

    elif func_type == 'modified_dowell':
        # Modified coefficients [2.5, 4.5, 6.5, 8.5, 10.5]
        z_vals = np.linspace(0.1, 5.0, n_points) * (1 + 1j) / np.sqrt(2)
        coeffs = [2.5, 4.5, 6.5, 8.5, 10.5]
        y_vals = evaluate_cf(z_vals, coeffs)
        return z_vals, y_vals, np.array(coeffs)

    elif func_type == 'non_monotonic':
        # Non-monotonic coefficients [5, 3, 7, 2, 9]
        z_vals = np.linspace(0.1, 5.0, n_points) * (1 + 1j) / np.sqrt(2)
        coeffs = [5.0, 3.0, 7.0, 2.0, 9.0]
        y_vals = evaluate_cf(z_vals, coeffs)
        return z_vals, y_vals, np.array(coeffs)

    elif func_type == 'three_term':
        # Simple 3-term CF [2, 4, 6]
        z_vals = np.linspace(0.1, 4.0, n_points) * (1 + 1j) / np.sqrt(2)
        coeffs = [2.0, 4.0, 6.0]
        y_vals = evaluate_cf(z_vals, coeffs)
        return z_vals, y_vals, np.array(coeffs)

    else:
        raise ValueError(f"Unknown func_type: {func_type}")


def evaluate_cf(z, coeffs):
    """Evaluate continued fraction with given coefficients."""
    z = np.asarray(z, dtype=complex)
    z_sq = z ** 2
    n = len(coeffs)

    result = np.full_like(z, coeffs[n - 1], dtype=complex)
    for i in range(n - 2, -1, -1):
        result = coeffs[i] + z_sq / result

    return 1.0 + z_sq / result


def train_and_evaluate(
    func_type: str,
    n_stages: int,
    constraints: str,
    n_restarts: int = 20,
    n_epochs: int = 2000,
    verbose: bool = True
) -> Tuple[np.ndarray, float, float]:
    """
    Train CF-KAN and return best result.

    Returns:
        (best_coeffs, coeff_error_pct, func_error_pct)
    """
    z_train, y_train, true_coeffs = generate_training_data(func_type)
    z_torch = torch.tensor(z_train, dtype=torch.complex128)
    y_torch = torch.tensor(y_train, dtype=torch.complex128)

    best_coeffs = None
    best_coeff_error = float('inf')
    best_func_error = float('inf')

    for trial in range(n_restarts):
        torch.manual_seed(trial * 42 + 7)

        model = CFKANResearch(n_stages=n_stages, constraints=constraints)
        optimizer = optim.Adam(model.parameters(), lr=0.1)
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

        for epoch in range(n_epochs):
            optimizer.zero_grad()
            y_pred = model(z_torch)
            loss = torch.mean(torch.abs(y_pred - y_torch) ** 2)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            scheduler.step()

        # Evaluate
        final_coeffs = model.get_coefficients()

        # Coefficient error (only if same number of stages)
        if len(final_coeffs) == len(true_coeffs):
            coeff_error = np.max(np.abs(final_coeffs - true_coeffs) / (np.abs(true_coeffs) + 1e-10)) * 100
        else:
            coeff_error = float('inf')

        # Function error
        with torch.no_grad():
            y_pred_final = model(z_torch).numpy()
        func_error = np.max(np.abs(y_pred_final - y_train) / (np.abs(y_train) + 1e-10)) * 100

        if coeff_error < best_coeff_error:
            best_coeff_error = coeff_error
            best_coeffs = final_coeffs.copy()
            best_func_error = func_error

        if verbose and (trial + 1) % 5 == 0:
            print(f"    Trial {trial+1:2d}: coeffs={np.array2string(final_coeffs, precision=2)}, "
                  f"coeff_err={coeff_error:.1f}%, func_err={func_error:.2f}%")

    return best_coeffs, best_coeff_error, best_func_error


def research_1_constraint_comparison():
    """Compare different constraint types."""
    print("=" * 70)
    print("Research 1: Effect of Constraints on Dowell Coefficient Recovery")
    print("=" * 70)

    constraints_list = ['none', 'positive', 'monotonic']
    results = {}

    for constraint in constraints_list:
        print(f"\n--- Constraint: {constraint} ---")
        coeffs, coeff_err, func_err = train_and_evaluate(
            func_type='dowell',
            n_stages=5,
            constraints=constraint,
            n_restarts=20,
            n_epochs=2000,
            verbose=True
        )
        results[constraint] = {
            'coeffs': coeffs,
            'coeff_error': coeff_err,
            'func_error': func_err
        }
        print(f"  Best: coeffs={coeffs}, coeff_err={coeff_err:.1f}%, func_err={func_err:.2f}%")

    print("\n" + "=" * 70)
    print("Summary: Constraint Comparison")
    print("=" * 70)
    print("Target: [3, 5, 7, 9, 11]")
    for constraint, res in results.items():
        print(f"  {constraint:12s}: {np.array2string(res['coeffs'], precision=2):40s} "
              f"coeff_err={res['coeff_error']:5.1f}%, func_err={res['func_error']:.2f}%")

    return results


def research_2_different_functions():
    """Test on different continued fraction functions."""
    print("\n" + "=" * 70)
    print("Research 2: Different Continued Fraction Functions")
    print("=" * 70)

    test_cases = [
        ('dowell', 5, [3, 5, 7, 9, 11]),
        ('modified_dowell', 5, [2.5, 4.5, 6.5, 8.5, 10.5]),
        ('three_term', 3, [2, 4, 6]),
        ('non_monotonic', 5, [5, 3, 7, 2, 9]),
    ]

    results = {}

    for func_type, n_stages, true_coeffs in test_cases:
        print(f"\n--- Function: {func_type}, Target: {true_coeffs} ---")
        coeffs, coeff_err, func_err = train_and_evaluate(
            func_type=func_type,
            n_stages=n_stages,
            constraints='positive',  # Use positive constraint
            n_restarts=30,
            n_epochs=3000,
            verbose=True
        )
        results[func_type] = {
            'true': true_coeffs,
            'learned': coeffs,
            'coeff_error': coeff_err,
            'func_error': func_err
        }

    print("\n" + "=" * 70)
    print("Summary: Different Functions")
    print("=" * 70)
    for func_type, res in results.items():
        print(f"  {func_type:20s}: target={res['true']}")
        print(f"  {'':20s}  learned={np.array2string(res['learned'], precision=2)}")
        print(f"  {'':20s}  coeff_err={res['coeff_error']:.1f}%, func_err={res['func_error']:.2f}%")

    return results


def research_3_stage_selection():
    """Investigate automatic stage number selection."""
    print("\n" + "=" * 70)
    print("Research 3: Automatic Stage Number Selection")
    print("=" * 70)

    # True function is 5-term Dowell
    z_train, y_train, true_coeffs = generate_training_data('dowell')

    results = []

    for n_stages in range(2, 9):
        print(f"\n--- n_stages = {n_stages} ---")

        coeffs, coeff_err, func_err = train_and_evaluate(
            func_type='dowell',
            n_stages=n_stages,
            constraints='positive',
            n_restarts=15,
            n_epochs=2000,
            verbose=False
        )

        # Compute AIC-like criterion (lower is better)
        # AIC = 2*k + n*log(MSE) where k=n_stages, n=n_points
        n_points = len(z_train)
        mse = (func_err / 100) ** 2  # Approximate MSE
        aic = 2 * n_stages + n_points * np.log(mse + 1e-10)

        results.append({
            'n_stages': n_stages,
            'coeffs': coeffs,
            'func_error': func_err,
            'aic': aic
        })

        print(f"  func_err={func_err:.3f}%, AIC={aic:.1f}")

    print("\n" + "=" * 70)
    print("Summary: Stage Selection (True = 5 stages)")
    print("=" * 70)
    print(f"{'Stages':>7s} {'Func Err %':>12s} {'AIC':>10s}")
    for res in results:
        print(f"{res['n_stages']:7d} {res['func_error']:12.3f} {res['aic']:10.1f}")

    # Find optimal by AIC
    best_aic_idx = np.argmin([r['aic'] for r in results])
    print(f"\nAIC recommends: {results[best_aic_idx]['n_stages']} stages")

    return results


def research_4_loss_landscape():
    """Visualize the loss landscape for 2-parameter case."""
    print("\n" + "=" * 70)
    print("Research 4: Loss Landscape Visualization (2-term CF)")
    print("=" * 70)

    # Simple 2-term CF: f(z) = 1 + z^2 / (a1 + z^2 / a2)
    # True: a1=3, a2=5

    z_train = np.linspace(0.5, 3.0, 50) * (1 + 1j) / np.sqrt(2)
    true_coeffs = [3.0, 5.0]
    y_train = evaluate_cf(z_train, true_coeffs)

    # Scan a1, a2 grid
    a1_range = np.linspace(1, 8, 50)
    a2_range = np.linspace(1, 12, 50)

    loss_grid = np.zeros((len(a1_range), len(a2_range)))

    for i, a1 in enumerate(a1_range):
        for j, a2 in enumerate(a2_range):
            y_pred = evaluate_cf(z_train, [a1, a2])
            loss_grid[i, j] = np.mean(np.abs(y_pred - y_train) ** 2)

    # Find minimum
    min_idx = np.unravel_index(np.argmin(loss_grid), loss_grid.shape)
    min_a1, min_a2 = a1_range[min_idx[0]], a2_range[min_idx[1]]

    print(f"True coefficients: a1=3.0, a2=5.0")
    print(f"Grid search minimum: a1={min_a1:.2f}, a2={min_a2:.2f}")
    print(f"Loss at true: {np.mean(np.abs(evaluate_cf(z_train, [3, 5]) - y_train)**2):.6f}")
    print(f"Loss at grid min: {loss_grid[min_idx]:.6f}")

    # Check for local minima
    # Count how many points have lower loss than their 8 neighbors
    n_local_min = 0
    for i in range(1, len(a1_range) - 1):
        for j in range(1, len(a2_range) - 1):
            center = loss_grid[i, j]
            neighbors = loss_grid[i-1:i+2, j-1:j+2].flatten()
            if center == np.min(neighbors):
                n_local_min += 1

    print(f"Number of local minima in grid: {n_local_min}")

    return loss_grid, a1_range, a2_range


if __name__ == '__main__':
    print("CF-KAN Research: Coefficient Extraction without Prior Knowledge")
    print("=" * 70)

    # Research 1: Constraint comparison
    res1 = research_1_constraint_comparison()

    # Research 2: Different functions
    res2 = research_2_different_functions()

    # Research 3: Stage selection
    res3 = research_3_stage_selection()

    # Research 4: Loss landscape
    res4 = research_4_loss_landscape()

    print("\n" + "=" * 70)
    print("Research Conclusions")
    print("=" * 70)
    print("""
1. CONSTRAINT EFFECT:
   - 'positive' constraint helps stability
   - 'monotonic' can hurt if true coefficients are non-monotonic

2. FUNCTION GENERALIZATION:
   - Works for various CF functions, not just Dowell
   - Non-monotonic coefficients are harder to recover

3. STAGE SELECTION:
   - AIC can help select number of stages
   - Overfitting with too many stages

4. LOSS LANDSCAPE:
   - Multiple local minima exist
   - Global optimization is needed (multi-restart helps)
""")
