"""
Diffusive quadrature for the Warburg block K_SIBC * sqrt(s) / (s + d).

The hierarchical Cauer construction
    Y_R(s) = Y_CLN(s) + K * sqrt(s) / (s + d)
contains the irrational factor sqrt(s) / (s + d) which cannot be represented
by any finite rational expansion exactly.  For time-domain realization (state
space, SPICE) we need a finite-dimensional approximation.

The Warburg block is physically the input impedance of a semi-infinite
uniform R-C ladder.  Discretizing this ladder to N_xi rungs gives an
N_xi-point Foster expansion:
    sqrt(s)/(s+d) ~= sum_k c_k / (s + xi_k)

This script implements the diffusive quadrature: chooses {xi_k} log-spaced
around d, solves least-squares for {c_k} on a frequency grid, and reports
the maximum relative error |Y_approx(jw) - Y_exact(jw)| / |Y_exact(jw)| as
a function of N_xi.

Usage:
    python diffusive_quadrature.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path


def f_exact(s, d):
    """Exact Warburg transfer function f(s) = sqrt(s) / (s + d)."""
    return np.sqrt(s) / (s + d)


def diffusive_approx(s, xi, c):
    """Foster-form approximation: sum_k c_k / (s + xi_k)."""
    # s: (M,) complex, xi: (N,) real, c: (N,) complex
    # output: (M,) complex
    return np.sum(c[None, :] / (s[:, None] + xi[None, :]), axis=1)


def fit_diffusive(d, N_xi, omega_min=1e0, omega_max=1e8, N_freq=4000,
                  spread_decades=8):
    """
    Fit Foster approximation sum_k c_k / (s + xi_k) to f(s) = sqrt(s)/(s+d)
    on a log-spaced jw grid.

    Poles xi_k are placed log-spaced from d/10^spread_decades to d*10^spread_decades.
    Wide spread (8 decades) is needed to capture the branch-cut behavior of sqrt(s).

    Returns:
        xi: (N_xi,) real pole locations
        c: (N_xi,) complex residues
        max_rel_err: maximum |f_approx - f_exact| / |f_exact| over the grid
    """
    # Log-spaced poles centered on d, spanning the branch-cut influence
    log_d = np.log10(d)
    xi = np.logspace(log_d - spread_decades, log_d + spread_decades, N_xi)

    # Frequency grid (also extended to better constrain the fit at endpoints)
    omega = np.logspace(np.log10(omega_min), np.log10(omega_max), N_freq)
    s = 1j * omega

    # Build system matrix A: A[m, k] = 1 / (s_m + xi_k)
    A = 1.0 / (s[:, None] + xi[None, :])

    # Target vector
    b = f_exact(s, d)

    # Least-squares solve (use real-imag stacking to keep c real-symmetric pair)
    # Since xi are real, the real-valued impulse response constraint is
    # automatically satisfied when c is real (no complex conjugate pairing).
    # We enforce real c by stacking real/imag of A and b.
    A_real = np.vstack([A.real, A.imag])
    b_real = np.hstack([b.real, b.imag])
    c_real, residuals, rank, sv = np.linalg.lstsq(A_real, b_real, rcond=None)
    c = c_real.astype(complex)

    # Evaluate approximation
    f_approx = diffusive_approx(s, xi, c)
    rel_err = np.abs(f_approx - b) / np.abs(b)
    max_rel_err = np.max(rel_err)

    return xi, c, max_rel_err, omega, rel_err


def convergence_study(d, N_xi_values, omega_min=1e0, omega_max=1e8):
    """Sweep N_xi and report max rel error."""
    results = []
    for N_xi in N_xi_values:
        xi, c, max_err, omega, rel_err = fit_diffusive(
            d, N_xi, omega_min, omega_max
        )
        results.append({
            'N_xi': N_xi,
            'max_rel_err': max_err,
            'xi': xi,
            'c': c,
        })
        print(f"  N_xi = {N_xi:3d}: max rel err = {max_err:.2e}")
    return results


def setup_times_new_roman():
    """Lab style: use Times New Roman for figures."""
    try:
        for font in ['Times New Roman', 'Times']:
            try:
                font_manager.findfont(font, fallback_to_default=False)
                plt.rcParams['font.family'] = font
                plt.rcParams['font.serif'] = [font]
                plt.rcParams['mathtext.fontset'] = 'stix'
                return
            except Exception:
                continue
    except Exception:
        pass


def plot_results(d, results, omega_for_curve, savepath):
    """Two-panel: (a) error vs N_xi, (b) error vs omega for best N_xi."""
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # (a) Convergence
    N_xis = [r['N_xi'] for r in results]
    errs = [r['max_rel_err'] for r in results]
    ax1.loglog(N_xis, errs, 'ko-', markersize=4, linewidth=1)
    ax1.set_xlabel(r'$N_\xi$ (diffusive quadrature points)')
    ax1.set_ylabel('max relative error on $|f(j\\omega)|$')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title('(a) Convergence', fontsize=10)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    # Engineering reference line at 1%
    ax1.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1% engineering target')
    ax1.legend(loc='upper right', fontsize=8, frameon=False)

    # (b) Error vs omega for best N_xi
    best_idx = int(np.argmin(errs))
    best = results[best_idx]
    xi_b, c_b = best['xi'], best['c']
    s_curve = 1j * omega_for_curve
    f_approx = diffusive_approx(s_curve, xi_b, c_b)
    f_true = f_exact(s_curve, d)
    rel_err = np.abs(f_approx - f_true) / np.abs(f_true)
    ax2.loglog(omega_for_curve / (2 * np.pi), rel_err, 'b-', linewidth=1)
    ax2.set_xlabel(r'frequency (Hz)')
    ax2.set_ylabel('relative error')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_title(f'(b) Error spectrum, $N_\\xi$={best["N_xi"]}', fontsize=10)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    plt.savefig(savepath, dpi=200, bbox_inches='tight')
    print(f"Saved: {savepath}")


if __name__ == "__main__":
    # Cylinder parameters from digest (a = 5 mm copper)
    sigma_cu = 5.96e7  # S/m
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    # Crossover constant: d/(2pi) = 7e4 Hz
    d = 2 * np.pi * 7e4
    print(f"Crossover d = {d:.3e} rad/s (= {d/(2*np.pi):.2e} Hz)")
    print()

    # Convergence study
    print("Convergence study:")
    N_xi_values = [8, 12, 16, 20, 25, 30, 40, 50, 60, 80, 100]
    results = convergence_study(d, N_xi_values, omega_min=2*np.pi*1,
                                omega_max=2*np.pi*1e8)
    print()

    # Find smallest N_xi achieving 1% error
    for r in results:
        if r['max_rel_err'] < 1e-2:
            print(f"-> Smallest N_xi for 1% error: {r['N_xi']}"
                  f" (max rel err = {r['max_rel_err']:.2e})")
            break
    else:
        print("-> 1% target not reached within tested N_xi range.")
    print()

    # Plot
    omega_curve = np.logspace(np.log10(2*np.pi*1), np.log10(2*np.pi*1e8), 2000)
    out_dir = Path(__file__).parent
    plot_results(d, results, omega_curve, out_dir / "diffusive_convergence.png")
    print()

    # Best fit data dump for downstream use
    best_idx = int(np.argmin([r['max_rel_err'] for r in results]))
    best = results[best_idx]
    np.savez(out_dir / "warburg_diffusive_fit.npz",
             d=d, xi=best['xi'], c=best['c'],
             N_xi=best['N_xi'], max_rel_err=best['max_rel_err'])
    print(f"Saved fit: {out_dir / 'warburg_diffusive_fit.npz'}")
