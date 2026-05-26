"""
Time-domain step response: hierarchical Cauer vs Bessel exact (cylinder).

Verifies that the hierarchical Cauer realization
    Y_R(s) = Y_CLN(s) + K_SIBC * sqrt(s) / (s + d)
when realized as finite LTI via diffusive quadrature, produces a step
response y(t) that tracks the exact cylinder admittance time-domain
solution.

Exact reference: long Cu cylinder of radius a, with
    Y_cyl(s) = pi a^2 sigma * 2 I_1(gamma a) / [gamma a I_0(gamma a)]
    gamma = sqrt(s mu sigma)
Step response y_cyl(t) is obtained by numerical inverse Laplace transform
(de Hoog / Talbot).  We use scipy.special for Bessel and a basic Talbot
method.

Hierarchical Cauer:
    Y_CLN: rank-10 R-terminated CLN (rational, easy state space)
    Warburg: diffusive Foster expansion from diffusive_quadrature.py
    Combined LTI: (10 + N_xi) states, all first-order ODEs

Usage:
    python time_domain_verify_cylinder.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy import special


# ----------------------------------------------------------------------
# Exact cylinder admittance Y_cyl(s)
# ----------------------------------------------------------------------
def Y_cyl_exact(s, a, sigma, mu):
    """Exact 2D cylinder admittance via Bessel I0/I1.  Uses asymptotic ratio
    for |gamma a| > 50 to avoid overflow in the inverse-Laplace contour."""
    s = np.atleast_1d(s) if not np.isscalar(s) else np.array([s])
    gamma_a = np.sqrt(s * mu * sigma) * a
    abs_ga = np.abs(gamma_a)
    ratio = np.empty_like(gamma_a, dtype=complex)
    large = abs_ga > 50
    # Asymptotic: I_1(z) / I_0(z) ~ 1 - 1/(2z) - 1/(8z^2) + 1/(8z^3) + ...
    z_large = gamma_a[large]
    ratio[large] = (1 - 1.0 / (2 * z_large)
                    - 1.0 / (8 * z_large**2)
                    + 1.0 / (8 * z_large**3))
    # Direct calc for moderate |gamma a|
    if np.any(~large):
        z_small = gamma_a[~large]
        i0 = special.iv(0, z_small)
        i1 = special.iv(1, z_small)
        ratio[~large] = i1 / i0
    result = np.pi * a**2 * sigma * 2 * ratio / gamma_a
    return result if result.size > 1 else result[0]


# ----------------------------------------------------------------------
# Rank-N CLN (R-terminated) from Foster spectrum (Bessel zeros)
# ----------------------------------------------------------------------
def foster_spectrum_cyl(N_mode, a, sigma, mu):
    """
    Foster decomposition Y_cyl(s) = sum_k a_k / (s + xi_k) (truncated).
    For cylinder, poles are at zeros of J_0:
        xi_k = j_{0,k}^2 / (a^2 mu sigma)
    and residues are
        a_k = 4 / (mu * a^2 * j_{0,k}^2 * ...).
    We compute by direct Foster fit using known formula.

    Actually, the closed-form Foster of Y_cyl:
        Y_cyl(s) = sum_k 4 pi sigma / [j_{0,k}^2 (1 + s * a^2 mu sigma / j_{0,k}^2)]
                 = sum_k 4 pi / (mu * j_{0,k}^2) / (s + xi_k) * 1
        with xi_k = j_{0,k}^2 / (a^2 mu sigma)
    """
    # Bessel function zeros j_{0,k}
    j0_zeros = special.jn_zeros(0, N_mode)
    xi = j0_zeros**2 / (a**2 * mu * sigma)
    # Residues: from cylinder Y(s) Foster, the amplitude is 4/(mu * j_{0,k}^2)
    # Per zero, this gives Y(s) = sum 4 pi a^2 sigma / j_{0,k}^2 * 1/(1 + s/xi_k)
    # in residue form: a_k = (4 pi a^2 sigma / j_{0,k}^2) * xi_k
    #                    = 4 pi / mu  (interesting -- all residues identical)
    residue = 4 * np.pi / mu  # residue per pole, in units (A/m) / (V/m)
    a_k = residue * np.ones_like(xi)
    return xi, a_k


def Y_CLN_RTerm(s, N, a, sigma, mu):
    """
    R-terminated CLN of N modes.  In Foster form, this is just the first N
    poles of the cylinder Foster expansion.  This is the "rational ROM" that
    saturates at the resistive plateau above the N-th pole.
    """
    xi, a_k = foster_spectrum_cyl(N, a, sigma, mu)
    return np.sum(a_k[None, :] / (s[:, None] + xi[None, :]), axis=1)


# ----------------------------------------------------------------------
# Hierarchical Cauer: Y_R = Y_CLN + K_SIBC * sqrt(s) / (s + d)
# Approximated via diffusive Foster for the Warburg block
# ----------------------------------------------------------------------
def Y_R_hierarchical(s, N_cln, a, sigma, mu, d, K_SIBC, xi_w, c_w):
    """
    Hierarchical Cauer: Y_CLN (N_cln modes) + Warburg diffusive Foster
    (N_xi rungs) approximating K_SIBC * sqrt(s)/(s+d).
    """
    Y_cln = Y_CLN_RTerm(s, N_cln, a, sigma, mu)
    # Diffusive Foster of K_SIBC * sqrt(s) / (s+d)
    Y_warburg = K_SIBC * np.sum(c_w[None, :] / (s[:, None] + xi_w[None, :]),
                                axis=1)
    return Y_cln + Y_warburg


# ----------------------------------------------------------------------
# Numerical inverse Laplace transform (Talbot method)
# ----------------------------------------------------------------------
def talbot_inverse(F_func, t, M=64):
    """
    Inverse Laplace transform via Talbot's method.
    F_func: callable F(s) for complex s.
    t: scalar time > 0.
    M: number of quadrature nodes (even).

    Returns f(t) such that F(s) = L{f}(s).
    """
    # Talbot parameters
    k = np.arange(M)
    theta = -np.pi + (k + 0.5) * 2 * np.pi / M

    # Contour: s = (M/(5t)) * (theta * cot(theta) - 1 + i*nu*theta)
    # Standard Talbot parameter set (Abate-Whitt):
    sigma_t = 0.6122
    mu_t = 0.5017
    nu_t = 0.2645

    # Talbot contour points
    cot_theta = 1.0 / np.tan(theta)
    sigma = (M / t) * (sigma_t * theta * cot_theta - mu_t + 1j * nu_t * theta)
    # Derivative ds/dtheta
    d_sigma = (M / t) * (
        sigma_t * (cot_theta - theta / np.sin(theta)**2)
        + 1j * nu_t
    )

    # Evaluate F on contour
    F_vals = F_func(sigma)

    # Trapezoidal sum on contour
    integrand = np.exp(sigma * t) * F_vals * d_sigma
    f_t = (1.0 / (2j * np.pi)) * (2 * np.pi / M) * np.sum(integrand)
    return f_t.real


def step_response(Y_func, t_array, M=64):
    """Step response: y(t) = L^-1[ Y(s) / s ](t)."""
    y = np.zeros_like(t_array)
    for i, t in enumerate(t_array):
        if t <= 0:
            y[i] = 0
            continue
        # H(s) = Y(s) / s, so step input gives integral
        F_step = lambda s: Y_func(s) / s
        y[i] = talbot_inverse(F_step, t, M=M)
    return y


# ----------------------------------------------------------------------
# Lab style helper
# ----------------------------------------------------------------------
def setup_times_new_roman():
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


# ----------------------------------------------------------------------
# Main verification
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Cylinder parameters
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    K_SIBC = 2 * np.pi * a * np.sqrt(sigma_cu / mu_0)
    d = 2 * np.pi * 7e4
    print(f"K_SIBC^cyl = {K_SIBC:.3e}")
    print(f"d = {d:.3e} rad/s")
    print()

    # Load Warburg diffusive fit
    fit_path = Path(__file__).parent / "warburg_diffusive_fit.npz"
    if not fit_path.exists():
        print(f"ERROR: run diffusive_quadrature.py first to produce {fit_path}")
        exit(1)
    data = np.load(fit_path)
    xi_w = data['xi']
    c_w = data['c']
    print(f"Loaded Warburg fit: N_xi = {len(xi_w)}, "
          f"max rel err = {data['max_rel_err']:.2e}")
    print()

    # Verify frequency-domain agreement first (sanity check)
    # Use N_cln = 30 so the bulk Foster covers up to ~730 kHz (well above
    # wall band d/(2pi) = 70 kHz); the Warburg block handles the SIBC tail.
    N_cln = 30
    omega_test = np.logspace(0, 8, 200) * 2 * np.pi
    s_test = 1j * omega_test
    Y_exact = Y_cyl_exact(s_test, a, sigma_cu, mu_0)
    Y_approx = Y_R_hierarchical(s_test, N_cln, a, sigma_cu, mu_0,
                                d, K_SIBC, xi_w, c_w)
    rel_err = np.abs(Y_approx - Y_exact) / np.abs(Y_exact)
    print(f"Frequency-domain: max rel err |Y|/|Y_exact| = {np.max(rel_err):.2e}")
    print(f"  at f = {omega_test[np.argmax(rel_err)]/(2*np.pi):.2e} Hz")
    print(f"  N_CLN = {N_cln} Foster modes, N_xi = {len(xi_w)} Warburg points")
    print()

    # Step response over multiple time scales
    # Characteristic times:
    #   t_d = 1/d  ~ 2.3e-6 s    (wall-band crossover)
    #   t_skin (at f=10MHz) ~ 1e-7 s
    #   t_DC ~ 1e-3 s
    t_array = np.logspace(-8, -3, 60)
    print("Computing step response (exact)...")
    Y_exact_func = lambda s: Y_cyl_exact(s, a, sigma_cu, mu_0)
    y_exact = step_response(Y_exact_func, t_array, M=64)

    print(f"Computing step response (hierarchical Cauer, N_CLN={N_cln})...")
    Y_approx_func = lambda s: Y_R_hierarchical(
        s, N_cln, a, sigma_cu, mu_0, d, K_SIBC, xi_w, c_w
    )
    y_approx = step_response(Y_approx_func, t_array, M=64)

    # Error metrics
    rel_err_t = np.abs(y_approx - y_exact) / np.maximum(np.abs(y_exact), 1e-30)
    max_t_err = np.nanmax(rel_err_t)
    print()
    print(f"Time-domain: max rel err on step response = {max_t_err:.2e}")
    print()

    # Plot
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    ax1.loglog(t_array, np.abs(y_exact), 'k-', linewidth=1.2,
               label='exact Bessel')
    ax1.loglog(t_array, np.abs(y_approx), 'r--', linewidth=1.0,
               label='hierarchical Cauer (Talbot)')
    ax1.axvline(1.0 / d, color='gray', linestyle=':', linewidth=0.6,
                label='$t = 1/d$')
    ax1.set_xlabel('time $t$ (s)')
    ax1.set_ylabel('step response $|y(t)|$')
    ax1.legend(fontsize=8, frameon=False, loc='lower left')
    ax1.grid(True, which='both', alpha=0.3)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    ax2.loglog(t_array, rel_err_t, 'b-', linewidth=1.0)
    ax2.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1% engineering')
    ax2.set_xlabel('time $t$ (s)')
    ax2.set_ylabel('relative error on $y(t)$')
    ax2.legend(fontsize=8, frameon=False, loc='upper left')
    ax2.grid(True, which='both', alpha=0.3)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "time_domain_step_response.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'time_domain_step_response.png'}")

    # Save data for downstream
    np.savez(out_dir / "step_response_data.npz",
             t=t_array, y_exact=y_exact, y_approx=y_approx,
             rel_err=rel_err_t, max_rel_err=max_t_err,
             K_SIBC=K_SIBC, d=d, N_cln=10, N_xi=len(xi_w))
    print(f"Saved: {out_dir / 'step_response_data.npz'}")
