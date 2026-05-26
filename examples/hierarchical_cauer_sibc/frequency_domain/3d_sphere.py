"""
3D extension verification: conducting sphere.

The 2D cylinder digest result generalizes to 3D via surface-area scaling
K_SIBC = S * sqrt(sigma/mu).  For a sphere of radius a:
    S = 4 pi a^2     (surface area)
    V = (4/3) pi a^3 (volume)
    K_SIBC^sph = 4 pi a^2 sqrt(sigma/mu)
    Y(0) = V sigma

Exact admittance via modified spherical Bessel:
    Y_sph(s) = (4/3) pi a^3 sigma * (3 / (gamma a)) * (coth(gamma a) - 1/(gamma a))

Mittag-Leffler (Foster) expansion:
    Y_sph(s) = sum_{k=1}^inf (8 pi a / mu) / (s + xi_k)
    xi_k = (k pi)^2 / (a^2 mu sigma)

Limits:
    s -> 0:    Y -> V sigma                      (DC)
    s -> inf:  Y -> K_SIBC^sph / sqrt(s)         (SIBC tail)

Verification:
    (i) Show K_SIBC^sph = 4 pi a^2 sqrt(sigma/mu) matches exact at high f.
    (ii) Hierarchical Cauer with PoU: tracks Y_sph_exact across 8 decades
         with wall-band error comparable to cylinder digest result.

Usage:
    python 3d_sphere.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path


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


def coth_complex(z):
    """Numerically stable coth(z) for complex z."""
    z = np.asarray(z, dtype=complex)
    result = np.empty_like(z, dtype=complex)
    big = np.real(z) > 30
    sml = np.real(z) < -30
    mid = ~(big | sml)
    result[big] = 1.0
    result[sml] = -1.0
    result[mid] = 1.0 + 2.0 / (np.exp(2 * z[mid]) - 1)
    return result


def Y_sph_exact(s, a, sigma, mu):
    """Exact 3D sphere admittance."""
    gamma_a = np.sqrt(s * mu * sigma) * a
    coth_ga = coth_complex(gamma_a)
    prefactor = (4.0 / 3.0) * np.pi * a**3 * sigma
    bessel_ratio = (3.0 / gamma_a) * (coth_ga - 1.0 / gamma_a)
    return prefactor * bessel_ratio


def sphere_pade(N_modes_for_moments, a, sigma, mu, m=5, n=6):
    """
    Construct Pade[m/n] approximation of Y_sph(s) at s=0 from numerical
    Taylor moments computed via Mittag-Leffler (Foster) expansion.

    Y_sph(s) = sum_k r / (s + xi_k), r = 8 pi a / mu, xi_k = (k pi)^2/(a^2 mu sigma)
    Taylor: c_n = r (-1)^n sum_k xi_k^{-(n+1)}

    Returns (num_coefs, den_coefs) of the Pade approximation, such that
    Y_pade(s) = sum_i num_coefs[i] s^i / sum_j den_coefs[j] s^j
    """
    r = 8 * np.pi * a / mu
    k = np.arange(1, N_modes_for_moments + 1)
    xi = (k * np.pi)**2 / (a**2 * mu * sigma)

    # Compute Taylor moments c_n for n = 0, 1, ..., m+n
    n_moments = m + n + 1
    c = np.zeros(n_moments)
    for nn in range(n_moments):
        c[nn] = r * (-1)**nn * np.sum(1.0 / xi**(nn + 1))

    # Pade[m/n]: solve linear system for denominator coefficients then numerator
    # Denominator: b_0 = 1, solve for b_1..b_n from
    #   sum_{j=0}^n b_j c_{i-j} = 0 for i = m+1, ..., m+n
    # where c_k = 0 for k < 0
    from scipy.interpolate import pade
    num_poly, den_poly = pade(c, n, m)
    return num_poly, den_poly


def Y_Pade_sphere(s, num_poly, den_poly):
    """Evaluate Pade approximation at complex s."""
    return num_poly(s) / den_poly(s)


def Y_R_PoU(s, num_poly, den_poly, d, K_SIBC, n_pou=3):
    """
    PoU-corrected hierarchical Cauer:
        Y_R^PoU(s) = w_d(s) Y_Pade(s) + w_s(s) K_SIBC / sqrt(s)
    with w_d = d^n / (s^n + d^n), w_s = s^n / (s^n + d^n).
    """
    w_d = d**n_pou / (s**n_pou + d**n_pou)
    w_s = s**n_pou / (s**n_pou + d**n_pou)
    return w_d * Y_Pade_sphere(s, num_poly, den_poly) + w_s * K_SIBC / np.sqrt(s)


def Y_R_additive(s, num_poly, den_poly, d, K_SIBC):
    """Bare additive Schur form: Y_R = Y_Pade + K_SIBC sqrt(s)/(s+d)."""
    return Y_Pade_sphere(s, num_poly, den_poly) + K_SIBC * np.sqrt(s) / (s + d)


if __name__ == "__main__":
    # Sphere parameters
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    K_SIBC = 4 * np.pi * a**2 * np.sqrt(sigma_cu / mu_0)
    V_sphere = (4.0 / 3.0) * np.pi * a**3
    Y_DC = V_sphere * sigma_cu

    print(f"3D Sphere parameters (Cu, a = {a*1e3:.1f} mm):")
    print(f"  V = {V_sphere:.3e} m^3, S = {4*np.pi*a**2:.3e} m^2")
    print(f"  Y_DC = V sigma = {Y_DC:.4e} S")
    print(f"  K_SIBC^sph = 4 pi a^2 sqrt(sigma/mu) = {K_SIBC:.4e}")
    print()

    # Build Pade[5/6] from Mittag-Leffler moments
    print("Building Pade[5/6] approximation of Y_sph from 1000 Foster modes...")
    num_poly, den_poly = sphere_pade(1000, a, sigma_cu, mu_0, m=5, n=6)
    print(f"  Pade num degree = {num_poly.order}, den degree = {den_poly.order}")
    print(f"  Y_Pade(0) = {Y_Pade_sphere(np.array([0+0j]), num_poly, den_poly)[0].real:.4e}")
    print(f"  expected Y_DC = {Y_DC:.4e}")
    print()

    # Frequency sweep
    omega = np.logspace(0, 8, 200) * 2 * np.pi
    s = 1j * omega
    Y_exact = Y_sph_exact(s, a, sigma_cu, mu_0)

    # Crossover d: optimize on PoU
    from scipy.optimize import minimize_scalar
    def err_of_d(log_d):
        d_try = 10**log_d
        Y_pou = Y_R_PoU(s, num_poly, den_poly, d_try, K_SIBC)
        return np.max(np.abs(Y_pou - Y_exact) / np.abs(Y_exact))

    res = minimize_scalar(err_of_d, bounds=(3, 8), method='bounded')
    d_opt = 10**res.x
    err_opt = res.fun
    print(f"Optimal d for PoU: d = {d_opt:.3e} rad/s "
          f"({d_opt/(2*np.pi):.2e} Hz), max err = {err_opt:.2e}")
    print()

    # Final evaluation
    Y_PoU = Y_R_PoU(s, num_poly, den_poly, d_opt, K_SIBC)
    Y_add = Y_R_additive(s, num_poly, den_poly, d_opt, K_SIBC)
    rel_err_PoU = np.abs(Y_PoU - Y_exact) / np.abs(Y_exact)
    rel_err_add = np.abs(Y_add - Y_exact) / np.abs(Y_exact)

    print(f"Wall-band error comparison (d = {d_opt/(2*np.pi):.2e} Hz):")
    print(f"  PoU corrected     : max err = {np.max(rel_err_PoU):.2e}")
    print(f"  Additive Schur    : max err = {np.max(rel_err_add):.2e}")
    print()

    # Plot
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 3.0))

    f_hz = omega / (2 * np.pi)
    ax1.loglog(f_hz, np.abs(Y_exact), 'k-', linewidth=1.4,
               label='exact (spherical Bessel)')
    ax1.loglog(f_hz, np.abs(Y_PoU), 'r--', linewidth=1.0,
               label='hierarchical Cauer + PoU')
    ax1.loglog(f_hz, K_SIBC / np.sqrt(omega), 'b:', linewidth=0.8, alpha=0.6,
               label=r'$K_{\mathrm{SIBC}}/\sqrt{\omega}$ asymptote')
    ax1.axhline(Y_DC, color='gray', linestyle=':', linewidth=0.6, alpha=0.6,
                label=r'$Y(0) = V\sigma$')
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'$|Y(j\omega)|$ (S)')
    ax1.legend(fontsize=8, frameon=False, loc='upper right')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title(f'(a) Cu sphere $a$={a*1e3:.1f}mm', fontsize=10)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    ax2.loglog(f_hz, rel_err_add, 'b-', linewidth=1.0,
               label=f'additive Schur ({np.max(rel_err_add)*100:.1f}\\% peak)')
    ax2.loglog(f_hz, rel_err_PoU, 'r--', linewidth=1.0,
               label=f'PoU corrected ({np.max(rel_err_PoU)*100:.2f}\\% peak)')
    ax2.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1\\% engineering')
    ax2.set_xlabel('frequency (Hz)')
    ax2.set_ylabel('relative error')
    ax2.legend(fontsize=8, frameon=False, loc='lower left')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_title('(b) Error', fontsize=10)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "3d_sphere_verification.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / '3d_sphere_verification.png'}")

    np.savez(out_dir / "3d_sphere_data.npz",
             omega=omega, Y_exact=Y_exact, Y_PoU=Y_PoU, Y_additive=Y_add,
             rel_err_PoU=rel_err_PoU, rel_err_additive=rel_err_add,
             K_SIBC=K_SIBC, d=d_opt, Y_DC=Y_DC)
    print(f"Saved: {out_dir / '3d_sphere_data.npz'}")
