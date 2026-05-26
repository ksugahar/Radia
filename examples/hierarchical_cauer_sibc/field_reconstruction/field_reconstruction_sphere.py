"""
Sphere step response field reconstruction.

For a conducting sphere of radius a in external uniform field H_ext(t) = H_0 u(t),
the diffusion equation in spherical coordinates (axisymmetric, l=0 angular mode)
gives the exact modal solution:

    H(r, t) / H_0 = 1 - sum_{k=1}^inf c_k * (sin(k pi r/a) / (k pi r/a)) * exp(-xi_k t)

with eigenvalues xi_k = (k pi)^2 / (a^2 mu sigma) and modal coefficients
    c_k = 2 (-1)^{k+1}

(derived from spherical Bessel orthogonality:
    int_0^a [sin(k pi r/a)/(k pi r/a)] r^2 dr = (-1)^{k+1} a^4 / (k^2 pi^2)
    and the eigenfunction norm a^3/(2 k^2 pi^2).)

This is the sphere analog of the cylinder Bessel J_0 modal expansion.  State
variables z_k(t) = c_k (1 - exp(-xi_k t)) ARE the spherical mode amplitudes.

Usage:
    python field_reconstruction_sphere.py
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


def H_sph_step(r, t, N_modes, a, sigma, mu):
    """
    H(r, t) / H_0 for sphere step input.
    r: (M,) array with r in [0, a], t: scalar.
    Special handling at r=0: sin(x)/x -> 1.
    """
    k = np.arange(1, N_modes + 1)
    xi = (k * np.pi)**2 / (a**2 * mu * sigma)
    c = 2.0 * (-1)**(k + 1)
    decay = np.exp(-xi * t)
    # spatial mode: sinc(kr/a) = sin(k pi r/a) / (k pi r/a)
    # at r=0: limit is 1
    eps = 1e-15
    r_safe = np.where(r > eps, r, eps)
    arg = k[None, :] * np.pi * r_safe[:, None] / a
    sinc_arr = np.sin(arg) / arg
    # At r=0: replace with 1
    if r[0] < eps:
        sinc_arr[0, :] = 1.0
    H = 1.0 - np.sum(c[None, :] * sinc_arr * decay[None, :], axis=1)
    return H


def J_sph_step(r, t, N_modes, a, sigma, mu):
    """
    Induced current density magnitude |J|(r, t) / H_0 for sphere step input.
    J = -dH/dr  (azimuthal current density, taking H = H_z axially)
    """
    k = np.arange(1, N_modes + 1)
    xi = (k * np.pi)**2 / (a**2 * mu * sigma)
    c = 2.0 * (-1)**(k + 1)
    decay = np.exp(-xi * t)
    eps = 1e-15
    r_safe = np.where(r > eps, r, eps)
    arg = k[None, :] * np.pi * r_safe[:, None] / a
    # d/dr [sin(x)/x] where x = k pi r/a:
    # = (k pi/a) * [cos(x)/x - sin(x)/x^2]
    deriv = (k[None, :] * np.pi / a) * (np.cos(arg) / arg
                                         - np.sin(arg) / arg**2)
    if r[0] < eps:
        deriv[0, :] = 0  # d/dr at r=0 is 0 by symmetry
    # H = 1 - sum c_k * sinc * decay  =>  dH/dr = -sum c_k * deriv * decay
    # J = -dH/dr = sum c_k * deriv * decay
    J = np.sum(c[None, :] * deriv * decay[None, :], axis=1)
    return J


if __name__ == "__main__":
    # Cu sphere
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    K_SIBC = 4 * np.pi * a**2 * np.sqrt(sigma_cu / mu_0)
    Y_DC = (4.0/3.0) * np.pi * a**3 * sigma_cu

    print(f"Sphere a = {a*1e3:.1f} mm Cu:")
    print(f"  K_SIBC^sph = {K_SIBC:.4e}")
    print(f"  Y_DC = V sigma = {Y_DC:.4e} S")
    print()

    # Wall-band time: t_d such that delta = sqrt(2/(d mu sigma)) ~ a
    #   d ~ 2/(a^2 mu sigma) for skin penetration through whole sphere
    d_sphere = 2.0 / (a**2 * mu_0 * sigma_cu)
    t_d = 1.0 / d_sphere
    print(f"Estimated wall-band:")
    print(f"  d = 2/(a^2 mu sigma) = {d_sphere:.3e} rad/s ({d_sphere/(2*np.pi):.2e} Hz)")
    print(f"  t_d = 1/d = {t_d:.3e} s")
    print()

    # Radial array
    r = np.linspace(0, a, 200)

    # Snapshots
    N_exact = 200
    N_reduced = 30
    t_snap = np.array([0.01 * t_d, 0.1 * t_d, 1.0 * t_d,
                       10.0 * t_d, 100.0 * t_d])
    t_labels = [r'$t=0.01/d$', r'$t=0.1/d$', r'$t=1/d$',
                r'$t=10/d$', r'$t=100/d$']

    H_exact = np.zeros((len(t_snap), len(r)))
    H_reduced = np.zeros((len(t_snap), len(r)))
    J_exact = np.zeros((len(t_snap), len(r)))
    J_reduced = np.zeros((len(t_snap), len(r)))
    for i, t in enumerate(t_snap):
        H_exact[i] = H_sph_step(r, t, N_exact, a, sigma_cu, mu_0)
        H_reduced[i] = H_sph_step(r, t, N_reduced, a, sigma_cu, mu_0)
        J_exact[i] = J_sph_step(r, t, N_exact, a, sigma_cu, mu_0)
        J_reduced[i] = J_sph_step(r, t, N_reduced, a, sigma_cu, mu_0)

    # Plot
    setup_times_new_roman()
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7.5, 6.0))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(t_snap)))

    for i, (lab, col) in enumerate(zip(t_labels, colors)):
        ax1.plot(r * 1e3, H_exact[i], '-', color=col, linewidth=1.2,
                 label=lab)
        ax1.plot(r * 1e3, H_reduced[i], '--', color=col, linewidth=0.7,
                 alpha=0.8)
    ax1.set_xlabel('radial position $r$ (mm)')
    ax1.set_ylabel(r'$H_z(r, t) / H_0$')
    ax1.legend(fontsize=7, frameon=False, loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'(a) Sphere $H_z$ (solid $N$={N_exact}, dashed $N$={N_reduced})',
                  fontsize=9)
    ax1.set_ylim(-0.1, 1.1)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    for i, (lab, col) in enumerate(zip(t_labels, colors)):
        ax2.plot(r * 1e3, J_exact[i] * a, '-', color=col, linewidth=1.2,
                 label=lab)
        ax2.plot(r * 1e3, J_reduced[i] * a, '--', color=col, linewidth=0.7,
                 alpha=0.8)
    ax2.set_xlabel('radial position $r$ (mm)')
    ax2.set_ylabel(r'$|J|(r, t) \cdot a / H_0$')
    ax2.legend(fontsize=7, frameon=False, loc='upper left')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('(b) Induced current $|J|$', fontsize=9)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    # Error in H
    H_err = np.abs(H_reduced - H_exact)
    for i, (lab, col) in enumerate(zip(t_labels, colors)):
        ax3.semilogy(r * 1e3, H_err[i] + 1e-12, color=col, linewidth=1.0,
                     label=lab)
    ax3.set_xlabel('radial position $r$ (mm)')
    ax3.set_ylabel(r'$|H^{\mathrm{red}} - H^{\mathrm{ex}}|$')
    ax3.legend(fontsize=7, frameon=False, loc='lower left')
    ax3.grid(True, which='both', alpha=0.3)
    ax3.set_title(f'(c) Truncation error in $H$ ($N$={N_reduced})', fontsize=9)
    for spine in ax3.spines.values():
        spine.set_linewidth(0.8)
    ax3.tick_params(direction='in', which='both')

    # Error in J
    J_err = np.abs(J_reduced - J_exact)
    J_scale = np.max(np.abs(J_exact))
    for i, (lab, col) in enumerate(zip(t_labels, colors)):
        ax4.semilogy(r * 1e3, J_err[i] / J_scale + 1e-12, color=col,
                     linewidth=1.0, label=lab)
    ax4.set_xlabel('radial position $r$ (mm)')
    ax4.set_ylabel(r'$|J^{\mathrm{red}} - J^{\mathrm{ex}}|/|J|_{\max}$')
    ax4.legend(fontsize=7, frameon=False, loc='lower left')
    ax4.grid(True, which='both', alpha=0.3)
    ax4.set_title(f'(d) Truncation error in $J$', fontsize=9)
    for spine in ax4.spines.values():
        spine.set_linewidth(0.8)
    ax4.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "field_reconstruction_sphere.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'field_reconstruction_sphere.png'}")
    print()

    for i, t in enumerate(t_snap):
        H_max_err = np.max(H_err[i])
        J_max_err_norm = np.max(J_err[i] / J_scale)
        print(f"  t = {t:.2e} s ({t/t_d:.3g}/d): H err = {H_max_err:.2e}, "
              f"J err (norm) = {J_max_err_norm:.2e}")

    np.savez(out_dir / "field_reconstruction_sphere_data.npz",
             r=r, t_snap=t_snap, H_exact=H_exact, H_reduced=H_reduced,
             J_exact=J_exact, J_reduced=J_reduced,
             a=a, d=d_sphere, N_exact=N_exact, N_reduced=N_reduced)
    print(f"Saved: {out_dir / 'field_reconstruction_sphere_data.npz'}")
