"""
Field reconstruction from modal state variables: cylinder step response.

Demonstrates that the hierarchical Cauer / Foster state variables
{z_k(t)} also encode the spatial field distribution: H_z(r, t) and
J_phi(r, t) inside the conductor are linear combinations of the modal
shapes (Bessel J_0, J_1) weighted by the time-domain state amplitudes.

For a long Cu cylinder of radius a in external axial field H_ext(t) = H_0 u(t)
(step), the diffusion equation gives the exact modal solution:

    H_z(r, t) / H_0 = 1 - sum_{k=1}^inf A_k J_0(j_{0,k} r/a) exp(-xi_k t)

with
    A_k = 2 / (j_{0,k} J_1(j_{0,k}))
    xi_k = j_{0,k}^2 / (a^2 mu sigma)

and induced azimuthal current density:
    J_phi(r, t) / H_0 = -d(H_z/H_0)/dr
                     = sum_{k=1}^inf B_k J_1(j_{0,k} r/a) exp(-xi_k t)
with B_k = 2 / (a J_1(j_{0,k})).

State-space interpretation: each "branch current" i_k(t) ~ (1 - exp(-xi_k t))
in the Foster network IS the modal amplitude, multiplied by the spatial
shape J_0 or J_1.  Field reconstruction is a post-process linear combination
of state variables.

This script compares:
  (i)  Exact: Foster sum with N = 200 modes (effectively infinite)
  (ii) Reduced: Foster sum with N = 30 modes
shows H_z and J_phi snapshots at several characteristic times and
quantifies the Gibbs ringing of the truncation.

Usage:
    python field_reconstruction_cylinder.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy import special


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


def cylinder_modal_data(N_modes, a, sigma, mu):
    """
    Returns Bessel zeros j_{0,k}, eigenvalues xi_k, modal amplitude A_k for
    step response, modal coefficient B_k for J_phi.
    """
    j0_zeros = special.jn_zeros(0, N_modes)
    j1_vals = special.j1(j0_zeros)  # = J_1(j_{0,k})
    xi = j0_zeros**2 / (a**2 * mu * sigma)
    A = 2.0 / (j0_zeros * j1_vals)
    B = 2.0 / (a * j1_vals)
    return j0_zeros, j1_vals, xi, A, B


def H_z_step(r, t, N_modes, a, sigma, mu):
    """
    H_z(r, t) / H_0 for step input H_ext = H_0 u(t).
    r: (M,) array, t: scalar, returns (M,) array.
    """
    j0_z, j1_v, xi, A, B = cylinder_modal_data(N_modes, a, sigma, mu)
    # Bessel J_0(j_{0,k} r/a) for each r and each mode k
    arg = j0_z[None, :] * r[:, None] / a
    J0_arr = special.j0(arg)
    decay = np.exp(-xi * t)
    H = 1.0 - np.sum(A[None, :] * J0_arr * decay[None, :], axis=1)
    return H


def J_phi_step(r, t, N_modes, a, sigma, mu):
    """
    J_phi(r, t) / H_0 for step input.  Returns (M,) array.
    """
    j0_z, j1_v, xi, A, B = cylinder_modal_data(N_modes, a, sigma, mu)
    arg = j0_z[None, :] * r[:, None] / a
    J1_arr = special.j1(arg)
    decay = np.exp(-xi * t)
    J = np.sum(B[None, :] * J1_arr * decay[None, :], axis=1)
    return J


if __name__ == "__main__":
    # Cu cylinder
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    d = 2 * np.pi * 7e4  # wall-band crossover (digest)
    t_d = 1.0 / d

    print(f"Cylinder a = {a*1e3:.1f} mm Cu:")
    print(f"  wall-band time t_d = 1/d = {t_d:.2e} s")
    print()

    # Mode count comparison
    N_exact = 200  # effectively converged
    N_reduced = 30  # reduced model
    print(f"Modes: exact = {N_exact}, reduced = {N_reduced}")
    print()

    # Radial array
    r = np.linspace(0, a, 200)  # avoid r=0 where derivatives may diverge

    # Characteristic times for snapshots
    t_snap = np.array([0.01 * t_d, 0.1 * t_d, 1.0 * t_d,
                       10.0 * t_d, 100.0 * t_d])
    t_labels = [r'$t=0.01 / d$', r'$t=0.1 / d$', r'$t=1 / d$',
                r'$t=10 / d$', r'$t=100 / d$']

    # Compute H_z and J_phi for each snapshot, both exact and reduced
    H_exact = np.zeros((len(t_snap), len(r)))
    H_reduced = np.zeros((len(t_snap), len(r)))
    J_exact = np.zeros((len(t_snap), len(r)))
    J_reduced = np.zeros((len(t_snap), len(r)))
    for i, t in enumerate(t_snap):
        H_exact[i] = H_z_step(r, t, N_exact, a, sigma_cu, mu_0)
        H_reduced[i] = H_z_step(r, t, N_reduced, a, sigma_cu, mu_0)
        J_exact[i] = J_phi_step(r, t, N_exact, a, sigma_cu, mu_0)
        J_reduced[i] = J_phi_step(r, t, N_reduced, a, sigma_cu, mu_0)

    # Plot snapshots
    setup_times_new_roman()
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7.5, 6.0))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(t_snap)))

    # (a) H_z(r, t)
    for i, (t, lab, col) in enumerate(zip(t_snap, t_labels, colors)):
        ax1.plot(r * 1e3, H_exact[i], '-', color=col, linewidth=1.2,
                 label=lab)
        ax1.plot(r * 1e3, H_reduced[i], '--', color=col, linewidth=0.7,
                 alpha=0.8)
    ax1.set_xlabel('radial position $r$ (mm)')
    ax1.set_ylabel(r'$H_z(r, t) / H_0$')
    ax1.legend(fontsize=7, frameon=False, loc='lower right', ncol=1)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'(a) $H_z$ (solid: $N$={N_exact}, dashed: $N$={N_reduced})',
                  fontsize=9)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')
    ax1.set_ylim(-0.1, 1.1)

    # (b) J_phi(r, t)
    for i, (t, lab, col) in enumerate(zip(t_snap, t_labels, colors)):
        ax2.plot(r * 1e3, J_exact[i] * a, '-', color=col, linewidth=1.2,
                 label=lab)
        ax2.plot(r * 1e3, J_reduced[i] * a, '--', color=col, linewidth=0.7,
                 alpha=0.8)
    ax2.set_xlabel('radial position $r$ (mm)')
    ax2.set_ylabel(r'$J_\varphi(r, t) \cdot a / H_0$')
    ax2.legend(fontsize=7, frameon=False, loc='upper left', ncol=1)
    ax2.grid(True, alpha=0.3)
    ax2.set_title(f'(b) $J_\\varphi$', fontsize=9)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    # (c) Error in H_z (reduced vs exact)
    H_err = np.abs(H_reduced - H_exact)
    for i, (t, lab, col) in enumerate(zip(t_snap, t_labels, colors)):
        ax3.semilogy(r * 1e3, H_err[i] + 1e-12, color=col, linewidth=1.0,
                     label=lab)
    ax3.set_xlabel('radial position $r$ (mm)')
    ax3.set_ylabel(r'$|H_z^{\mathrm{reduced}} - H_z^{\mathrm{exact}}|$')
    ax3.legend(fontsize=7, frameon=False, loc='lower left', ncol=1)
    ax3.grid(True, which='both', alpha=0.3)
    ax3.set_title(f'(c) Truncation error in $H_z$ ($N$={N_reduced})',
                  fontsize=9)
    for spine in ax3.spines.values():
        spine.set_linewidth(0.8)
    ax3.tick_params(direction='in', which='both')

    # (d) Error in J_phi
    J_err = np.abs(J_reduced - J_exact)
    J_scale = np.max(np.abs(J_exact))
    for i, (t, lab, col) in enumerate(zip(t_snap, t_labels, colors)):
        ax4.semilogy(r * 1e3, J_err[i] / J_scale + 1e-12, color=col,
                     linewidth=1.0, label=lab)
    ax4.set_xlabel('radial position $r$ (mm)')
    ax4.set_ylabel(r'$|J_\varphi^{\mathrm{red}} - J_\varphi^{\mathrm{ex}}| / |J|_{\max}$')
    ax4.legend(fontsize=7, frameon=False, loc='lower left', ncol=1)
    ax4.grid(True, which='both', alpha=0.3)
    ax4.set_title(f'(d) Truncation error in $J_\\varphi$', fontsize=9)
    for spine in ax4.spines.values():
        spine.set_linewidth(0.8)
    ax4.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "field_reconstruction_cylinder.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'field_reconstruction_cylinder.png'}")
    print()

    # Diagnostic: where is the truncation worst?
    for i, t in enumerate(t_snap):
        H_max_err = np.max(H_err[i])
        H_max_pos = r[np.argmax(H_err[i])] * 1e3
        J_max_err = np.max(J_err[i] / J_scale)
        J_max_pos = r[np.argmax(J_err[i])] * 1e3
        print(f"  t = {t:.2e} s:")
        print(f"    H_z max err = {H_max_err:.3e} at r = {H_max_pos:.2f} mm")
        print(f"    J_phi max err (norm) = {J_max_err:.3e} at r = {J_max_pos:.2f} mm")

    np.savez(out_dir / "field_reconstruction_data.npz",
             r=r, t_snap=t_snap, H_exact=H_exact, H_reduced=H_reduced,
             J_exact=J_exact, J_reduced=J_reduced,
             N_exact=N_exact, N_reduced=N_reduced,
             a=a, d=d)
    print(f"Saved: {out_dir / 'field_reconstruction_data.npz'}")
