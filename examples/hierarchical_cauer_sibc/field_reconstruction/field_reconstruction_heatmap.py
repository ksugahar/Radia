"""
Field evolution heatmap: H_z(r, t) and J_phi(r, t) for cylinder step response.

Companion to field_reconstruction_cylinder.py.  This script visualizes the
SPATIO-TEMPORAL evolution of the reconstructed field as a 2D heatmap on
the (r, t) plane, illustrating the skin-depth penetration dynamics.

Usage:
    python field_reconstruction_heatmap.py
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


def H_z_field(r_arr, t_arr, N_modes, a, sigma, mu):
    """
    H_z(r, t) / H_0 for step input, evaluated on a 2D grid (t, r).
    Returns array shape (len(t_arr), len(r_arr)).
    """
    j0_zeros = special.jn_zeros(0, N_modes)
    j1_vals = special.j1(j0_zeros)
    xi = j0_zeros**2 / (a**2 * mu * sigma)
    A = 2.0 / (j0_zeros * j1_vals)
    # Compute J_0(j_{0,k} r/a): shape (M_r, K_modes)
    arg = j0_zeros[None, :] * r_arr[:, None] / a
    J0_arr = special.j0(arg)
    # exp(-xi_k t): shape (M_t, K_modes)
    decay = np.exp(-xi[None, :] * t_arr[:, None])
    # H = 1 - sum_k A_k * J0(j_{0,k} r/a) * exp(-xi_k t)
    # Shape: (M_t, M_r, K) -> sum over K
    H = np.zeros((len(t_arr), len(r_arr)))
    for ti in range(len(t_arr)):
        H[ti] = 1.0 - np.sum(A[None, :] * J0_arr * decay[ti, None, :], axis=1)
    return H


def J_phi_field(r_arr, t_arr, N_modes, a, sigma, mu):
    """J_phi(r, t) / H_0 on (t, r) grid."""
    j0_zeros = special.jn_zeros(0, N_modes)
    j1_vals = special.j1(j0_zeros)
    xi = j0_zeros**2 / (a**2 * mu * sigma)
    B = 2.0 / (a * j1_vals)
    arg = j0_zeros[None, :] * r_arr[:, None] / a
    J1_arr = special.j1(arg)
    decay = np.exp(-xi[None, :] * t_arr[:, None])
    J = np.zeros((len(t_arr), len(r_arr)))
    for ti in range(len(t_arr)):
        J[ti] = np.sum(B[None, :] * J1_arr * decay[ti, None, :], axis=1)
    return J


def skin_depth_estimator(omega, sigma, mu):
    """Classical skin depth delta(omega) = sqrt(2/(omega mu sigma))."""
    return np.sqrt(2.0 / (omega * mu * sigma))


if __name__ == "__main__":
    # Cu cylinder parameters
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    a = 5e-3
    d = 2 * np.pi * 7e4
    t_d = 1.0 / d

    # Grid for heatmap
    r_arr = np.linspace(0, a, 150)
    t_arr = np.logspace(np.log10(t_d * 1e-3), np.log10(t_d * 1e3), 200)

    print(f"Computing field on {len(t_arr)} x {len(r_arr)} grid...")
    N = 200
    H = H_z_field(r_arr, t_arr, N, a, sigma_cu, mu_0)
    J = J_phi_field(r_arr, t_arr, N, a, sigma_cu, mu_0)
    print(f"  done.  H range: [{H.min():.3e}, {H.max():.3e}]")
    print(f"          J range: [{J.min():.3e}, {J.max():.3e}]")
    print()

    # Plot heatmap
    setup_times_new_roman()
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.4))

    R, T = np.meshgrid(r_arr * 1e3, t_arr)

    # H_z heatmap
    im1 = axes[0].pcolormesh(R, T, H, cmap='RdBu_r', shading='auto',
                              vmin=0, vmax=1)
    axes[0].set_yscale('log')
    axes[0].set_xlabel('radial position $r$ (mm)')
    axes[0].set_ylabel('time $t$ (s)')
    axes[0].set_title('(a) $H_z(r, t) / H_0$', fontsize=10)
    axes[0].axhline(t_d, color='black', linestyle='--', linewidth=0.6,
                    alpha=0.6, label=r'$t = 1/d$')
    axes[0].axvline(a * 1e3, color='gray', linewidth=0.4, alpha=0.5)
    axes[0].legend(loc='lower right', fontsize=8, frameon=False,
                   labelcolor='black')
    cbar1 = plt.colorbar(im1, ax=axes[0], shrink=0.85)
    cbar1.ax.tick_params(direction='in')

    # J_phi heatmap (signed, log-magnitude)
    J_max = np.max(np.abs(J))
    im2 = axes[1].pcolormesh(R, T, J / J_max, cmap='RdBu_r', shading='auto',
                              vmin=-1, vmax=1)
    axes[1].set_yscale('log')
    axes[1].set_xlabel('radial position $r$ (mm)')
    axes[1].set_ylabel('time $t$ (s)')
    axes[1].set_title('(b) $J_\\varphi(r, t)$ (normalised)', fontsize=10)
    axes[1].axhline(t_d, color='black', linestyle='--', linewidth=0.6,
                    alpha=0.6, label=r'$t = 1/d$')
    axes[1].axvline(a * 1e3, color='gray', linewidth=0.4, alpha=0.5)
    axes[1].legend(loc='lower right', fontsize=8, frameon=False)
    cbar2 = plt.colorbar(im2, ax=axes[1], shrink=0.85)
    cbar2.ax.tick_params(direction='in')

    for ax in axes:
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "field_reconstruction_heatmap.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / 'field_reconstruction_heatmap.png'}")

    # Also compute and plot the time evolution of skin depth boundary
    # (where H_z drops below 1/e of surface value, treating r->a as surface)
    print()
    print("Skin penetration analysis:")
    for t_check in [0.01 * t_d, 0.1 * t_d, t_d, 10 * t_d]:
        ti = np.argmin(np.abs(t_arr - t_check))
        H_slice = H[ti]
        # Find where H drops to 1/e of H(r=0)
        if H_slice[0] > 0:
            target = H_slice[0] * (1 - 1/np.e) + 0.0  # not quite skin depth def
        else:
            target = 1 - 1/np.e  # threshold from surface
        # Find penetration depth from surface: r where H crosses (1 - 1/e)
        idx = np.where(H_slice > 1 - 1/np.e)[0]
        if len(idx) > 0:
            r_pen = r_arr[idx[0]]
            pen = (a - r_pen) * 1e3  # depth from surface in mm
            print(f"  t = {t_check:.2e} s: H drops below (1 - 1/e) at depth {pen:.2f} mm from surface")
        else:
            print(f"  t = {t_check:.2e} s: H not penetrated yet")

    np.savez(out_dir / "field_heatmap_data.npz",
             r=r_arr, t=t_arr, H=H, J=J, N_modes=N, a=a, d=d)
    print(f"Saved: {out_dir / 'field_heatmap_data.npz'}")
