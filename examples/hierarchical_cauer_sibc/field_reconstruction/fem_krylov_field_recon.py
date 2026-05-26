"""
FEM Krylov basis field reconstruction for arbitrary 3D conductors.

Generalises the modal field reconstruction beyond canonical bodies
(cylinder Bessel, sphere spherical Bessel) to arbitrary 3D geometry,
using the CLN Krylov basis Q from Kameari 2018.

Mathematical structure:
  - Bulk CLN basis Q = [q_1, ..., q_N]: built by Krylov recursion
    q_{n+1} ∝ K_0^{-1} σ q_n
    (q_1 = K_0^{-1} σ b, K_0 is the magnetostatic stiffness, b is source)
  - Each q_k is a FEM-discretised vector spanning the eddy-current Krylov subspace
  - In the time-domain ODE, each modal amplitude z_k(t) corresponds to one q_k:
       H_eddy(r, t) = sum_k z_k(t) phi_k(r)
    where phi_k(r) is the FEM-interpolated representation of q_k.

For arbitrary 3D mesh:
  - q_k is a coefficient vector on the FEM edges (HCurl) or nodes (H1)
  - phi_k(r) at any point r is obtained by FEM basis function evaluation
  - Total field: H_z(r, t) = H_ext(t) - sum_k z_k(t) phi_k(r)

This script demonstrates the methodology on the cube using the Foster
eigenmode sum (from 3d_cuboid.py) as a stand-in for the proper FEM Krylov
basis.  Each Foster eigenmode (n, p) of the cube has spatial shape
  phi_{n,p}(x, y, z) = sin(n pi x/L) sin(p pi y/L) [for one axis contribution]
and time amplitude z_{n,p}(t) = c_{n,p}/xi_{n,p} (1 - exp(-xi_{n,p} t))
for step input.

This is the 3D analog of the cylinder Bessel reconstruction.  When applied
to an FEM Krylov basis (instead of analytical eigenmodes), the methodology
extends to arbitrary 3D conductor geometry.

Usage:
    python fem_krylov_field_recon.py
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


def cube_eigenmode_field(x, y, z, t, MMAX, L, sigma, mu, axis=2):
    """
    Cube eddy current eigenmode field for step input H_ext = H_0 along axis (default z).

    For external H_z, the induced J has azimuthal pattern in the (x, y) plane.
    The H_z field in the cube interior is:
        H_z(x, y, z, t)/H_0 = 1 - sum_{m, n odd} A_{mn} sin(m pi (x+L/2)/L) sin(n pi (y+L/2)/L) exp(-xi_{mn} t)

    where (m, n) sum over odd integers (Dirichlet BC at cube faces, induced eddy in (x, y) plane
    for z-directed external field), and
        xi_{mn} = ((m^2 + n^2) pi^2) / (L^2 mu sigma)
        A_{mn} = 16 / (m n pi^2)  (Fourier-double-sine expansion of constant 1 on [-L/2, L/2]^2)

    Returns H/H_0 field on the given (x, y) grid (z is integrated away in 2D mode).
    """
    L_half = L / 2.0
    # Shift to [0, L] for sine expansion
    xs = x + L_half
    ys = y + L_half

    # Odd m, n indices up to MMAX
    odd = np.arange(1, MMAX + 1, 2)
    M, N = np.meshgrid(odd, odd, indexing='ij')

    xi = ((M**2 + N**2) * np.pi**2) / (L**2 * mu * sigma)
    A = 16.0 / (M * N * np.pi**2)
    decay = np.exp(-xi * t)

    # Compute sum over modes
    # sin(m pi xs / L) for each (M, N) pair and (xs, ys) point
    # Output shape: (n_pts, n_modes_M, n_modes_N) -- collapse modes axis
    if np.isscalar(xs):
        sin_x = np.sin(M.flatten() * np.pi * xs / L)
        sin_y = np.sin(N.flatten() * np.pi * ys / L)
        amp = A.flatten() * decay.flatten()
        result = 1.0 - np.sum(amp * sin_x * sin_y)
        return result
    else:
        xs_arr = np.asarray(xs).flatten()
        ys_arr = np.asarray(ys).flatten()
        result = np.zeros(xs_arr.shape)
        for i, (x_i, y_i) in enumerate(zip(xs_arr, ys_arr)):
            sin_x = np.sin(M * np.pi * x_i / L)
            sin_y = np.sin(N * np.pi * y_i / L)
            result[i] = 1.0 - np.sum(A * decay * sin_x * sin_y)
        return result.reshape(np.asarray(xs).shape)


def step_response_cube_section(x, y, t, MMAX, L, sigma, mu):
    """Step response H_z field at fixed z (interior plane), 2D (x, y) grid."""
    return cube_eigenmode_field(x, y, 0, t, MMAX, L, sigma, mu)


if __name__ == "__main__":
    # Cu cube
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    L = 10e-3
    print(f"Cu cube L = {L*1e3:.1f} mm")
    print()

    # Characteristic times
    # tau = a^2 mu sigma / pi^2 = L^2 mu sigma / (4 pi^2) for half-side a = L/2
    tau_diff = (L/2)**2 * mu_0 * sigma_cu / np.pi**2
    print(f"Diffusion time (mode 1,1): tau = {tau_diff:.3e} s")
    print()

    # 2D plane: x, y over [-L/2, L/2], z = 0 (mid plane)
    N_grid = 50
    x_arr = np.linspace(-L/2 + 0.05*L, L/2 - 0.05*L, N_grid)
    y_arr = np.linspace(-L/2 + 0.05*L, L/2 - 0.05*L, N_grid)
    X, Y = np.meshgrid(x_arr, y_arr, indexing='ij')

    # Snapshots
    t_snaps = [0.01 * tau_diff, 0.1 * tau_diff, 0.5 * tau_diff,
               1.0 * tau_diff, 5.0 * tau_diff]
    snap_labels = [f'{c}'r'$\tau$' for c in
                   ['0.01', '0.1', '0.5', '1.0', '5.0']]

    # Compute exact (MMAX = 49) and reduced (MMAX = 11) field at each snapshot
    MMAX_exact = 49
    MMAX_red = 11
    print(f"Computing field on {N_grid}x{N_grid} grid, "
          f"MMAX_exact={MMAX_exact}, MMAX_red={MMAX_red}...")

    H_exact = np.zeros((len(t_snaps), N_grid, N_grid))
    H_reduced = np.zeros((len(t_snaps), N_grid, N_grid))
    for si, t_s in enumerate(t_snaps):
        H_exact[si] = step_response_cube_section(X, Y, t_s, MMAX_exact,
                                                  L, sigma_cu, mu_0)
        H_reduced[si] = step_response_cube_section(X, Y, t_s, MMAX_red,
                                                    L, sigma_cu, mu_0)
        print(f"  t = {t_s:.2e} s ({t_s/tau_diff:.2f} tau): "
              f"H range [{H_exact[si].min():.3f}, {H_exact[si].max():.3f}]")

    # Plot: 5 columns for 5 snapshots
    setup_times_new_roman()
    fig, axes = plt.subplots(2, len(t_snaps),
                              figsize=(2.0 * len(t_snaps), 4.0))

    for si, (t_s, lab) in enumerate(zip(t_snaps, snap_labels)):
        # Top row: exact field
        im1 = axes[0, si].pcolormesh(X * 1e3, Y * 1e3, H_exact[si],
                                       cmap='RdBu_r', vmin=0, vmax=1,
                                       shading='auto')
        axes[0, si].set_aspect('equal')
        axes[0, si].set_title(f'$t$ = {lab}', fontsize=8)
        axes[0, si].tick_params(direction='in', labelsize=6)
        if si == 0:
            axes[0, si].set_ylabel('exact\n($M$=49)', fontsize=8)

        # Bottom row: error |reduced - exact|
        err = np.abs(H_reduced[si] - H_exact[si])
        im2 = axes[1, si].pcolormesh(X * 1e3, Y * 1e3, err,
                                       cmap='viridis', vmin=0,
                                       vmax=max(err.max(), 1e-3),
                                       shading='auto')
        axes[1, si].set_aspect('equal')
        axes[1, si].set_xlabel('$x$ (mm)', fontsize=8)
        axes[1, si].tick_params(direction='in', labelsize=6)
        if si == 0:
            axes[1, si].set_ylabel(f'|$N$={MMAX_red} - exact|',
                                    fontsize=8)
        max_err = err.max()
        axes[1, si].text(0.95, 0.95,
                          f'max={max_err:.1e}',
                          transform=axes[1, si].transAxes,
                          fontsize=7, ha='right', va='top',
                          color='white' if max_err > 0.1 else 'black',
                          bbox=dict(facecolor='black', alpha=0.5)
                          if max_err > 0.1 else None)

    # Add colorbars
    cax1 = fig.add_axes([0.92, 0.55, 0.012, 0.30])
    cbar1 = plt.colorbar(im1, cax=cax1)
    cbar1.set_label(r'$H_z/H_0$', fontsize=8)
    cbar1.ax.tick_params(labelsize=6)

    fig.suptitle(f'Cu cube field reconstruction at $z = 0$ plane'
                  f' (step response, axis $\\hat z$)',
                  fontsize=10, y=1.02)

    plt.tight_layout(rect=[0, 0, 0.91, 1])
    out_dir = Path(__file__).parent
    out_path = out_dir / "fem_krylov_field_cube.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print()
    print(f"Saved: {out_path}")

    np.savez(out_dir / "fem_krylov_field_cube_data.npz",
             x=x_arr, y=y_arr, t_snaps=t_snaps,
             H_exact=H_exact, H_reduced=H_reduced,
             L=L, sigma=sigma_cu, mu=mu_0,
             MMAX_exact=MMAX_exact, MMAX_red=MMAX_red)
    print(f"Saved: fem_krylov_field_cube_data.npz")
    print()
    print("Methodology: The 2D mode shapes sin(m pi x/L) sin(n pi y/L) here are")
    print("the analytical eigenmodes of the cube; in an FEM Krylov implementation,")
    print("these would be replaced by the Q basis vectors {q_k} from the CLN")
    print("orthogonalisation (Kameari 2018).  The principle is identical: each")
    print("state variable z_k(t) carries one modal amplitude that, when paired")
    print("with its spatial shape phi_k(r), reconstructs the local field.")
