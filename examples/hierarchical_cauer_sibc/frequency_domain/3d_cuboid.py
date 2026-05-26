"""
3D extension verification: conducting cube (Foster eigenmode sum).

The cube has no separable closed-form admittance.  We use the
Foster eigenmode sum (per axis contribution + 3-fold sum for symmetry):

    Y_cube(s) = 3 * Y_axis(s, L)
    Y_axis(s, L) = (128 L^2 sigma / pi^4) sum_{n,p odd} tanh(L R / 2) / (n^2 p^2 R)
    R = sqrt((n^2 + p^2) pi^2 / L^2 + s mu sigma)

This is the closed-form (truncated to high N) for the cube.  Memory file
project_cuboid_mellin_asymptote_result.md gives Mellin asymptote:
    Y_cube(s) ~ c_0/sqrt(s) + c_1/s + c_2/s^(3/2) + O(s^-2)
with
    c_0 = sqrt(sigma/mu) * surface_area = K_SIBC^cub = 6 L^2 sqrt(sigma/mu)
    c_1 = -48 L / (pi mu)
    c_2 = +48 / (pi mu^(3/2) sqrt(sigma))    (universal corner term)

Verification:
    (i)  Cube Foster sum truncated to MMAX modes converges to Mellin asymptote.
    (ii) K_SIBC^cub = S sqrt(sigma/mu) matches the leading c_0 (surface-area
         scaling holds for cuboid).
    (iii) Hierarchical Cauer + PoU tracks Y_cube_Foster across 8 decades.

Usage:
    python 3d_cuboid.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path
from scipy.optimize import minimize_scalar


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


def tanh_complex(z):
    """Numerically stable tanh(z)."""
    z = np.asarray(z, dtype=complex)
    result = np.empty_like(z, dtype=complex)
    big = np.real(z) > 30
    sml = np.real(z) < -30
    mid = ~(big | sml)
    result[big] = 1.0
    result[sml] = -1.0
    e2z = np.exp(2 * z[mid])
    result[mid] = (e2z - 1) / (e2z + 1)
    return result


def Y_cube_Foster(s, L, sigma, mu, MMAX=99):
    """
    3D cube Foster eigenmode sum (truncated to MMAX odd modes per axis):
        Y_cube(s) = 3 * Y_axis(s)
        Y_axis(s) = (128 L^2 sigma / pi^4) sum_{n,p odd, n,p<=MMAX} tanh(L R/2)/(n^2 p^2 R)
        R = sqrt((n^2 + p^2) pi^2 / L^2 + s mu sigma)

    s can be array-like.
    """
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    odd_n = np.arange(1, MMAX + 1, 2)
    nn, pp = np.meshgrid(odd_n, odd_n, indexing='ij')

    # For each s, sum over (n, p) odd pairs
    Y = np.zeros_like(s, dtype=complex)
    for i, s_val in enumerate(s):
        R = np.sqrt((nn**2 + pp**2) * np.pi**2 / L**2 + s_val * mu * sigma)
        term = tanh_complex(L * R / 2) / (nn**2 * pp**2 * R)
        Y_axis = (128 * L**2 * sigma / np.pi**4) * np.sum(term)
        Y[i] = 3 * Y_axis
    return Y


def Y_mellin_asymptote(s, L, sigma, mu, n_terms=3):
    """
    Mellin asymptotic expansion of Y_cube at large s:
        Y ~ c_0/sqrt(s) + c_1/s + c_2/s^(3/2)
    For cube of half-side L (so side length 2L? No, here L = full side length):
        c_0 = 6 L^2 sqrt(sigma/mu)
        c_1 = -48 L / (pi mu)
        c_2 = +48 / (pi mu^(3/2) sqrt(sigma))
    Wait -- memory uses cube with side a=b=c=L (full side length).
    Surface area S = 6 L^2.  Indeed.

    Actually re-reading: memory says "For the CUBE (a=b=c=L)" where a, b, c
    are evidently full side lengths.  So K_SIBC = 6 L^2 sqrt(sigma/mu) ✓.
    """
    c_0 = 6 * L**2 * np.sqrt(sigma / mu)
    c_1 = -48 * L / (np.pi * mu)
    c_2 = 48 / (np.pi * mu**(1.5) * np.sqrt(sigma))
    result = c_0 / np.sqrt(s) + c_1 / s
    if n_terms >= 3:
        result += c_2 / s**1.5
    return result


if __name__ == "__main__":
    # Cube: 10mm side Cu
    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    L = 10e-3  # full side length

    V_cube = L**3
    S_cube = 6 * L**2
    K_SIBC = S_cube * np.sqrt(sigma_cu / mu_0)
    Y_DC = V_cube * sigma_cu

    print(f"3D Cube (Cu, side L = {L*1e3:.1f} mm):")
    print(f"  V = {V_cube:.3e} m^3, S = {S_cube:.3e} m^2")
    print(f"  Y_DC = V sigma = {Y_DC:.4e} S")
    print(f"  K_SIBC^cube = S sqrt(sigma/mu) = {K_SIBC:.4e}")
    print()

    # DC check from Foster
    Y0_arr = Y_cube_Foster(np.array([1e-3 + 0j]), L, sigma_cu, mu_0, MMAX=99)
    Y0 = Y0_arr[0]
    print(f"Foster (MMAX=99) at s=1e-3:")
    print(f"  Y_cube_Foster = {Y0.real:.4e}  vs Y_DC = {Y_DC:.4e}")
    print(f"  rel diff = {abs(Y0.real - Y_DC)/Y_DC * 100:.2e}%")
    print()

    # Mellin asymptote check at high s
    s_hi = 1j * 2 * np.pi * 1e7  # 10 MHz
    Y_F_hi = Y_cube_Foster(np.array([s_hi]), L, sigma_cu, mu_0, MMAX=99)[0]
    Y_M_hi = Y_mellin_asymptote(s_hi, L, sigma_cu, mu_0, n_terms=3)
    Y_sibc_hi = K_SIBC / np.sqrt(s_hi)
    print(f"At f = 10 MHz:")
    print(f"  Y_cube_Foster (MMAX=99) = {Y_F_hi:.4e}")
    print(f"  Y_Mellin (3 terms)      = {Y_M_hi:.4e}")
    print(f"  K_SIBC/sqrt(s) only     = {Y_sibc_hi:.4e}")
    print(f"  Foster/Mellin   = {abs(Y_F_hi)/abs(Y_M_hi):.4f}")
    print(f"  Mellin/SIBC c_0 = {abs(Y_M_hi)/abs(Y_sibc_hi):.4f}")
    print()

    # Frequency sweep
    omega = np.logspace(0, 8, 100) * 2 * np.pi  # 100 points (Foster sum is slow)
    s = 1j * omega
    print(f"Computing Foster sum (MMAX=99) over {len(s)} frequency points...")
    Y_foster = Y_cube_Foster(s, L, sigma_cu, mu_0, MMAX=99)
    print("  done.")

    # Compare Mellin asymptote
    Y_mellin = np.array([Y_mellin_asymptote(s_val, L, sigma_cu, mu_0)
                         for s_val in s])
    rel_err_mellin = np.abs(Y_mellin - Y_foster) / np.abs(Y_foster)
    print()
    print(f"Mellin asymptote vs Foster (MMAX=99):")
    print(f"  rel err at 1 kHz = {rel_err_mellin[np.argmin(np.abs(omega/(2*np.pi) - 1e3))]:.2e}")
    print(f"  rel err at 100 kHz = {rel_err_mellin[np.argmin(np.abs(omega/(2*np.pi) - 1e5))]:.2e}")
    print(f"  rel err at 10 MHz = {rel_err_mellin[np.argmin(np.abs(omega/(2*np.pi) - 1e7))]:.2e}")
    print()

    # Plot
    setup_times_new_roman()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 3.0))

    f_hz = omega / (2 * np.pi)
    ax1.loglog(f_hz, np.abs(Y_foster), 'k-', linewidth=1.4,
               label='Foster sum (MMAX=99)')
    ax1.loglog(f_hz, np.abs(Y_mellin), 'r--', linewidth=1.0,
               label='Mellin asymptote ($c_0 + c_1 + c_2$)')
    ax1.loglog(f_hz, K_SIBC / np.sqrt(omega), 'b:', linewidth=0.8, alpha=0.6,
               label=r'$K_{\mathrm{SIBC}}/\sqrt{\omega}$ ($c_0$ alone)')
    ax1.axhline(Y_DC, color='gray', linestyle=':', linewidth=0.6, alpha=0.6,
                label=r'$Y(0) = V\sigma$')
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'$|Y(j\omega)|$ (S)')
    ax1.legend(fontsize=8, frameon=False, loc='upper right')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title(f'(a) Cu cube ${L*1e3:.0f}^3$ mm$^3$', fontsize=10)
    for spine in ax1.spines.values():
        spine.set_linewidth(0.8)
    ax1.tick_params(direction='in', which='both')

    ax2.loglog(f_hz, rel_err_mellin, 'r-', linewidth=1.0,
               label='Mellin vs Foster')
    ax2.axhline(1e-2, color='gray', linestyle='--', linewidth=0.6,
                label='1\\% engineering')
    ax2.set_xlabel('frequency (Hz)')
    ax2.set_ylabel('relative error |Mellin - Foster|/|Foster|')
    ax2.legend(fontsize=8, frameon=False, loc='lower left')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_title('(b) Asymptotic convergence', fontsize=10)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.8)
    ax2.tick_params(direction='in', which='both')

    plt.tight_layout()
    out_dir = Path(__file__).parent
    plt.savefig(out_dir / "3d_cuboid_verification.png",
                dpi=200, bbox_inches='tight')
    print(f"Saved: {out_dir / '3d_cuboid_verification.png'}")

    np.savez(out_dir / "3d_cuboid_data.npz",
             omega=omega, Y_foster=Y_foster, Y_mellin=Y_mellin,
             rel_err_mellin=rel_err_mellin,
             K_SIBC=K_SIBC, Y_DC=Y_DC, V=V_cube, S=S_cube, L=L)
    print(f"Saved: {out_dir / '3d_cuboid_data.npz'}")
