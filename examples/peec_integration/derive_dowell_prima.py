"""
Dowell式からPRIMAパラメータを導出

Dowell式: Z = R_dc * F_R(xi) + j*omega * L_int_dc * F_L(xi)

これをPRIMA I形式で表現できるか検討。

問題: 1D拡散式 Z = s*mu*(2/(k*d))*tan(k*d/2)*d は
      DC極限で Z -> j*omega*mu*d （抵抗成分なし）

      一方、Dowell式はDC極限で Z -> R_dc + j*omega*L_int_dc

この違いは何か？
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


# =============================================
# 1D拡散方程式の解（両面対称境界条件）
# =============================================
def calc_1d_diffusion_symmetric(d, sigma, mu, freqs):
    """
    両面対称境界条件の1D拡散
    厚さdの導体板、両面からH場が侵入

    Z = s*mu*d * coth(gamma*d/2) / (gamma*d/2)
    where gamma = sqrt(j*omega*mu*sigma)

    DC極限: Z -> 2/(sigma*d) = 2*R_dc (両面から見た抵抗)
    """
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        s = 1j * omega
        gamma = np.sqrt(s * mu * sigma)
        gd2 = gamma * d / 2

        if np.abs(gd2) < 1e-10:
            # DC極限: coth(x)/x -> 1/x^2 + 1/3 for small x
            # Z -> s*mu*d * (2/(gamma*d)^2 + 1/3)
            #    = s*mu*d * (2/(s*mu*sigma*d^2) + 1/3)
            #    = 2/(sigma*d) + s*mu*d/3
            Z[i] = 2.0 / (sigma * d) + s * mu * d / 3.0
        else:
            Z[i] = s * mu * d * np.cosh(gd2) / (gd2 * np.sinh(gd2))

    return Z


# =============================================
# 1D拡散方程式の解（片面境界条件 - 元の式）
# =============================================
def calc_1d_diffusion_onesided(d, sigma, mu, freqs):
    """
    片面境界条件の1D拡散

    Z = s*mu*(2/(k*d))*tan(k*d/2)*d
    where k = sqrt(-s*mu*sigma)

    DC極限: Z -> j*omega*mu*d (抵抗なし)
    """
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        s = 1j * 2 * np.pi * f
        k = np.sqrt(-s * mu * sigma)
        kd = k * d

        if np.abs(kd) < 1e-10:
            Z[i] = s * mu * d
        else:
            Z[i] = s * mu * (2.0 / kd) * np.tan(kd / 2.0) * d

    return Z


# =============================================
# Dowell式
# =============================================
def dowell_F_R(xi):
    if xi < 0.001:
        return 1.0 + xi**4 / 45.0
    sinh_xi = np.sinh(xi)
    cosh_xi = np.cosh(xi)
    sin_xi = np.sin(xi)
    cos_xi = np.cos(xi)
    M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
    return xi * M


def dowell_F_L(xi):
    if xi < 0.001:
        return 1.0 - xi**4 / 45.0
    xi2 = 2 * xi
    sinh_2xi = np.sinh(xi2)
    cosh_2xi = np.cosh(xi2)
    sin_2xi = np.sin(xi2)
    cos_2xi = np.cos(xi2)
    denom = xi**3 * (cosh_2xi - cos_2xi)
    if abs(denom) < 1e-30:
        return 1.0
    return 3.0 * (sinh_2xi - sin_2xi) / denom


def calc_dowell_impedance(d, sigma, mu, freqs):
    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        F_R = dowell_F_R(xi)
        F_L = dowell_F_L(xi)
        Z[i] = R_dc * F_R + 1j * omega * L_int_dc * F_L

    return Z


def main():
    print("="*70)
    print("1D拡散式とDowell式の比較")
    print("="*70)

    d = 0.1e-3
    sigma = 5.8e7
    mu = MU_0

    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\n物理パラメータ:")
    print(f"  d = {d*1e6:.1f} um")
    print(f"  sigma = {sigma:.2e} S/m")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = mu*d/3 = {L_int_dc:.6e} H*m^2")
    print(f"  mu*d = {mu*d:.6e} H*m^2")

    freqs = np.logspace(0, 8, 200)

    Z_onesided = calc_1d_diffusion_onesided(d, sigma, mu, freqs)
    Z_symmetric = calc_1d_diffusion_symmetric(d, sigma, mu, freqs)
    Z_dowell = calc_dowell_impedance(d, sigma, mu, freqs)

    print(f"\n{'Freq':>12} {'|Z_1side|':>14} {'|Z_symm|':>14} {'|Z_Dowell|':>14}")
    print("-"*60)

    test_freqs = [1, 10, 100, 1e3, 10e3, 100e3, 1e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {np.abs(Z_onesided[idx]):>14.6e} {np.abs(Z_symmetric[idx]):>14.6e} {np.abs(Z_dowell[idx]):>14.6e}")

    # Dowell vs Symmetric の比較
    print(f"\nDowell vs Symmetric (両面対称) 比較:")
    print(f"{'Freq':>12} {'Re(Symm)':>14} {'Re(Dow)':>14} {'Im(Symm)':>14} {'Im(Dow)':>14} {'Err':>10}")
    print("-"*80)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_sym = Z_symmetric[idx]
        z_dow = Z_dowell[idx]
        err = np.abs(z_sym - z_dow) / np.abs(z_dow) * 100
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {z_sym.real:>14.6e} {z_dow.real:>14.6e} {z_sym.imag:>14.6e} {z_dow.imag:>14.6e} {err:>9.2f}%")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z_onesided), 'b-', linewidth=2, label='1D Diffusion (one-sided)')
    ax1.loglog(freqs, np.abs(Z_symmetric), 'g--', linewidth=2, label='1D Diffusion (symmetric)')
    ax1.loglog(freqs, np.abs(Z_dowell), 'r:', linewidth=2, label='Dowell')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    ax2 = axes[0, 1]
    err_sym = np.abs(Z_symmetric - Z_dowell) / np.abs(Z_dowell) * 100
    ax2.semilogx(freqs, err_sym, 'g-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Relative Error [%]')
    ax2.set_title('Symmetric vs Dowell Error')
    ax2.grid(True, which='both', alpha=0.3)

    ax3 = axes[1, 0]
    ax3.loglog(freqs, np.real(Z_symmetric), 'g--', linewidth=2, label='Symmetric')
    ax3.loglog(freqs, np.real(Z_dowell), 'r:', linewidth=2, label='Dowell')
    ax3.axhline(y=R_dc, color='k', linestyle=':', alpha=0.5, label='R_dc')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z)')
    ax3.set_title('Real Part')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    ax4 = axes[1, 1]
    ax4.loglog(freqs, np.imag(Z_symmetric), 'g--', linewidth=2, label='Symmetric')
    ax4.loglog(freqs, np.imag(Z_dowell), 'r:', linewidth=2, label='Dowell')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Im(Z)')
    ax4.set_title('Imaginary Part')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('derive_dowell_prima.png', dpi=150)
    print(f"\nSaved: derive_dowell_prima.png")
    plt.close()


if __name__ == '__main__':
    main()
