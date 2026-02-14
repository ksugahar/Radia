"""
Dowell式と1D拡散式の直接比較

両者は同じ物理現象（導体中の表皮効果）を表しているので、
正しく定式化すれば一致するはず。
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    """表皮深さ delta = sqrt(2/(omega*mu*sigma))"""
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


def calc_1d_diffusion_impedance(d, sigma, mu, freqs):
    """
    1D拡散方程式の解析解（表面インピーダンス）

    Z = s*mu*(2/(k*d))*tan(k*d/2)*d
    where k = sqrt(-s*mu*sigma)

    これは単位面積あたりのインピーダンス [Ohm*m^2]
    """
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        s = 1j * 2 * np.pi * f
        k = np.sqrt(-s * mu * sigma)
        kd = k * d
        if np.abs(kd) < 1e-10:
            # DC極限
            Z[i] = d / sigma  # R_dc
        else:
            Z[i] = s * mu * (2.0 / kd) * np.tan(kd / 2.0) * d
    return Z


def dowell_F_R(xi):
    """Dowell's resistance factor F_R"""
    if np.isscalar(xi):
        if xi < 0.01:
            return 1.0
    else:
        xi = np.asarray(xi)
        result = np.ones_like(xi)
        mask = xi >= 0.01
        xi_calc = xi[mask]

        sinh_xi = np.sinh(xi_calc)
        cosh_xi = np.cosh(xi_calc)
        sin_xi = np.sin(xi_calc)
        cos_xi = np.cos(xi_calc)

        M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
        result[mask] = xi_calc * M
        return result

    # scalar case
    sinh_xi = np.sinh(xi)
    cosh_xi = np.cosh(xi)
    sin_xi = np.sin(xi)
    cos_xi = np.cos(xi)

    M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
    return xi * M


def dowell_F_L(xi):
    """Dowell's inductance factor F_L"""
    if np.isscalar(xi):
        if xi < 0.01:
            return 1.0
    else:
        xi = np.asarray(xi)
        result = np.ones_like(xi)
        mask = xi >= 0.01
        xi_calc = xi[mask]

        xi2 = 2 * xi_calc
        sinh_2xi = np.sinh(xi2)
        cosh_2xi = np.cosh(xi2)
        sin_2xi = np.sin(xi2)
        cos_2xi = np.cos(xi2)

        denom = xi_calc**3 * (cosh_2xi - cos_2xi)
        result[mask] = 3.0 * (sinh_2xi - sin_2xi) / denom
        return result

    # scalar case
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
    """
    Dowell式によるインピーダンス

    Z(f) = R_dc * F_R(xi) + j*omega * L_int_dc * F_L(xi)

    where:
      R_dc = 1/(sigma*d)  (DC抵抗、単位面積あたり) [Ohm*m^2]
      L_int_dc = mu * d / 3  (DC内部インダクタンス、単位面積あたり) [H*m^2]
      xi = d / delta = d * sqrt(omega*mu*sigma/2)
    """
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
    print("Dowell式 vs 1D拡散式 直接比較")
    print("="*70)

    # 物理パラメータ
    d = 0.1e-3      # 導体厚さ 0.1 mm
    sigma = 5.8e7   # 銅の導電率
    mu = MU_0

    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\n物理パラメータ:")
    print(f"  導体厚さ d = {d*1e3:.3f} mm")
    print(f"  導電率 sigma = {sigma:.2e} S/m")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.4e} Ohm*m^2")
    print(f"  L_int_dc = mu*d/3 = {L_int_dc:.4e} H*m^2")

    # 周波数範囲
    freqs = np.logspace(1, 8, 200)  # 10 Hz to 100 MHz

    # インピーダンス計算
    Z_diffusion = calc_1d_diffusion_impedance(d, sigma, mu, freqs)
    Z_dowell = calc_dowell_impedance(d, sigma, mu, freqs)

    # 比較
    print(f"\n{'Freq':>12} {'d/delta':>8} {'|Z_diff|':>14} {'|Z_dowell|':>14} {'Error':>10}")
    print("-"*65)

    test_freqs = [100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta

        z_diff = Z_diffusion[idx]
        z_dow = Z_dowell[idx]
        err = np.abs(z_diff - z_dow) / np.abs(z_diff) * 100

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {xi:>8.3f} {np.abs(z_diff):>14.6e} {np.abs(z_dow):>14.6e} {err:>9.4f}%")

    # 実部・虚部の比較
    print(f"\n実部・虚部の詳細比較:")
    print(f"{'Freq':>12} {'Re(Z_diff)':>14} {'Re(Z_dow)':>14} {'Im(Z_diff)':>14} {'Im(Z_dow)':>14}")
    print("-"*75)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_diff = Z_diffusion[idx]
        z_dow = Z_dowell[idx]

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {z_diff.real:>14.6e} {z_dow.real:>14.6e} {z_diff.imag:>14.6e} {z_dow.imag:>14.6e}")

    # DC極限の確認
    print(f"\n DC極限の確認:")
    print(f"   R_dc (理論値) = {R_dc:.6e} Ohm*m^2")
    print(f"   Z_diffusion(10Hz).real = {Z_diffusion[0].real:.6e}")
    print(f"   Z_dowell(10Hz).real = {Z_dowell[0].real:.6e}")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: インピーダンス振幅
    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z_diffusion), 'b-', linewidth=2, label='1D Diffusion')
    ax1.loglog(freqs, np.abs(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot 2: 相対誤差
    ax2 = axes[0, 1]
    err = np.abs(Z_diffusion - Z_dowell) / np.abs(Z_diffusion) * 100
    ax2.semilogx(freqs, err, 'b-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Relative Error [%]')
    ax2.set_title('|Z_diff - Z_dowell| / |Z_diff|')
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: 実部
    ax3 = axes[1, 0]
    ax3.loglog(freqs, np.real(Z_diffusion), 'b-', linewidth=2, label='1D Diffusion')
    ax3.loglog(freqs, np.real(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax3.axhline(y=R_dc, color='k', linestyle=':', label=f'R_dc = {R_dc:.2e}')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z) [Ohm*m^2]')
    ax3.set_title('Real Part (Resistance)')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Plot 4: 虚部
    ax4 = axes[1, 1]
    ax4.loglog(freqs, np.imag(Z_diffusion), 'b-', linewidth=2, label='1D Diffusion')
    ax4.loglog(freqs, np.imag(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Im(Z) [Ohm*m^2]')
    ax4.set_title('Imaginary Part (Reactance)')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('verify_dowell_diffusion.png', dpi=150)
    print(f"\nSaved: verify_dowell_diffusion.png")
    plt.close()

    # 最大誤差
    max_err = np.max(err)
    print(f"\n最大相対誤差: {max_err:.4f}%")


if __name__ == '__main__':
    main()
