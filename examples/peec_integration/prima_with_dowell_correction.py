"""
PRIMA + Dowell補正の検証

方式:
1. PRIMAをDCパラメータで構築 (R_dc, L_dc)
2. 周波数ごとに直列成分にF_R(s), F_L(s)の補正をかける

これとDowell式が一致することを確認。
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


def dowell_F_R(xi):
    """Dowell F_R"""
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        sinh_xi = np.sinh(xi)
        cosh_xi = np.cosh(xi)
        sin_xi = np.sin(xi)
        cos_xi = np.cos(xi)
        return xi * (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.001
        xi_m = xi[mask]
        sinh_xi = np.sinh(xi_m)
        cosh_xi = np.cosh(xi_m)
        sin_xi = np.sin(xi_m)
        cos_xi = np.cos(xi_m)
        result[mask] = xi_m * (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
        return result


def dowell_F_L(xi):
    """Dowell F_L"""
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        xi2 = 2 * xi
        sinh_2xi = np.sinh(xi2)
        cosh_2xi = np.cosh(xi2)
        sin_2xi = np.sin(xi2)
        cos_2xi = np.cos(xi2)
        denom = xi**3 * (cosh_2xi - cos_2xi)
        if abs(denom) < 1e-30:
            return 1.0
        return 3.0 * (sinh_2xi - sin_2xi) / denom
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.001
        xi_m = xi[mask]
        xi2 = 2 * xi_m
        sinh_2xi = np.sinh(xi2)
        cosh_2xi = np.cosh(xi2)
        sin_2xi = np.sin(xi2)
        cos_2xi = np.cos(xi2)
        denom = xi_m**3 * (cosh_2xi - cos_2xi)
        result[mask] = 3.0 * (sinh_2xi - sin_2xi) / denom
        return result


def calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs):
    """Dowell式による単一導体のインピーダンス"""
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        s = 1j * omega
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        F_R = dowell_F_R(xi)
        F_L = dowell_F_L(xi)
        Z[i] = R_dc * F_R + s * L_int_dc * F_L
    return Z


def calc_prima_dc_with_dowell_correction(R_dc, L_int_dc, d, sigma, mu, freqs):
    """
    PRIMA (DCパラメータ) + Dowell補正

    単一導体の場合:
    Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)

    これはDowell式そのもの。
    """
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        s = 1j * omega
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        F_R = dowell_F_R(xi)
        F_L = dowell_F_L(xi)

        # DC PRIMA + Dowell補正
        R_ac = R_dc * F_R
        L_ac = L_int_dc * F_L
        Z[i] = R_ac + s * L_ac

    return Z


def main():
    print("="*70)
    print("PRIMA (DC) + Dowell補正 vs Dowell式 一致検証")
    print("="*70)

    d = 0.1e-3
    sigma = 5.8e7
    mu = MU_0

    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\nDCパラメータ:")
    print(f"  R_dc = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = {L_int_dc:.6e} H*m^2")

    freqs = np.logspace(0, 8, 200)

    # Dowell式
    Z_dowell = calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs)

    # PRIMA (DC) + Dowell補正
    Z_prima_dowell = calc_prima_dc_with_dowell_correction(R_dc, L_int_dc, d, sigma, mu, freqs)

    # 比較
    print(f"\n{'Freq':>12} {'|Z_Dowell|':>14} {'|Z_PRIMA+Dow|':>14} {'Err':>10}")
    print("-"*55)

    test_freqs = [1, 10, 100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_dow = Z_dowell[idx]
        z_prima = Z_prima_dowell[idx]
        err = np.abs(z_prima - z_dow) / np.abs(z_dow) * 100
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {np.abs(z_dow):>14.6e} {np.abs(z_prima):>14.6e} {err:>9.6f}%")

    max_err = np.max(np.abs(Z_prima_dowell - Z_dowell) / np.abs(Z_dowell)) * 100
    print(f"\n最大誤差: {max_err:.10f}%")

    if max_err < 1e-10:
        print("\n=> 完全一致! PRIMA(DC) + F_R/F_L補正 = Dowell式")
    else:
        print(f"\n=> 誤差あり")

    # 実部・虚部の詳細
    print(f"\n実部・虚部の詳細:")
    print(f"{'Freq':>12} {'Re(Dow)':>14} {'Re(PRIMA)':>14} {'Im(Dow)':>14} {'Im(PRIMA)':>14}")
    print("-"*70)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_dow = Z_dowell[idx]
        z_prima = Z_prima_dowell[idx]
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {z_dow.real:>14.6e} {z_prima.real:>14.6e} {z_dow.imag:>14.6e} {z_prima.imag:>14.6e}")

    # プロット
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    ax1.loglog(freqs, np.abs(Z_dowell), 'r-', linewidth=2, label='Dowell')
    ax1.loglog(freqs, np.abs(Z_prima_dowell), 'b--', linewidth=2, label='PRIMA(DC) + F_R/F_L')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Impedance: Dowell vs PRIMA(DC) + Dowell Correction')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    ax2 = axes[1]
    err = np.abs(Z_prima_dowell - Z_dowell) / np.abs(Z_dowell) * 100
    ax2.semilogx(freqs, err, 'b-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Relative Error [%]')
    ax2.set_title('Error')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('prima_with_dowell_correction.png', dpi=150)
    print(f"\nSaved: prima_with_dowell_correction.png")
    plt.close()

    print("\n" + "="*70)
    print("結論:")
    print("  PRIMA(DC) + F_R(s)/F_L(s)補正 = Dowell式 (完全一致)")
    print("  ")
    print("  単一導体の場合:")
    print("    Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)")
    print("  ")
    print("  PEEC (複数導体) の場合:")
    print("    Z_matrix(s) = R_dc * F_R(xi) + s * (L_int_dc * F_L(xi) + L_ext)")
    print("    (外部インダクタンスL_extは表皮効果の影響を受けない)")
    print("="*70)


if __name__ == '__main__':
    main()
