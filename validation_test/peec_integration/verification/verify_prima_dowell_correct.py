"""
PRIMA vs Dowell 正しい比較

Dowell式:
  Z(f) = R_dc * F_R(xi) + j*omega * L_int_dc * F_L(xi)

PRIMA I 式:
  表皮効果による内部インピーダンスを表現
  DC極限で Z -> R_dc + j*omega*L_int_dc

両者が一致するか確認。
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    """表皮深さ"""
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


def skin_effect_prima_params(d, sigma, mu=MU_0, n_stages=7):
    """PRIMA I パラメータ"""
    R_prima = np.zeros(n_stages)
    L_prima = np.zeros(n_stages)

    for n in range(1, n_stages + 1):
        if n == 1:
            R_prima[n-1] = 1e-16
        else:
            R_prima[n-1] = (4*n - 5) * 4.0 / (sigma * d)
        L_prima[n-1] = d * mu / (4*n - 3)

    return R_prima, L_prima


def build_tridiagonal(L):
    """三重対角行列を構築"""
    n = len(L)
    U = np.eye(n)
    for i in range(n - 1):
        U[i, i + 1] = -1
    return U.T @ np.diag(L) @ U


def calc_prima_impedance(R_prima, L_prima, freqs):
    """PRIMA I インピーダンス"""
    RR = np.diag(R_prima)
    LL = build_tridiagonal(L_prima)
    n = len(R_prima)
    V = np.zeros(n)
    V[0] = 1.0

    Z = np.zeros(len(freqs), dtype=complex)
    for i, freq in enumerate(freqs):
        s = 1j * 2 * np.pi * freq
        ZZ = RR + s * LL
        I = np.linalg.solve(ZZ, V)
        Z[i] = 1.0 / I[0]

    return Z


def dowell_F_R(xi):
    """Dowell F_R (単層導体)"""
    if xi < 0.001:
        # Taylor展開: F_R ~ 1 + xi^4/45 + ...
        return 1.0 + xi**4 / 45.0
    sinh_xi = np.sinh(xi)
    cosh_xi = np.cosh(xi)
    sin_xi = np.sin(xi)
    cos_xi = np.cos(xi)
    M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
    return xi * M


def dowell_F_L(xi):
    """Dowell F_L (単層導体の内部インダクタンス)"""
    if xi < 0.001:
        # Taylor展開: F_L ~ 1 - xi^4/45 + ...
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


def calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs):
    """Dowell式によるインピーダンス"""
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
    print("PRIMA vs Dowell 正しい比較")
    print("="*70)

    # 物理パラメータ
    d = 0.1e-3      # 0.1 mm
    sigma = 5.8e7   # 銅
    mu = MU_0

    # DC パラメータ（単位面積あたり）
    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\n物理パラメータ:")
    print(f"  d = {d*1e3:.3f} mm")
    print(f"  sigma = {sigma:.2e} S/m")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = mu*d/3 = {L_int_dc:.6e} H*m^2")

    # PRIMA パラメータ
    n_stages = 20  # 多めに取る
    R_prima, L_prima = skin_effect_prima_params(d, sigma, mu, n_stages)

    # PRIMAのDC極限を確認
    # DC極限では全てのLは短絡、Rは直列
    R_total_prima = np.sum(R_prima)
    L_total_prima = L_prima[0]  # 最初のLが支配的（他は並列で小さくなる）

    print(f"\nPRIMA DC極限:")
    print(f"  sum(R_prima) = {R_total_prima:.6e} (R_dcと比較)")
    print(f"  L_prima[0] = {L_prima[0]:.6e} (L_int_dcと比較: {L_int_dc:.6e})")

    # L_prima[0] = d*mu/(4*1-3) = d*mu/1 = d*mu
    # L_int_dc = mu*d/3
    # 比: L_prima[0] / L_int_dc = 3
    print(f"  L_prima[0] / L_int_dc = {L_prima[0] / L_int_dc:.2f}")

    # 周波数範囲
    freqs = np.logspace(1, 8, 200)

    # インピーダンス計算
    Z_prima = calc_prima_impedance(R_prima, L_prima, freqs)
    Z_dowell = calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs)

    # 比較
    print(f"\n{'Freq':>12} {'xi':>8} {'|Z_PRIMA|':>14} {'|Z_Dow|':>14} {'Err':>10}")
    print("-"*65)

    test_freqs = [10, 100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta

        z_prima = Z_prima[idx]
        z_dow = Z_dowell[idx]
        err = np.abs(z_prima - z_dow) / np.abs(z_dow) * 100

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {xi:>8.3f} {np.abs(z_prima):>14.6e} {np.abs(z_dow):>14.6e} {err:>9.2f}%")

    # 実部・虚部の比較
    print(f"\n実部・虚部の詳細:")
    print(f"{'Freq':>12} {'Re(PRIMA)':>14} {'Re(Dow)':>14} {'Im(PRIMA)':>14} {'Im(Dow)':>14}")
    print("-"*75)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_prima = Z_prima[idx]
        z_dow = Z_dowell[idx]

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {z_prima.real:>14.6e} {z_dow.real:>14.6e} {z_prima.imag:>14.6e} {z_dow.imag:>14.6e}")

    # DC極限での確認
    print(f"\nDC極限での確認:")
    print(f"  Z_dowell(10Hz) = {Z_dowell[0]:.6e}")
    print(f"  Z_prima(10Hz) = {Z_prima[0]:.6e}")
    print(f"  R_dc = {R_dc:.6e}")
    print(f"  jw*L_int_dc(10Hz) = {1j * 2*np.pi*10 * L_int_dc:.6e}")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: 振幅
    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z_prima), 'b-', linewidth=2, label='PRIMA')
    ax1.loglog(freqs, np.abs(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot 2: 誤差
    ax2 = axes[0, 1]
    err = np.abs(Z_prima - Z_dowell) / np.abs(Z_dowell) * 100
    ax2.semilogx(freqs, err, 'b-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Relative Error [%]')
    ax2.set_title('PRIMA vs Dowell Error')
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: 実部
    ax3 = axes[1, 0]
    ax3.loglog(freqs, np.real(Z_prima), 'b-', linewidth=2, label='PRIMA')
    ax3.loglog(freqs, np.real(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax3.axhline(y=R_dc, color='k', linestyle=':', label='R_dc')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z)')
    ax3.set_title('Real Part (Resistance)')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Plot 4: 虚部
    ax4 = axes[1, 1]
    ax4.loglog(freqs, np.imag(Z_prima), 'b-', linewidth=2, label='PRIMA')
    ax4.loglog(freqs, np.imag(Z_dowell), 'r--', linewidth=2, label='Dowell')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Im(Z)')
    ax4.set_title('Imaginary Part (Reactance)')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('verify_prima_dowell_correct.png', dpi=150)
    print(f"\nSaved: verify_prima_dowell_correct.png")
    plt.close()


if __name__ == '__main__':
    main()
