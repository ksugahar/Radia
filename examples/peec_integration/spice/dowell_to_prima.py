"""
Dowell式をPRIMA形式で表現

Dowell式:
  Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)
  xi = d/delta = d * sqrt(omega*mu*sigma/2)

F_R(xi)とF_L(xi)をsの関数として展開し、PRIMA I形式で表現できるか検証。
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


def dowell_F_R(xi):
    """F_R = xi * (sinh(xi) + sin(xi)) / (cosh(xi) + cos(xi))"""
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0 + xi**4 / 45.0
        sinh_xi = np.sinh(xi)
        cosh_xi = np.cosh(xi)
        sin_xi = np.sin(xi)
        cos_xi = np.cos(xi)
        M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
        return xi * M
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.001
        xi_m = xi[mask]
        sinh_xi = np.sinh(xi_m)
        cosh_xi = np.cosh(xi_m)
        sin_xi = np.sin(xi_m)
        cos_xi = np.cos(xi_m)
        M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
        result[mask] = xi_m * M
        result[~mask] = 1.0 + xi[~mask]**4 / 45.0
        return result


def dowell_F_L(xi):
    """F_L = 3*(sinh(2*xi) - sin(2*xi)) / (xi^3 * (cosh(2*xi) - cos(2*xi)))"""
    if np.isscalar(xi):
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
        result[~mask] = 1.0 - xi[~mask]**4 / 45.0
        return result


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


def dowell_prima_params(d, sigma, mu=MU_0, n_stages=7):
    """
    Dowell式に基づくPRIMAパラメータ

    Dowell DC parameters:
      R_dc = 1/(sigma*d)
      L_int_dc = mu*d/3

    PRIMA I DC極限で Z -> R_dc + j*omega*L_int_dc となるように設定
    """
    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    R_prima = np.zeros(n_stages)
    L_prima = np.zeros(n_stages)

    # Stage 1: DC抵抗とDC内部インダクタンス
    R_prima[0] = R_dc
    L_prima[0] = L_int_dc

    # Stage 2以降: 高周波補正項
    # 表皮効果による追加のR, L
    for n in range(2, n_stages + 1):
        # 経験的なスケーリング (要調整)
        R_prima[n-1] = R_dc * (4*n - 5)
        L_prima[n-1] = L_int_dc / (4*n - 3)

    return R_prima, L_prima


def main():
    print("="*70)
    print("Dowell式のPRIMA表現")
    print("="*70)

    d = 0.1e-3
    sigma = 5.8e7
    mu = MU_0

    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\nDowell DC parameters:")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = mu*d/3 = {L_int_dc:.6e} H*m^2")

    freqs = np.logspace(0, 8, 200)

    # Dowell式
    Z_dowell = calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs)

    # PRIMA (Dowell DC parameters)
    n_stages = 10
    R_prima, L_prima = dowell_prima_params(d, sigma, mu, n_stages)
    Z_prima = calc_prima_impedance(R_prima, L_prima, freqs)

    print(f"\nPRIMA parameters ({n_stages} stages):")
    print(f"  R_prima[:3] = {R_prima[:3]}")
    print(f"  L_prima[:3] = {L_prima[:3]}")

    # DC極限確認
    print(f"\nDC極限 (f=1Hz):")
    print(f"  Z_dowell = {Z_dowell[0]:.6e}")
    print(f"  Z_prima = {Z_prima[0]:.6e}")
    print(f"  R_dc + j*w*L_dc = {R_dc + 1j*2*np.pi*1*L_int_dc:.6e}")

    # 比較
    print(f"\n{'Freq':>12} {'|Z_Dowell|':>14} {'|Z_PRIMA|':>14} {'Err':>10}")
    print("-"*55)

    test_freqs = [1, 100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_dow = Z_dowell[idx]
        z_prima = Z_prima[idx]
        err = np.abs(z_prima - z_dow) / np.abs(z_dow) * 100
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {np.abs(z_dow):>14.6e} {np.abs(z_prima):>14.6e} {err:>9.2f}%")

    # PRIMAから逆算したF_R, F_L
    print(f"\nPRIMAから逆算したF_R, F_L:")
    print(f"{'Freq':>12} {'xi':>8} {'F_R(PRIMA)':>12} {'F_R(Dow)':>12} {'F_L(PRIMA)':>12} {'F_L(Dow)':>12}")
    print("-"*75)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        omega = 2 * np.pi * f
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta

        z_prima = Z_prima[idx]
        F_R_prima = z_prima.real / R_dc
        F_L_prima = z_prima.imag / (omega * L_int_dc) if omega * L_int_dc > 1e-30 else 1.0

        F_R_dow = dowell_F_R(xi)
        F_L_dow = dowell_F_L(xi)

        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {xi:>8.3f} {F_R_prima:>12.4f} {F_R_dow:>12.4f} {F_L_prima:>12.4f} {F_L_dow:>12.4f}")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z_dowell), 'r-', linewidth=2, label='Dowell')
    ax1.loglog(freqs, np.abs(Z_prima), 'b--', linewidth=2, label='PRIMA')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    ax2 = axes[0, 1]
    err = np.abs(Z_prima - Z_dowell) / np.abs(Z_dowell) * 100
    ax2.semilogx(freqs, err, 'b-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Relative Error [%]')
    ax2.set_title('PRIMA vs Dowell Error')
    ax2.grid(True, which='both', alpha=0.3)

    ax3 = axes[1, 0]
    ax3.loglog(freqs, np.real(Z_dowell), 'r-', linewidth=2, label='Dowell')
    ax3.loglog(freqs, np.real(Z_prima), 'b--', linewidth=2, label='PRIMA')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z)')
    ax3.set_title('Real Part (Resistance)')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    ax4 = axes[1, 1]
    ax4.loglog(freqs, np.imag(Z_dowell), 'r-', linewidth=2, label='Dowell')
    ax4.loglog(freqs, np.imag(Z_prima), 'b--', linewidth=2, label='PRIMA')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Im(Z)')
    ax4.set_title('Imaginary Part (Reactance)')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('dowell_to_prima.png', dpi=150)
    print(f"\nSaved: dowell_to_prima.png")
    plt.close()


if __name__ == '__main__':
    main()
