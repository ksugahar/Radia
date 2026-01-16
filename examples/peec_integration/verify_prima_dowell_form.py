"""
PRIMAインピーダンスをDowell形式に変換して比較

PRIMA I のインピーダンス Z(s) を計算し、それを
  Z(s) = R_dc * F_R(s) + j*omega * L_int_dc * F_L(s)
の形式に分解して、F_R と F_L を抽出し、Dowell の F_R, F_L と比較する。
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
    """PRIMA I パラメータ（1D拡散方程式の連分数展開）"""
    R_prima = np.zeros(n_stages)
    L_prima = np.zeros(n_stages)

    for n in range(1, n_stages + 1):
        if n == 1:
            R_prima[n-1] = 1e-16  # ほぼゼロ
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
    """PRIMA I インピーダンスを計算"""
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
    """Dowell の F_R"""
    if xi < 0.01:
        return 1.0
    sinh_xi = np.sinh(xi)
    cosh_xi = np.cosh(xi)
    sin_xi = np.sin(xi)
    cos_xi = np.cos(xi)
    M = (sinh_xi + sin_xi) / (cosh_xi + cos_xi)
    return xi * M


def dowell_F_L(xi):
    """Dowell の F_L"""
    if xi < 0.01:
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


def main():
    print("="*70)
    print("PRIMA -> Dowell形式変換による比較")
    print("="*70)

    # 物理パラメータ
    d = 0.1e-3      # 0.1 mm
    sigma = 5.8e7   # 銅
    mu = MU_0

    # DC パラメータ
    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    print(f"\n物理パラメータ:")
    print(f"  d = {d*1e3:.3f} mm")
    print(f"  sigma = {sigma:.2e} S/m")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = mu*d/3 = {L_int_dc:.6e} H*m^2")

    # PRIMA パラメータ
    n_stages = 10
    R_prima, L_prima = skin_effect_prima_params(d, sigma, mu, n_stages)

    print(f"\nPRIMA I パラメータ ({n_stages}段):")
    print(f"  R_prima[0:3] = {R_prima[:3]}")
    print(f"  L_prima[0:3] = {L_prima[:3]}")

    # 周波数範囲
    freqs = np.logspace(1, 8, 200)

    # PRIMA インピーダンス計算
    Z_prima = calc_prima_impedance(R_prima, L_prima, freqs)

    # PRIMAから F_R, F_L を抽出
    # Z(f) = R_dc * F_R + j*omega * L_int_dc * F_L
    # Re(Z) = R_dc * F_R  =>  F_R = Re(Z) / R_dc
    # Im(Z) = omega * L_int_dc * F_L  =>  F_L = Im(Z) / (omega * L_int_dc)

    F_R_prima = np.real(Z_prima) / R_dc
    F_L_prima = np.zeros(len(freqs))
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        if omega * L_int_dc > 1e-30:
            F_L_prima[i] = np.imag(Z_prima[i]) / (omega * L_int_dc)
        else:
            F_L_prima[i] = 1.0

    # Dowell の F_R, F_L
    F_R_dowell = np.zeros(len(freqs))
    F_L_dowell = np.zeros(len(freqs))
    xi_arr = np.zeros(len(freqs))

    for i, f in enumerate(freqs):
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        xi_arr[i] = xi
        F_R_dowell[i] = dowell_F_R(xi)
        F_L_dowell[i] = dowell_F_L(xi)

    # 比較表示
    print(f"\n{'Freq':>12} {'xi':>8} {'F_R(PRIMA)':>12} {'F_R(Dow)':>12} {'F_L(PRIMA)':>12} {'F_L(Dow)':>12}")
    print("-"*75)

    test_freqs = [100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        xi = xi_arr[idx]

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {xi:>8.3f} {F_R_prima[idx]:>12.6f} {F_R_dowell[idx]:>12.6f} "
              f"{F_L_prima[idx]:>12.6f} {F_L_dowell[idx]:>12.6f}")

    # F_R, F_L の誤差
    err_F_R = np.abs(F_R_prima - F_R_dowell) / np.maximum(np.abs(F_R_dowell), 1e-10) * 100
    err_F_L = np.abs(F_L_prima - F_L_dowell) / np.maximum(np.abs(F_L_dowell), 1e-10) * 100

    print(f"\n誤差 (PRIMA vs Dowell):")
    print(f"{'Freq':>12} {'xi':>8} {'Err F_R':>12} {'Err F_L':>12}")
    print("-"*50)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        xi = xi_arr[idx]

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1e3:
            freq_str = f"{f/1e3:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"  {freq_str:>10} {xi:>8.3f} {err_F_R[idx]:>11.4f}% {err_F_L[idx]:>11.4f}%")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: F_R 比較
    ax1 = axes[0, 0]
    ax1.semilogx(freqs, F_R_prima, 'b-', linewidth=2, label='F_R (PRIMA)')
    ax1.semilogx(freqs, F_R_dowell, 'r--', linewidth=2, label='F_R (Dowell)')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('F_R')
    ax1.set_title('Resistance Factor F_R')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot 2: F_L 比較
    ax2 = axes[0, 1]
    ax2.semilogx(freqs, F_L_prima, 'b-', linewidth=2, label='F_L (PRIMA)')
    ax2.semilogx(freqs, F_L_dowell, 'r--', linewidth=2, label='F_L (Dowell)')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('F_L')
    ax2.set_title('Inductance Factor F_L')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: F_R vs xi
    ax3 = axes[1, 0]
    ax3.plot(xi_arr, F_R_prima, 'b-', linewidth=2, label='F_R (PRIMA)')
    ax3.plot(xi_arr, F_R_dowell, 'r--', linewidth=2, label='F_R (Dowell)')
    ax3.set_xlabel('xi = d/delta')
    ax3.set_ylabel('F_R')
    ax3.set_title('F_R vs xi')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, 10])

    # Plot 4: F_L vs xi
    ax4 = axes[1, 1]
    ax4.plot(xi_arr, F_L_prima, 'b-', linewidth=2, label='F_L (PRIMA)')
    ax4.plot(xi_arr, F_L_dowell, 'r--', linewidth=2, label='F_L (Dowell)')
    ax4.set_xlabel('xi = d/delta')
    ax4.set_ylabel('F_L')
    ax4.set_title('F_L vs xi')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, 10])

    plt.tight_layout()
    plt.savefig('verify_prima_dowell_form.png', dpi=150)
    print(f"\nSaved: verify_prima_dowell_form.png")
    plt.close()


if __name__ == '__main__':
    main()
