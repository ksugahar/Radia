"""
PEEC + Dowell vs PRIMA + F_R/F_L の一致検証

PEEC with Dowell:
  Z_matrix(s) = R_dc * F_R(s) + s * L_matrix * F_L(s)
  (内部インダクタンスのみF_Lで変化、外部インダクタンスは変化しない)

PRIMA with F_R/F_L:
  R_prima = R_dc として、F_R(s)で周波数依存性を表現
  L_prima を F_L(s) で変化させる

単一導体の場合で一致を確認。
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


# =============================================
# PEEC + Dowell (周波数ごとにF_R, F_Lを適用)
# =============================================
def calc_peec_dowell_impedance(R_dc, L_int_dc, L_ext, d, sigma, mu, freqs):
    """
    PEEC with Dowell correction

    Z(s) = R_dc * F_R(xi) + s * (L_int_dc * F_L(xi) + L_ext)

    Note: L_ext (外部インダクタンス) は表皮効果の影響を受けない
    """
    Z = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2 * np.pi * f
        s = 1j * omega
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        F_R = dowell_F_R(xi)
        F_L = dowell_F_L(xi)

        # 内部インダクタンスのみF_Lで変化
        L_total = L_int_dc * F_L + L_ext
        Z[i] = R_dc * F_R + s * L_total

    return Z


# =============================================
# 単純Dowell式 (内部インピーダンスのみ、外部Lなし)
# =============================================
def calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs):
    """
    Pure Dowell (internal impedance only)
    Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)
    """
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


# =============================================
# PRIMA I form (周波数非依存のR, L)
# =============================================
def build_tridiagonal(L):
    n = len(L)
    U = np.eye(n)
    for i in range(n - 1):
        U[i, i + 1] = -1
    return U.T @ np.diag(L) @ U


def calc_prima_impedance(R_prima, L_prima, freqs):
    """PRIMA I impedance: (R + sL)I = V, Z = 1/I[0]"""
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


def skin_effect_prima_params(d, sigma, mu=MU_0, n_stages=7):
    """元のPRIMA I パラメータ (1D拡散式由来)"""
    R_prima = np.zeros(n_stages)
    L_prima = np.zeros(n_stages)

    for n in range(1, n_stages + 1):
        if n == 1:
            R_prima[n-1] = 1e-16
        else:
            R_prima[n-1] = (4*n - 5) * 4.0 / (sigma * d)
        L_prima[n-1] = d * mu / (4*n - 3)

    return R_prima, L_prima


def main():
    print("="*70)
    print("PEEC + Dowell vs PRIMA 一致検証")
    print("="*70)

    d = 0.1e-3      # 導体厚さ
    sigma = 5.8e7   # 銅
    mu = MU_0

    # Dowell DC parameters
    R_dc = 1.0 / (sigma * d)
    L_int_dc = mu * d / 3.0

    # 外部インダクタンス (仮に内部の10倍)
    L_ext = L_int_dc * 10

    print(f"\n物理パラメータ:")
    print(f"  d = {d*1e6:.1f} um")
    print(f"  sigma = {sigma:.2e} S/m")
    print(f"  R_dc = {R_dc:.6e} Ohm*m^2")
    print(f"  L_int_dc = {L_int_dc:.6e} H*m^2")
    print(f"  L_ext = {L_ext:.6e} H*m^2")

    freqs = np.logspace(0, 8, 200)

    # PEEC + Dowell
    Z_peec_dowell = calc_peec_dowell_impedance(R_dc, L_int_dc, L_ext, d, sigma, mu, freqs)

    # Pure Dowell (no external L)
    Z_dowell_only = calc_dowell_impedance(R_dc, L_int_dc, d, sigma, mu, freqs)

    # PRIMA (元の1D拡散式由来)
    n_stages = 10
    R_prima, L_prima = skin_effect_prima_params(d, sigma, mu, n_stages)
    Z_prima = calc_prima_impedance(R_prima, L_prima, freqs)

    print(f"\n--- Dowell式のみ (外部Lなし) vs PRIMA ---")
    print(f"{'Freq':>12} {'|Z_Dowell|':>14} {'|Z_PRIMA|':>14} {'Err':>10}")
    print("-"*55)

    test_freqs = [1, 10, 100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        z_dow = Z_dowell_only[idx]
        z_prima = Z_prima[idx]
        err = np.abs(z_prima - z_dow) / np.abs(z_dow) * 100
        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {np.abs(z_dow):>14.6e} {np.abs(z_prima):>14.6e} {err:>9.2f}%")

    # F_R, F_L の比較
    print(f"\n--- F_R, F_L 比較 ---")
    print(f"PRIMAからF_R, F_Lを逆算: F_R = Re(Z)/R_dc, F_L = Im(Z)/(omega*L_int_dc)")
    print(f"{'Freq':>12} {'xi':>8} {'F_R(PRIMA)':>12} {'F_R(Dow)':>12} {'F_L(PRIMA)':>12} {'F_L(Dow)':>12}")
    print("-"*75)

    for f in test_freqs:
        idx = np.argmin(np.abs(freqs - f))
        omega = 2 * np.pi * f
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta

        z_prima_val = Z_prima[idx]
        F_R_prima = z_prima_val.real / R_dc if R_dc > 0 else 1.0
        F_L_prima = z_prima_val.imag / (omega * L_int_dc) if omega * L_int_dc > 1e-30 else 1.0

        F_R_dow = dowell_F_R(xi)
        F_L_dow = dowell_F_L(xi)

        freq_str = f"{f:.0f} Hz" if f < 1000 else f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10} {xi:>8.4f} {F_R_prima:>12.4f} {F_R_dow:>12.4f} {F_L_prima:>12.4f} {F_L_dow:>12.4f}")

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z_dowell_only), 'r-', linewidth=2, label='Dowell')
    ax1.loglog(freqs, np.abs(Z_prima), 'b--', linewidth=2, label='PRIMA (1D diff)')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Internal Impedance: Dowell vs PRIMA')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    ax2 = axes[0, 1]
    # F_R, F_L vs xi
    xi_arr = np.array([d / calc_skin_depth(f, sigma, mu) for f in freqs])
    F_R_prima_arr = np.real(Z_prima) / R_dc
    F_L_prima_arr = np.imag(Z_prima) / (2 * np.pi * freqs * L_int_dc)
    F_R_dow_arr = dowell_F_R(xi_arr)
    F_L_dow_arr = dowell_F_L(xi_arr)

    ax2.plot(xi_arr, F_R_prima_arr, 'b-', linewidth=2, label='F_R (PRIMA)')
    ax2.plot(xi_arr, F_R_dow_arr, 'r--', linewidth=2, label='F_R (Dowell)')
    ax2.set_xlabel('xi = d/delta')
    ax2.set_ylabel('F_R')
    ax2.set_title('F_R vs xi')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, 10])

    ax3 = axes[1, 0]
    ax3.plot(xi_arr, F_L_prima_arr, 'b-', linewidth=2, label='F_L (PRIMA)')
    ax3.plot(xi_arr, F_L_dow_arr, 'r--', linewidth=2, label='F_L (Dowell)')
    ax3.set_xlabel('xi = d/delta')
    ax3.set_ylabel('F_L')
    ax3.set_title('F_L vs xi')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, 10])

    ax4 = axes[1, 1]
    err = np.abs(Z_prima - Z_dowell_only) / np.abs(Z_dowell_only) * 100
    ax4.semilogx(freqs, err, 'b-', linewidth=2)
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Relative Error [%]')
    ax4.set_title('PRIMA vs Dowell Error')
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('verify_peec_prima_dowell.png', dpi=150)
    print(f"\nSaved: verify_peec_prima_dowell.png")
    plt.close()


if __name__ == '__main__':
    main()
