"""
Dowell表皮効果モデルのSPICE化

=== 2つのアプローチ ===

1. DC PRIMA + F_R/F_L補正 (推奨、Dowell式と完全一致)
   - DC値: R_dc = 1/(sigma*d), L_dc = mu*d/3
   - 周波数補正: Z(s) = R_dc * F_R(xi) + s * L_dc * F_L(xi)
   - SPICEでは周波数依存素子として実装
   - 誤差: < 1e-10% (数値精度限界)

2. RLラダー回路 (1D拡散固有モード、Dowellとは異なるモデル)
   - PRIMAラダー: L_n = mu*d/(4n-3), R_n = (4n-5)*4/(sigma*d)
   - 1D拡散方程式: Z = s*mu*(2/k/d)*tan(k*d/2)*d
   - Dowellのz*coth(z)とは異なる (tan vs coth)
   - 高周波極限では類似するが、低〜中周波で差異あり

=== 結論 ===
Dowell式をSPICE化するには、方法1(DC PRIMA + F_R/F_L補正)を使用。
方法2(RLラダー)はtan型の1D拡散モデルであり、Dowellのcoth型とは異なる。

=== SPICE実装 ===
方法1をSPICE回路で実装するには:

.SUBCKT DOWELL_SKIN in out d={d} sigma={sigma} mu={mu}
* 周波数依存R, L素子
* R(f) = R_dc * F_R(xi), xi = d / delta(f)
* L(f) = L_dc * F_L(xi)
R1 in out R='1/(sigma*d) * FR(d/delta(frequency))'
L1 in out L='mu*d/3 * FL(d/delta(frequency))'
.ENDS

ただし、標準SPICEはF_R, F_L関数をサポートしないため、
Verilog-AまたはPythonコシミュレーションが必要。
"""

import numpy as np
import matplotlib.pyplot as plt


MU_0 = 4 * np.pi * 1e-7


def calc_skin_depth(freq, sigma, mu=MU_0):
    """表皮深さ delta = sqrt(2 / (omega * mu * sigma))"""
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)
    return np.sqrt(2.0 / (omega * mu * sigma))


def dowell_F_R(xi):
    """Dowell F_R = xi * (sinh(xi) + sin(xi)) / (cosh(xi) + cos(xi))"""
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        return xi * (np.sinh(xi) + np.sin(xi)) / (np.cosh(xi) + np.cos(xi))
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.001
        xi_m = xi[mask]
        result[mask] = xi_m * (np.sinh(xi_m) + np.sin(xi_m)) / (np.cosh(xi_m) + np.cos(xi_m))
        return result


def dowell_F_L(xi):
    """Dowell F_L = 3*(sinh(2*xi) - sin(2*xi)) / (xi^3 * (cosh(2*xi) - cos(2*xi)))"""
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        xi2 = 2 * xi
        denom = xi**3 * (np.cosh(xi2) - np.cos(xi2))
        if abs(denom) < 1e-30:
            return 1.0
        return 3.0 * (np.sinh(xi2) - np.sin(xi2)) / denom
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.001
        xi_m = xi[mask]
        xi2 = 2 * xi_m
        denom = xi_m**3 * (np.cosh(xi2) - np.cos(xi2))
        result[mask] = 3.0 * (np.sinh(xi2) - np.sin(xi2)) / denom
        return result


def dowell_impedance(R_dc, L_dc, d, sigma, mu, freq):
    """Dowell式: Z = R_dc * F_R + j*omega * L_dc * F_L"""
    omega = 2 * np.pi * freq
    delta = calc_skin_depth(freq, sigma, mu)
    xi = d / delta
    F_R = dowell_F_R(xi)
    F_L = dowell_F_L(xi)
    return R_dc * F_R + 1j * omega * L_dc * F_L


def diffusion_impedance(d, sigma, mu, freq):
    """1D拡散: Z = s*mu*(2/k/d)*tan(k*d/2)*d"""
    s = 1j * 2 * np.pi * freq

    if np.abs(s) < 1e-20:
        return s * mu * d / 3.0

    k = np.sqrt(s * mu * sigma)
    kd2 = k * d / 2.0

    if np.abs(kd2) < 0.01:
        return s * mu * d / 3.0 * (1 + kd2**2 / 3.0)

    return s * mu * (2.0 / k / d) * np.tan(kd2) * d


def main():
    print("=" * 70)
    print("Dowell vs 1D Diffusion: 表皮効果モデル比較")
    print("=" * 70)

    # Conductor parameters
    d = 0.1e-3       # 0.1 mm thickness
    sigma = 5.8e7    # Copper conductivity [S/m]
    mu = MU_0        # Permeability

    # DC parameters
    R_dc = 1.0 / (sigma * d)
    L_dc = mu * d / 3.0

    print(f"\nConductor parameters:")
    print(f"  d     = {d*1e3:.3f} mm")
    print(f"  sigma = {sigma:.2e} S/m")
    print(f"  mu    = {mu:.3e} H/m")
    print(f"\nDC parameters:")
    print(f"  R_dc = 1/(sigma*d) = {R_dc:.6e} Ohm*m^2")
    print(f"  L_dc = mu*d/3 = {L_dc:.6e} H*m^2")

    # Frequency response
    frequencies = np.logspace(0, 8, 200)

    Z_dowell = np.array([dowell_impedance(R_dc, L_dc, d, sigma, mu, f) for f in frequencies])
    Z_diff = np.array([diffusion_impedance(d, sigma, mu, f) for f in frequencies])

    # Comparison
    print(f"\n--- Dowell vs 1D Diffusion Comparison ---")
    print(f"{'Frequency':>12} {'|Z_Dowell|':>14} {'|Z_Diff|':>14} {'Diff %':>12}")
    print("-" * 55)

    test_freqs = [1, 100, 1e3, 10e3, 100e3, 1e6, 10e6, 100e6]
    for f in test_freqs:
        idx = np.argmin(np.abs(frequencies - f))
        z_dow = Z_dowell[idx]
        z_dif = Z_diff[idx]
        diff = (np.abs(z_dif) - np.abs(z_dow)) / np.abs(z_dow) * 100

        freq_str = f"{f:.0f} Hz" if f < 1000 else \
                   f"{f/1e3:.0f} kHz" if f < 1e6 else \
                   f"{f/1e6:.0f} MHz"
        print(f"{freq_str:>12} {np.abs(z_dow):>14.6e} {np.abs(z_dif):>14.6e} {diff:>11.2f}%")

    # Skin depth info
    print(f"\n--- Skin Depth Reference ---")
    for f in test_freqs:
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        freq_str = f"{f:.0f} Hz" if f < 1000 else \
                   f"{f/1e3:.0f} kHz" if f < 1e6 else \
                   f"{f/1e6:.0f} MHz"
        print(f"  {freq_str:>10}: delta = {delta*1e6:.2f} um, xi = d/delta = {xi:.3f}")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Magnitude comparison
    ax1 = axes[0, 0]
    ax1.loglog(frequencies, np.abs(Z_dowell), 'r-', linewidth=2, label='Dowell: z*coth(z)')
    ax1.loglog(frequencies, np.abs(Z_diff), 'b--', linewidth=2, label='1D Diff: tan(kd/2)')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm*m^2]')
    ax1.set_title('Surface Impedance: Dowell vs 1D Diffusion')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Relative difference
    ax2 = axes[0, 1]
    rel_diff = (np.abs(Z_diff) - np.abs(Z_dowell)) / np.abs(Z_dowell) * 100
    ax2.semilogx(frequencies, rel_diff, 'b-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('(|Z_diff| - |Z_Dowell|) / |Z_Dowell| [%]')
    ax2.set_title('1D Diffusion vs Dowell Difference')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

    # Real part (resistance)
    ax3 = axes[1, 0]
    ax3.loglog(frequencies, np.real(Z_dowell), 'r-', linewidth=2, label='Dowell')
    ax3.loglog(frequencies, np.real(Z_diff), 'b--', linewidth=2, label='1D Diff')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Re(Z) [Ohm*m^2]')
    ax3.set_title('AC Resistance')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Imaginary part (reactance)
    ax4 = axes[1, 1]
    ax4.loglog(frequencies, np.abs(np.imag(Z_dowell)), 'r-', linewidth=2, label='Dowell')
    ax4.loglog(frequencies, np.abs(np.imag(Z_diff)), 'b--', linewidth=2, label='1D Diff')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('|Im(Z)| [Ohm*m^2]')
    ax4.set_title('AC Reactance')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_dowell_spice.png', dpi=150)
    print(f"\nSaved: demo_dowell_spice.png")
    plt.close()

    # SPICE output for DC + F_R/F_L (frequency-dependent SPICE model)
    print(f"\n--- SPICE Model (DC + F_R/F_L) ---")
    print(f"""
* Dowell Surface Impedance Model
* Use with frequency-dependent element capability (Verilog-A or Python)
*
* Parameters:
*   d = {d} m
*   sigma = {sigma} S/m
*   mu = {mu} H/m
*   R_dc = {R_dc:.6e} Ohm*m^2
*   L_dc = {L_dc:.6e} H*m^2
*
* Frequency-dependent formulas:
*   xi = d / delta(f) = d * sqrt(pi*f*mu*sigma)
*   F_R(xi) = xi * (sinh(xi) + sin(xi)) / (cosh(xi) + cos(xi))
*   F_L(xi) = 3 * (sinh(2*xi) - sin(2*xi)) / (xi^3 * (cosh(2*xi) - cos(2*xi)))
*   R(f) = R_dc * F_R(xi)
*   L(f) = L_dc * F_L(xi)
*
* For fixed-element SPICE, use lookup table at discrete frequencies.
""")

    # Create frequency lookup table
    print("* Lookup table for fixed-element SPICE:")
    print("* Freq [Hz]     R [Ohm*m^2]    L [H*m^2]")
    for f in [1, 10, 100, 1e3, 10e3, 100e3, 1e6, 10e6, 100e6]:
        delta = calc_skin_depth(f, sigma, mu)
        xi = d / delta
        R_f = R_dc * dowell_F_R(xi)
        L_f = L_dc * dowell_F_L(xi)
        freq_str = f"{f:.0e}"
        print(f"* {freq_str:>10}  {R_f:.6e}  {L_f:.6e}")

    print(f"\n" + "=" * 70)
    print("結論:")
    print("  1. Dowell式 (z*coth(z)) と 1D拡散 (tan(kd/2)) は異なるモデル")
    print("  2. 低周波では両者の差が大きい (1D拡散はDCで Im(Z)=0)")
    print("  3. 高周波 (xi >> 1) では両者とも Z ~ (1+j)*R_dc*xi/sqrt(2) に収束")
    print("  4. SPICE実装にはDowell式 (DC + F_R/F_L補正) を推奨")
    print("=" * 70)


if __name__ == '__main__':
    main()
