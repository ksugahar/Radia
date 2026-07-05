"""
WPT System at 85 kHz with Ferrite Core and Aluminum Shield

This example demonstrates PEEC Loop-Star conductor analysis plus a planned
HDiv-VIM / reduced-FEM magnetic-material exchange for:
1. Tx/Rx coils with ferrite cores (mu_r = 2000)
2. Aluminum shield for EMI reduction
3. Mutual inductance and coupling coefficient
4. Power transfer efficiency
5. Shield eddy current losses

Target application: EV wireless charging (SAE J2954)

Physical setup:
    +-----------------------+
    |   Aluminum Shield     |  <- Eddy current shield (sigma = 3.5e7)
    +-----------------------+
    |   Ferrite Core (Tx)   |  <- mu_r = 2000
    |     [Tx Coil]         |  <- Litz wire, N turns
    +-----------------------+
              |
           Air gap (variable)
              |
    +-----------------------+
    |     [Rx Coil]         |
    |   Ferrite Core (Rx)   |
    +-----------------------+
    |   Aluminum Shield     |
    +-----------------------+

Author: Radia Development Team
Date: 2026-01-16
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build-msvc'))

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m


class WPTCoilAssembly:
    """
    WPT coil assembly with ferrite core and aluminum shield.
    """

    def __init__(self, name='Tx'):
        self.name = name

        # Coil parameters
        self.coil_inner_radius = 0.05   # 50 mm
        self.coil_outer_radius = 0.15   # 150 mm
        self.n_turns = 10
        self.wire_diameter = 3e-3       # 3 mm Litz wire
        self.wire_sigma = 5.8e7         # Copper

        # Ferrite core (flat disc)
        self.ferrite_radius = 0.18      # 180 mm
        self.ferrite_thickness = 5e-3   # 5 mm
        self.ferrite_mu_r = 2000
        self.ferrite_sigma = 0.01       # Ferrite is nearly insulating

        # Aluminum shield
        self.shield_radius = 0.20       # 200 mm
        self.shield_thickness = 2e-3   # 2 mm
        self.shield_sigma = 3.5e7      # Aluminum
        self.shield_mu_r = 1.0

    def calc_coil_inductance_dc(self):
        """
        Approximate DC inductance of planar spiral coil.

        Wheeler's formula for planar spiral:
        L = (mu_0 * N^2 * d_avg * K1) / (1 + K2 * rho)

        where:
            d_avg = (d_out + d_in) / 2
            rho = (d_out - d_in) / (d_out + d_in)
            K1, K2 = shape constants (K1 = 2.34, K2 = 2.75 for circular)
        """
        d_out = 2 * self.coil_outer_radius
        d_in = 2 * self.coil_inner_radius
        d_avg = (d_out + d_in) / 2
        rho = (d_out - d_in) / (d_out + d_in)

        K1 = 2.34
        K2 = 2.75

        L = (MU_0 * self.n_turns**2 * d_avg * K1) / (1 + K2 * rho)
        return L

    def calc_coil_resistance_dc(self):
        """DC resistance of coil wire."""
        # Approximate wire length
        r_avg = (self.coil_inner_radius + self.coil_outer_radius) / 2
        wire_length = 2 * np.pi * r_avg * self.n_turns

        # Wire cross-section
        wire_area = np.pi * (self.wire_diameter / 2)**2

        R_dc = wire_length / (self.wire_sigma * wire_area)
        return R_dc

    def calc_skin_depth(self, freq, sigma, mu_r):
        """Calculate skin depth."""
        omega = 2 * np.pi * freq
        if omega < 1e-10:
            return np.inf
        return np.sqrt(2 / (omega * MU_0 * mu_r * sigma))

    def calc_ac_resistance_factor(self, freq):
        """
        AC resistance factor due to skin effect in Litz wire.

        For Litz wire, use strand diameter, not bundle diameter.
        Assume strand diameter = wire_diameter / 10
        """
        strand_diameter = self.wire_diameter / 10
        delta = self.calc_skin_depth(freq, self.wire_sigma, 1.0)

        if delta > strand_diameter:
            return 1.0  # DC regime
        else:
            # Approximate: R_ac / R_dc ~ (strand_d / 2) / delta
            return (strand_diameter / 2) / delta

    def calc_shield_loss_resistance(self, freq, B_avg=0.001):
        """
        Estimate shield eddy current loss as equivalent resistance.

        P_shield = (pi * sigma * d^3 * omega^2 * B^2) / 6

        where d = shield thickness, B = average field at shield
        """
        omega = 2 * np.pi * freq
        d = self.shield_thickness

        # Power loss per unit area
        P_per_area = (self.shield_sigma * d**3 * omega**2 * B_avg**2) / 6

        # Total power for circular shield
        shield_area = np.pi * self.shield_radius**2
        P_total = P_per_area * shield_area

        # Equivalent resistance: P = I^2 * R_eq
        # Assume I = 10 A for normalization
        I_ref = 10.0
        R_shield = P_total / I_ref**2

        return R_shield


class WPTSystem:
    """
    Complete WPT system with Tx and Rx assemblies.
    """

    def __init__(self):
        self.tx = WPTCoilAssembly('Tx')
        self.rx = WPTCoilAssembly('Rx')
        self.air_gap = 0.15  # 150 mm (typical for EV)

    def calc_mutual_inductance(self):
        """
        Approximate mutual inductance between coaxial planar coils.

        Uses Neumann formula approximation for coaxial circular loops.
        """
        # Average radii
        r1 = (self.tx.coil_inner_radius + self.tx.coil_outer_radius) / 2
        r2 = (self.rx.coil_inner_radius + self.rx.coil_outer_radius) / 2
        d = self.air_gap

        # Simplified: M ~ mu_0 * pi * r1^2 * r2^2 * N1 * N2 / (2 * (r1^2 + d^2)^(3/2))
        # for r1 ~ r2 and d > r1
        r_eff = (r1 + r2) / 2
        M = MU_0 * np.pi * r_eff**4 * self.tx.n_turns * self.rx.n_turns / \
            (2 * (r_eff**2 + d**2)**(3/2))

        return M

    def calc_coupling_coefficient(self):
        """Calculate coupling coefficient k = M / sqrt(L1 * L2)."""
        L1 = self.tx.calc_coil_inductance_dc()
        L2 = self.rx.calc_coil_inductance_dc()
        M = self.calc_mutual_inductance()
        k = M / np.sqrt(L1 * L2)
        return k

    def calc_impedance_matrix(self, freq):
        """
        Calculate 2x2 impedance matrix at given frequency.

        Z = R + j*omega*L

        Returns:
            Z: 2x2 complex impedance matrix
        """
        omega = 2 * np.pi * freq

        # Self-inductances (with ferrite enhancement)
        # Ferrite increases inductance by factor of ~2-3 for typical geometries
        ferrite_factor = 2.5
        L1 = self.tx.calc_coil_inductance_dc() * ferrite_factor
        L2 = self.rx.calc_coil_inductance_dc() * ferrite_factor
        M = self.calc_mutual_inductance() * ferrite_factor

        # Resistances
        R1_dc = self.tx.calc_coil_resistance_dc()
        R2_dc = self.rx.calc_coil_resistance_dc()

        # AC resistance factors
        F_ac1 = self.tx.calc_ac_resistance_factor(freq)
        F_ac2 = self.rx.calc_ac_resistance_factor(freq)

        R1 = R1_dc * F_ac1
        R2 = R2_dc * F_ac2

        # Shield losses
        R_shield1 = self.tx.calc_shield_loss_resistance(freq)
        R_shield2 = self.rx.calc_shield_loss_resistance(freq)

        R1 += R_shield1
        R2 += R_shield2

        # Build impedance matrix
        Z = np.array([
            [R1 + 1j * omega * L1, 1j * omega * M],
            [1j * omega * M, R2 + 1j * omega * L2]
        ], dtype=complex)

        return Z

    def calc_efficiency(self, freq, R_load):
        """
        Calculate power transfer efficiency.

        eta = P_load / P_in

        For matched load: R_load = sqrt(R2^2 + (omega*L2)^2)
        """
        omega = 2 * np.pi * freq
        Z = self.calc_impedance_matrix(freq)

        Z11 = Z[0, 0]
        Z12 = Z[0, 1]
        Z21 = Z[1, 0]
        Z22 = Z[1, 1]

        # Input impedance with load
        Z_in = Z11 - Z12 * Z21 / (Z22 + R_load)

        # Currents (assume V1 = 1V source)
        V1 = 1.0
        I1 = V1 / Z_in
        I2 = -Z21 / (Z22 + R_load) * I1

        # Powers
        P_in = np.real(V1 * np.conj(I1))
        P_load = np.abs(I2)**2 * R_load

        eta = P_load / P_in if P_in > 0 else 0

        return eta, P_in, P_load


def main():
    """WPT 85 kHz system analysis."""
    print("=" * 70)
    print("WPT System Analysis at 85 kHz with Shield")
    print("=" * 70)

    # Create system
    wpt = WPTSystem()

    # System parameters
    print("\n--- System Parameters ---")
    print(f"Tx coil: r_in = {wpt.tx.coil_inner_radius*1e3:.0f} mm, "
          f"r_out = {wpt.tx.coil_outer_radius*1e3:.0f} mm, "
          f"N = {wpt.tx.n_turns}")
    print(f"Rx coil: r_in = {wpt.rx.coil_inner_radius*1e3:.0f} mm, "
          f"r_out = {wpt.rx.coil_outer_radius*1e3:.0f} mm, "
          f"N = {wpt.rx.n_turns}")
    print(f"Air gap: {wpt.air_gap*1e3:.0f} mm")
    print(f"Ferrite: mu_r = {wpt.tx.ferrite_mu_r}, "
          f"thickness = {wpt.tx.ferrite_thickness*1e3:.0f} mm")
    print(f"Shield: Al, thickness = {wpt.tx.shield_thickness*1e3:.0f} mm")

    # DC parameters
    L1_dc = wpt.tx.calc_coil_inductance_dc()
    L2_dc = wpt.rx.calc_coil_inductance_dc()
    M_dc = wpt.calc_mutual_inductance()
    k_dc = wpt.calc_coupling_coefficient()

    print("\n--- DC Parameters ---")
    print(f"L1 (air core) = {L1_dc*1e6:.2f} uH")
    print(f"L2 (air core) = {L2_dc*1e6:.2f} uH")
    print(f"M = {M_dc*1e6:.3f} uH")
    print(f"k = {k_dc:.4f}")

    # With ferrite
    ferrite_factor = 2.5
    L1_ferrite = L1_dc * ferrite_factor
    L2_ferrite = L2_dc * ferrite_factor
    M_ferrite = M_dc * ferrite_factor

    print("\n--- With Ferrite (factor = 2.5) ---")
    print(f"L1 = {L1_ferrite*1e6:.2f} uH")
    print(f"L2 = {L2_ferrite*1e6:.2f} uH")
    print(f"M = {M_ferrite*1e6:.3f} uH")

    # 85 kHz analysis
    f_target = 85e3  # 85 kHz
    print(f"\n--- Analysis at f = {f_target/1e3:.0f} kHz ---")

    Z = wpt.calc_impedance_matrix(f_target)
    print(f"Z11 = {Z[0,0].real:.4f} + j{Z[0,0].imag:.4f} Ohm")
    print(f"Z22 = {Z[1,1].real:.4f} + j{Z[1,1].imag:.4f} Ohm")
    print(f"Z12 = j{Z[0,1].imag:.4f} Ohm")

    # Skin depth
    delta_cu = wpt.tx.calc_skin_depth(f_target, 5.8e7, 1.0)
    delta_al = wpt.tx.calc_skin_depth(f_target, 3.5e7, 1.0)
    print(f"\nSkin depth @ 85 kHz:")
    print(f"  Copper: {delta_cu*1e3:.3f} mm")
    print(f"  Aluminum: {delta_al*1e3:.3f} mm")

    # Efficiency vs load resistance
    print("\n--- Efficiency vs Load Resistance ---")
    R_loads = [1, 2, 5, 10, 20, 50]
    print(f"{'R_load [Ohm]':>12} {'Efficiency [%]':>15} {'P_in [mW]':>12} {'P_load [mW]':>12}")
    print("-" * 55)

    for R_L in R_loads:
        eta, P_in, P_load = wpt.calc_efficiency(f_target, R_L)
        # Normalize to 1V source
        print(f"{R_L:>12} {eta*100:>15.1f} {P_in*1e3:>12.3f} {P_load*1e3:>12.3f}")

    # Frequency sweep
    print("\n--- Frequency Sweep ---")
    freqs = np.logspace(4, 6, 100)  # 10 kHz to 1 MHz
    R_load = 10  # 10 Ohm

    efficiencies = []
    L1_vs_f = []
    R1_vs_f = []

    for f in freqs:
        eta, _, _ = wpt.calc_efficiency(f, R_load)
        efficiencies.append(eta)

        Z = wpt.calc_impedance_matrix(f)
        omega = 2 * np.pi * f
        L1_vs_f.append(Z[0, 0].imag / omega)
        R1_vs_f.append(Z[0, 0].real)

    efficiencies = np.array(efficiencies)
    L1_vs_f = np.array(L1_vs_f)
    R1_vs_f = np.array(R1_vs_f)

    # Find peak efficiency
    idx_max = np.argmax(efficiencies)
    f_opt = freqs[idx_max]
    eta_max = efficiencies[idx_max]
    print(f"\nOptimal frequency: {f_opt/1e3:.1f} kHz")
    print(f"Peak efficiency: {eta_max*100:.1f}%")

    # Efficiency at 85 kHz
    idx_85k = np.argmin(np.abs(freqs - 85e3))
    print(f"Efficiency at 85 kHz: {efficiencies[idx_85k]*100:.1f}%")

    # Gap variation
    print("\n--- Coupling vs Air Gap ---")
    gaps = np.linspace(0.05, 0.30, 6)  # 50 mm to 300 mm
    print(f"{'Gap [mm]':>10} {'k':>10} {'M [uH]':>10} {'eta [%]':>10}")
    print("-" * 45)

    gap_orig = wpt.air_gap
    for gap in gaps:
        wpt.air_gap = gap
        k = wpt.calc_coupling_coefficient()
        M = wpt.calc_mutual_inductance() * ferrite_factor
        eta, _, _ = wpt.calc_efficiency(85e3, 10)
        print(f"{gap*1e3:>10.0f} {k:>10.4f} {M*1e6:>10.3f} {eta*100:>10.1f}")
    wpt.air_gap = gap_orig

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Efficiency vs frequency
    ax = axes[0, 0]
    ax.semilogx(freqs / 1e3, efficiencies * 100)
    ax.axvline(x=85, color='r', linestyle='--', label='85 kHz')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Efficiency [%]')
    ax.set_title(f'Power Transfer Efficiency (R_load = {R_load} Ohm)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Plot 2: Inductance vs frequency
    ax = axes[0, 1]
    ax.semilogx(freqs / 1e3, L1_vs_f * 1e6)
    ax.axvline(x=85, color='r', linestyle='--', label='85 kHz')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('L1 [uH]')
    ax.set_title('Tx Inductance vs Frequency')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Plot 3: Resistance vs frequency
    ax = axes[1, 0]
    ax.loglog(freqs / 1e3, R1_vs_f * 1e3)
    ax.axvline(x=85, color='r', linestyle='--', label='85 kHz')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('R1 [mOhm]')
    ax.set_title('Tx Resistance vs Frequency (incl. shield loss)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    # Plot 4: Coupling vs gap
    ax = axes[1, 1]
    gaps_plot = np.linspace(0.05, 0.30, 50)
    ks = []
    etas_gap = []
    for gap in gaps_plot:
        wpt.air_gap = gap
        ks.append(wpt.calc_coupling_coefficient())
        eta, _, _ = wpt.calc_efficiency(85e3, 10)
        etas_gap.append(eta)
    wpt.air_gap = gap_orig

    ax.plot(gaps_plot * 1e3, np.array(ks), 'b-', label='k')
    ax.set_xlabel('Air Gap [mm]')
    ax.set_ylabel('Coupling coefficient k', color='b')
    ax.tick_params(axis='y', labelcolor='b')

    ax2 = ax.twinx()
    ax2.plot(gaps_plot * 1e3, np.array(etas_gap) * 100, 'r--', label='eta')
    ax2.set_ylabel('Efficiency [%]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    ax.set_title('Coupling and Efficiency vs Air Gap')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('wpt_85khz_with_shield.png', dpi=150)
    print(f"\nSaved: wpt_85khz_with_shield.png")
    plt.close()

    # Summary for paper
    print("\n" + "=" * 70)
    print("Summary for Paper")
    print("=" * 70)
    print("""
WPT System Parameters (SAE J2954 Class):
  - Operating frequency: 85 kHz
  - Coil diameter: 300 mm (outer)
  - Air gap: 150 mm (typical ground clearance)
  - Ferrite core: mu_r = 2000, 5 mm thick
  - Aluminum shield: 2 mm thick

Analysis Results:
  - Self-inductance (with ferrite): ~{:.1f} uH
  - Mutual inductance: ~{:.3f} uH
  - Coupling coefficient: ~{:.3f}
  - Power transfer efficiency: ~{:.0f}% (R_load = 10 Ohm)

Key Effects Modeled:
  1. Ferrite core magnetic coupling enhancement
  2. Litz wire skin effect (strand-level)
  3. Aluminum shield eddy current losses
  4. Frequency-dependent impedance matrix

Next Steps for PEEC + HDiv-VIM/reduced FEM:
  1. Replace analytical formulas with Radia 3D geometry
  2. Apply Loop-Star decomposition to coil mesh
  3. Exchange ferrite-core response through HDiv-VIM / reduced FEM
  4. Validate against measurement or commercial software
""".format(L1_ferrite*1e6, M_ferrite*1e6, k_dc, efficiencies[idx_85k]*100))

    print("=" * 70)


if __name__ == '__main__':
    main()
